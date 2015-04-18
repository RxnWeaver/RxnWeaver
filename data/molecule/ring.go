package molecule

import (
	"fmt"
	"math"

	bits "github.com/willf/bitset"

	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// Ring represents a simple cycle in a molecule.
//
// A ring holds information of the atoms and the bonds it comprises.
// It also knows its neighbouring rings.
//
// Rings are supposed to be immutable: once completed, their
// composition should never change.
//
// The atom IDs held by rings are their input IDs, to match those held
// by bonds.  That makes using atoms and bonds together easier, when
// doing ring detection, etc.
type _Ring struct {
	mol  *Molecule // Containing molecule of this ring.
	id   uint8     // A unique identifier for this ring.
	rsId uint8     // ID of the ring system to which this ring belongs.

	atoms []uint16 // List of atoms participating in this ring.
	bonds []uint16 // List of bonds participating in this ring.
	nbrs  []uint8  // List of rings neighbouring this ring.

	atomBitSet *bits.BitSet // For faster comparison.
	bondBitSet *bits.BitSet // For faster comparison.

	isAro    bool // Is this ring aromatic?
	isHetAro bool // Is this an aromatic ring with at least one hetero atom?

	isComplete bool // Has this ring been finalised?
}

// newRing creates and initialises a new ring.
func newRing(mol *Molecule, id uint8) {
	r := new(_Ring)
	r.mol = mol
	r.id = id

	r.atoms = make([]uint16, 0, cmn.ListSizeSmall)
	r.bonds = make([]uint16, 0, cmn.ListSizeSmall)
	r.nbrs = make([]uint8, 0, cmn.ListSizeSmall)

	r.atomBitSet = bits.New(cmn.ListSizeSmall)
	r.bondBitSet = bits.New(cmn.ListSizeSmall)
}

// size answers the size of this ring.  It is equivalently the number
// of atoms or the number of bonds participating in this ring.
func (r *_Ring) size() int {
	return len(r.atoms)
}

// hasAtom answers if this ring includes the given atom.
func (r *_Ring) hasAtom(aid uint16) bool {
	return r.atomBitSet.Test(uint(aid))
}

// atomIndex answers the index of the given atom in this ring, if it
// is found.  Answers `-1` otherwise.
//
// Note that the answer may or may not be idempotent, depending on the
// normalisation status of the ring.  Should the ring be normalised
// between two invocations of this method, the answers could vary.
func (r *_Ring) atomIndex(aid uint16) int {
	if !r.hasAtom(aid) {
		return -1
	}

	for i, id := range r.atoms {
		if id == aid {
			return i
		}
	}

	panic("Should never be here!")
}

// hasBond answers if this ring includes the given bond.
func (r *_Ring) hasBond(bid uint16) bool {
	return r.bondBitSet.Test(uint(bid))
}

// addAtom adds the given atom to this ring.
//
// This method errors if the ring is already completed.
//
// It checks to see that a bond exists between the most-recently-added
// atom, if one such exists, and the given atom.  Answers a non-nil
// error otherwise.
//
// This method is idempotent: the given atom is ignored if it is
// already a member of this ring.
func (r *_Ring) addAtom(aid uint16) error {
	if r.isComplete {
		return fmt.Errorf("Ring already complete.  ID : %d.", r.id)
	}

	if r.hasAtom(aid) {
		return nil
	}

	size := len(r.atoms)
	if size == 0 {
		r.atoms = append(r.atoms, aid)
		return nil
	}

	prev := r.atoms[size-1]
	b := r.mol.bondBetween(prev, aid)
	if b == nil {
		return fmt.Errorf("No bond exists between atom %d and atom %d.", prev, aid)
	}

	r.bonds = append(r.bonds, b.id)
	r.atoms = append(r.atoms, aid)
	r.bondBitSet.Set(uint(b.id))
	r.atomBitSet.Set(uint(aid))
	return nil
}

// complete closes the link between the last atom in the ring and the
// first.  This operation effectively freezes the ring.
//
// This method is idempotent.
func (r *_Ring) complete() error {
	if r.isComplete {
		return nil
	}

	size := len(r.atoms)
	if size < 3 {
		return fmt.Errorf("A ring must have at least 3 atoms.  This ring has only %d.", size)
	}

	aid1 := r.atoms[0]
	aid2 := r.atoms[size-1]
	b := r.mol.bondBetween(aid1, aid2)
	if b == nil {
		return fmt.Errorf("No bond between first atom %d and last atom %d.", aid1, aid2)
	}

	r.bonds = append(r.bonds, b.id)

	r.isComplete = true
	return nil
}

// isAromatic answers if this ring is aromatic.
//
// The actual aromaticity determination happens when
// `determineAromaticity` is called.  This method merely answers the
// set flag.
func (r *_Ring) isAromatic() bool {
	return r.isAro
}

// isHeteroAromatic answers if this ring is aromatic with at least one
// hetero atom.
//
// The actual aromaticity determination happens when
// `determineAromaticity` is called.  This method merely answers the
// set flag.
func (r *_Ring) isHeteroAromatic() bool {
	return r.isHetAro
}

// normalise transforms the ring into a `standard' representation, in
// which the ring logically begins with that atom which has the lowest
// normalised ID.
func (r *_Ring) normalise() error {
	l := len(r.atoms)
	if l == 0 {
		return fmt.Errorf("Cannot normalise an empty ring!")
	}

	nids := make([]uint16, l, l)

	mol := r.mol
	for i, aiid := range r.atoms {
		a := mol.atomWithIid(aiid)
		nids[i] = a.nId
	}

	min := uint16(math.MaxUint16)
	idx := -1
	for i, nid := range nids {
		if nid < min {
			idx = i
			min = nid
		}
	}

	// Rotate the ring so that the atom at `idx` becomes the first.
	r.atoms = append(r.atoms[idx:], r.atoms[:idx]...)
	return nil
}

// piElectronCount answers the total number of pi-electrons in this
// ring.
func (r *_Ring) piElectronCount() (int, bool) {
	n := 0
	mol := r.mol
	for _, aiid := range r.atoms {
		a := mol.atomWithIid(aiid)
		if c, ok := a.piElectronCount(); ok {
			n += c
		} else {
			return 0, false
		}
	}
	return n, true
}

// determineAromaticity examines this ring to see if it is aromatic in
// nature.
//
// TODO(js): May have to take exceptions into account, as we make
// progress.
func (r *_Ring) determineAromaticity() {
	n, ok := r.piElectronCount()
	if !ok { // Some condition preventing this ring from becoming aromatic.
		return
	}

	mol := r.mol

	// First, we apply Huckel's rule.
	if (n-2)%4 != 0 {
		// TODO(js): Take exceptions into account.
		return
	}

	for _, aiid := range r.atoms {
		a := mol.atomWithIid(aiid)
		if a.atNum == 6 {
			if a.unsaturation == cmn.UnsaturationNone {
				return // There should not be any sp3 C atoms.
			}
		}
	}

	// TODO(js): Take exceptions into account.

	// If we have come this far, this is an aromatic ring.
	r.isAro = true

	for _, aiid := range r.atoms {
		a := mol.atomWithIid(aiid)
		a.isInAroRing = true
		if a.atNum != 6 {
			r.isHetAro = true
		}
	}
	for _, bid := range r.bonds {
		b := mol.bondWithId(bid)
		b.isAro = true
	}
}

// commonAtoms answers a list of the atoms that participate in both
// this ring and the given ring.  The representation is a bitset.
func (r *_Ring) commonAtoms(other *_Ring) *bits.BitSet {
	return r.atomBitSet.Intersection(other.atomBitSet)
}

// commonBonds answers a list of the bonds that participate in both
// this ring and the given ring.  The representation is a bitset.
func (r *_Ring) commonBonds(other *_Ring) *bits.BitSet {
	return r.bondBitSet.Intersection(other.bondBitSet)
}

// distanceBetweenAtoms answers the shorter distance in the ring,
// between the two given atoms.
func (r *_Ring) distanceBetweenAtoms(aid1, aid2 uint16) (int, error) {
	if !r.hasAtom(aid1) {
		return 0, fmt.Errorf("Atom %d is not a member of this ring.", aid1)
	}
	if !r.hasAtom(aid2) {
		return 0, fmt.Errorf("Atom %d is not a member of this ring.", aid2)
	}

	i1, i2 := -1, -1
	c := 0
	for i, aid := range r.atoms {
		switch {
		case aid == aid1:
			i1 = i
			c++
		case aid == aid2:
			i2 = i
			c++
		}
		if c == 2 {
			break
		}
	}

	d1 := i1 - i2
	if d1 < 0 {
		d1 = -d1
	}
	d2 := r.size() - d1
	if d1 < d2 {
		return d1, nil
	}
	return d2, nil
}

// isSemiAromaticOfSize6 answers if this ring satisfies the following
// constraints:
//
//     number of aromatic atoms + number of atoms in double bonds +
//         number of carbons with exocyclic double bonds to hetero atoms +
//         number of NH nitrogens
//     == 6
//
//     AND
//
//     number of carbons with exocyclic double bonds to hetero atoms
//     == number of NH nitrogens.
func (r *_Ring) isSemiAromaticOfSize6() bool {
	if r.size() != 6 || r.isAro {
		return false
	}

	nAro := r.aromaticAtomCount()
	nDbly := r.doubleBondCount() * 2

	nNH := 0
	nExo := 0
	mol := r.mol
	for _, aiid := range r.atoms {
		a := mol.atomWithIid(aiid)

		switch a.atNum {
		case 6:
			for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
				b := mol.bondWithId(uint16(bid))
				if !b.isCyclic() && b.bType == cmn.BondTypeDouble {
					nbrId := b.otherAtomIid(aiid)
					nbr := mol.atomWithIid(nbrId)
					if nbr.atNum != 6 { // Exocyclic hetero neighbour.
						nExo++
						break // Carbon can have only one such.
					}
				}
			}

		case 7:
			if a.hCount == 1 {
				nNH++
			}
		}
	}

	// All constituents are ready.
	sum := nAro + nDbly + nNH + nExo
	return sum == 6 && nNH == nExo
}

// aromaticAtomCount answers the number of atoms in this ring that are
// marked as aromatic.
//
// Note that it is possible for a non-aromatic ring to contain some
// aromatic atoms.
func (r *_Ring) aromaticAtomCount() int {
	c := 0
	mol := r.mol
	for _, aiid := range r.atoms {
		a := mol.atomWithIid(aiid)
		if a.isInAroRing {
			c++
		}
	}
	return c
}

// doubleBondCount answers the number of double bonds in this ring.
func (r *_Ring) doubleBondCount() int {
	c := 0
	mol := r.mol
	for _, bid := range r.bonds {
		b := mol.bondWithId(bid)
		if b.bType == cmn.BondTypeDouble {
			c++
		}
	}
	return c
}

// hasAdjacentAtomsSatisfying answers if this ring has two adjacent
// atoms, both of which satisfy the given constraint.
//
// Upon success, it also answers the index of the first of the two
// atoms.
func (r *_Ring) hasAdjacentAtomsSatisfying(f func(*_Atom) bool) (bool, int) {
	found := false
	mol := r.mol
	for i, aiid := range r.atoms {
		if a := mol.atomWithIid(aiid); f(a) {
			if found {
				return true, i - 1
			}
			found = true
		} else {
			found = false
		}
	}

	// If the most-recently-found carbonyl C is the last atom in the
	// ring, we wrap around, and check the first atom again.
	if found {
		if a := mol.atomWithIid(r.atoms[0]); f(a) {
			return true, len(r.atoms) - 1
		}
	}

	return false, -1
}

// hasAdjacentCarbonyls answers if this ring has two adjacent atoms,
// both of which are carbonyl carbons.
func (r *_Ring) hasAdjacentCarbonyls() (bool, int) {
	return r.hasAdjacentAtomsSatisfying(func(a *_Atom) bool {
		return a.isCarbonylC()
	})
}

// hasAdjacentSaturatedCC answers if this ring has two adjacent atoms,
// both of which are saturated carbons.
func (r *_Ring) hasAdjacentSaturatedCC() (bool, int) {
	return r.hasAdjacentAtomsSatisfying(func(a *_Atom) bool {
		return a.isSaturatedC()
	})
}

// hasAdjacentCHCH answers if this ring has two adjacent atoms, both
// of which are saturated carbons with at least one hydrogen atom each
// bound to them.
func (r *_Ring) hasAdjacentCHCH() (bool, int) {
	return r.hasAdjacentAtomsSatisfying(func(a *_Atom) bool {
		return a.isSaturatedC() && a.hCount > 0
	})
}
