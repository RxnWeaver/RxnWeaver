package molecule

import (
	"fmt"
	"math"

	bits "github.com/willf/bitset"

	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// Atom represents a chemical atom.
//
// It has various physical and chemical properties, apart from the
// transient information attached to it during its participation in
// reactions.
type Atom struct {
	mol    *Molecule // Containing molecule of this atom.
	atNum  uint8     // Atomic number of this atom's element.
	symbol string    // Symbol, in case of a different isotope.
	iId    uint16    // Serial input ID of this atom.
	nId    uint16    // Normalised ID of this atom.

	X float32 // X-coordinate of this atom.
	Y float32 // Y-coordinate of this atom.
	Z float32 // Z-coordinate of this atom.

	hCount  uint8 // Number of implicit + explicit H atoms attached to this atom.
	charge  int8  // Residual net charge of this atom.
	valence int8  // Current valence configuration of this atom.

	unsaturation cmn.Unsaturation // Current composite state of this atom.

	pHash uint64 // A pseudo-hash of this atom, using some attributes.
	sHash uint64 // A pseudo-hash of this atom, using some attributes.

	bonds           *bits.BitSet // Bitmap of bonds of this atom.
	nbrs            []uint16     // Expanded list of neighbours of this atom.
	singleBondCount uint8        // Number of single bonds this atom has.
	doubleBondCount uint8        // Number of double bonds this atom has.
	tripleBondCount uint8        // Number of triple bonds this atom has.

	rings *bits.BitSet // Bitmap of IDs of rings this atom participates in.
	// Does this atom participate in at least one aromatic ring?
	isInAroRing bool
	// Is this atom a bridgehead of a bicyclic system of rings?
	isBridgeHead bool
	// Is this atom the sole common atom of all of its rings?
	isSpiro bool

	// The functional groups substituted on this atom.  They are listed in
	// descending order of importance.  The first is the primary feature.
	features []uint16

	// Number of electron-donating neighbours.
	edNbrCount int
	// Number of unsaturated electron-withdrawing neighbours.
	unsatEwNbrCount int
	// Number of saturated electron-withdrawing neighbours.
	satEwNbrCount int
}

// newAtom constructs and initialises a new atom of the given element
// type, and belonging to the given molecule.
func newAtom(mol *Molecule, atNum uint8) *Atom {
	atom := new(Atom)
	atom.mol = mol
	atom.atNum = atNum
	atom.valence = cmn.PeriodicTable[cmn.ElementSymbols[atNum]].Valence

	atom.bonds = bits.New(cmn.MaxBonds)
	atom.nbrs = make([]uint16, 0, cmn.MaxBonds)
	atom.rings = bits.New(cmn.MaxRings)

	atom.features = make([]uint16, 0, cmn.MaxFeatures)

	return atom
}

// AtomicNumber answers the atomic number of this atom.
func (a *Atom) AtomicNumber() uint8 {
	return a.atNum
}

// Parent answers the parent molecule of this atom.
func (a *Atom) Parent() *Molecule {
	return a.mol
}

// determineUnsaturation computes a composite metric that reflects the
// current state of the atom.
//
// This method is expected to be invoked during a molecule's
// normalisation only.
//
// Note that this method underlies - directly or indirectly - major
// parts of RxnWeaver's decision rules.  Exercise great caution should
// you need to modify this in any manner!
func (a *Atom) determineUnsaturation() error {
	nb := int(a.bonds.Count())
	nn := len(a.nbrs)

	// Atom has a residual charge.
	if a.charge != 0 {
		a.unsaturation = cmn.UnsaturationCharged
		return nil
	}

	// For an uncharged atom, valence should be sane.
	if a.hCount > 0 {
		os := int8(nn) + int8(a.hCount)
		if ok, err := cmn.IsValidOxidationState(a.atNum, os); !ok {
			return err
		}
	}

	// Case of all single bonds.
	if nb == nn {
		a.unsaturation = cmn.UnsaturationNone
		return nil
	}

	// Double or triple bonds exist.
	ndb := 0
	nhdb := 0
	ntb := 0
	nhtb := 0
	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid - 1))
		oaid := b.otherAtomIid(a.iId)
		oa := mol.atomWithIid(oaid)
		switch b.bType {
		case cmn.BondTypeDouble:
			ndb++
			if oa.atNum != 6 {
				nhdb++
			}
		case cmn.BondTypeTriple:
			ntb++
			if oa.atNum != 6 {
				nhtb++
			}
		}
	}

	if ntb > 0 {
		if nhtb > 0 {
			a.unsaturation = cmn.UnsaturationTripleBondW
		} else {
			a.unsaturation = cmn.UnsaturationTripleBondC
		}
		// TODO(js): Assess the exhaustiveness.
		return nil
	}

	if ndb > 0 {
		switch {
		case ndb == 1 && nhdb == 0:
			a.unsaturation = cmn.UnsaturationDoubleBondC
		case ndb == 1 && nhdb == 1:
			a.unsaturation = cmn.UnsaturationDoubleBondW
		case ndb == 2 && nhdb == 0:
			a.unsaturation = cmn.UnsaturationDoubleBondCC
		case ndb == 2 && nhdb == 1:
			a.unsaturation = cmn.UnsaturationDoubleBondCW
		case ndb == 2 && nhdb == 2:
			a.unsaturation = cmn.UnsaturationDoubleBondWW
		}
		// TODO(js): Assess the exhaustiveness.
	}
	return nil
}

// piElectronCount answers the number of delocalised pi electrons
// contributed by this atom.  This number is important for calculating
// the aromaticity of the rings this atom participates in.
func (a *Atom) piElectronCount() int {
	mol := a.mol
	wtSum := 100*int16(a.doubleBondCount) + 10*int16(a.singleBondCount) + int16(a.charge)

	switch a.atNum {
	case 6:
		switch wtSum {
		case 19:
			return 2
		case 110:
			return 1
		case 120:
			var b *Bond
			for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
				b = mol.bondWithId(uint16(bid - 1))
				if b.bType == cmn.BondTypeDouble {
					break
				}
			}
			if b.isCyclic() {
				return 1
			}
			return 0
		default:
			return 0
		}

	case 7:
		switch wtSum {
		case 20:
			return 2
		case 30:
			return 2
		case 110:
			return 1
		case 121:
			return 1
		default:
			return 0
		}

	case 8:
		switch wtSum {
		case 20:
			return 2
		case 111:
			return 1
		default:
			return 0
		}

	case 16:
		switch wtSum {
		case 20:
			return 2
		case 111:
			return 1
		case 120:
			var b *Bond
			for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
				b = mol.bondWithId(uint16(bid - 1))
				if b.bType == cmn.BondTypeDouble {
					break
				}
			}
			oaid := b.otherAtomIid(a.iId)
			oa := mol.atomWithIid(oaid)
			if oa.atNum == 8 && !oa.isCyclic() {
				return 2
			}
			return 0
		case 220:
			c := 0
			for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
				b := mol.bondWithId(uint16(bid - 1))
				if b.bType == cmn.BondTypeDouble {
					oaid := b.otherAtomIid(a.iId)
					if !mol.atomWithIid(oaid).isCyclic() {
						c++
					}
				}
			}
			if c > 1 {
				return -1
			}
			return 0
		default:
			return 0
		}
	}

	return 0
}

// isCyclic answers if this atom participates in at least one ring.
func (a *Atom) isCyclic() bool {
	return a.rings.Count() > 0
}

// isJunction answers if this atom has more than 2 distinct
// neighbours.
func (a *Atom) isJunction() bool {
	return a.bonds.Count() > 2
}

// addBond adds the given bond to this atom, if it is not already
// present.  It also adjusts the list of its neighbours appropriately.
//
// Note that it does NOT check to see if the addition conforms to this
// atom's current valence configuration.
func (a *Atom) addBond(b *Bond) {
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		if uint16(bid) == b.id {
			return
		}
	}

	a.bonds.Set(uint(b.id))
	nbrId := b.otherAtomIid(a.iId)
	n := int(b.bType)
	for i := 0; i < n; i++ {
		a.nbrs = append(a.nbrs, nbrId)
	}

	switch n {
	case 1:
		a.singleBondCount++
	case 2:
		a.doubleBondCount++
	case 3:
		a.tripleBondCount++
	}
}

// removeBond removes the given bond from this atom, and adjusts the
// list of neighbours appropriately.
//
// Note that it does NOT check to see if the removal conforms to this
// atom's current valence configuration.
func (a *Atom) removeBond(b *Bond) {
	nbrId := b.otherAtomIid(a.iId)

	switch b.bType {
	case 1:
		a.singleBondCount--
	case 2:
		a.doubleBondCount--
	case 3:
		a.tripleBondCount--
	}

	wid := 0
	for _, nid := range a.nbrs {
		if nid == nbrId {
			continue
		}
		a.nbrs[wid] = nid
		wid++
	}
	a.nbrs = a.nbrs[:wid]

	a.bonds.Clear(uint(b.id))
}

// bondTo answers the bond that binds this atom to the given atom, if
// one such bond exists.  Answers `nil` otherwise.
func (a *Atom) bondTo(other uint16) *Bond {
	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid - 1))
		if b.otherAtomIid(a.iId) == other {
			return b
		}
	}

	return nil
}

// firstDoublyBondedNbr answers this atom's doubly-bonded neighbour
// having the highest priority.
//
// This method assumes that the molecule is already normalised!
// Calling it on a molecule that has not be normalised yet, leads to
// incorrect results.
func (a *Atom) firstDoublyBondedNbr() uint16 {
	if a.doubleBondCount == 0 {
		return 0
	}

	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid - 1))
		if b.bType == cmn.BondTypeDouble {
			return b.otherAtomIid(a.iId)
		}
	}

	panic("Should never be here!")
}

// firstMultiplyBondedNbr answers this atom's doubly-bonded neighbour
// having the highest priority.
//
// This method assumes that the molecule is already normalised!
// Calling it on a molecule that has not be normalised yet, leads to
// incorrect results.
func (a *Atom) firstMultiplyBondedNbr() uint16 {
	if a.doubleBondCount == 0 && a.tripleBondCount == 0 {
		return 0
	}

	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid - 1))
		if b.bType >= cmn.BondTypeDouble {
			return b.otherAtomIid(a.iId)
		}
	}

	panic("Should never be here!")
}

// isInRingOfSize answers if this atom participates in at least one
// ring of the given size.
func (a *Atom) isInRingOfSize(n int) bool {
	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		if r.size() == n {
			return true
		}
	}

	return false
}

// isInRingLargerThan answers if this atom participates in at least
// one ring that is larger than the given number.
func (a *Atom) isInRingLargerThan(n int) bool {
	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		if r.size() > n {
			return true
		}
	}

	return false
}

// smallestRing answers the smallest unique ring in which this atom
// participates.  If no such unique ring exists, an error is answered.
func (a *Atom) smallestRing() (uint8, error) {
	if !a.isCyclic() {
		return 0, fmt.Errorf("Atom not cyclic.")
	}

	min := int(math.MaxUint8)
	c := 0
	var ret uint8

	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		size := r.size()
		if size == min {
			c++
		} else if size < min {
			ret = uint8(rid)
			c = 1
		}
	}

	if c > 1 {
		return 0, fmt.Errorf("Smallest ring size: %d, number of smallest rings: %d.", min, c)
	}

	return ret, nil
}

// isAromatic answers if this atom is part of an aromatic ring.
//
// Note that the actual aromaticity determination is handled by
// `Ring`.  This method merely answers the set flag.
func (a *Atom) isAromatic() bool {
	return a.isInAroRing
}

// isInHeteroAromaticRing answers if this atom is part of an aromatic
// ring with at least one hetero atom.
func (a *Atom) isInHeteroAromaticRing() bool {
	if a.isInAroRing && a.atNum != 6 {
		return true // Simplest case!
	}

	mol := a.mol
	for rid, ok := a.rings.NextSet(0); ok; rid, ok = a.rings.NextSet(rid + 1) {
		r := mol.ringWithId(uint8(rid))
		if r.isHeteroAromatic() {
			return true
		}
	}

	return false
}

// haveCommonRings answers if this atom has at least one ring in
// common with the given atom.
func (a *Atom) haveCommonRings(aiid uint16) bool {
	other := a.mol.atomWithIid(aiid)

	return a.rings.IntersectionCardinality(other.rings) > 0
}

// isInSameRingsAs answers if this atom participates in exactly the
// same rings as the given atom.
func (a *Atom) isInSameRingsAs(aiid uint16) bool {
	other := a.mol.atomWithIid(aiid)

	return a.rings.Equal(other.rings)
}

// isInAllRingsOf answers if this atom participates in every ring in
// which the given atom does.
//
// Note that this atom may participate in more rings, as well.
func (a *Atom) isInAllRingsOf(aiid uint16) bool {
	other := a.mol.atomWithIid(aiid)

	return other.rings.DifferenceCardinality(a.rings) == 0
}

// addRing adds the given ring to the list of this atom's rings.
func (a *Atom) addRing(r *Ring) {
	a.rings.Set(uint(r.id))
}

// removeRing removes the given ring from the list of this atom's
// rings.
//
// This method should, usually, not be called.  The composition of a
// ring never changes once it is completed.  However, during a
// reaction, a ring may get broken, leading to its death.  The
// constituent atoms are then notified of that death by calling this
// method.
func (a *Atom) removeRing(r *Ring) {
	a.rings.Clear(uint(r.id))
}

// functionalGroup answers the primary feature of this atom, if one is
// present.  Answers `0` otherwise.
func (a *Atom) functionalGroup() uint16 {
	if len(a.features) == 0 {
		return 0
	}

	return a.features[0]
}

// addFeature adds the given feature to this atom's list of features.
func (a *Atom) addFeature(fid uint16) {
	a.features = append(a.features, fid)
}

// removeFeature removes the first instance of the given feature from
// this atom's list of features, if it exists in it.
//
// Answers `true` upon a successful removal; `false` otherwise.
func (a *Atom) removeFeature(fid uint16) bool {
	idx := -1
	for i, f := range a.features {
		if f == fid {
			idx = i
			break
		}
	}

	if idx == -1 {
		return false
	}

	a.features = append(a.features[:idx], a.features[idx+1:]...)
	return true
}

// featureCount answers the number of features substituted on this
// atom.
func (a *Atom) featureCount() int {
	return len(a.features)
}

// hasFeature answers if this atom has at least one substituent
// matching the given group.
func (a *Atom) hasFeature(fid uint16) bool {
	for _, f := range a.features {
		if f == fid {
			return true
		}
	}

	return false
}

// isFunctional answers if this atom can play an active role in a
// reaction, in a substituting position.
//
// Note that an atom can yet be a reaction centre without being
// functional.
func (a *Atom) isFunctional() bool {
	if a.atNum != 6 {
		return true
	}
	if len(a.features) > 0 {
		return true
	}
	if a.unsaturation > cmn.UnsaturationNone {
		return true
	}

	return false
}

// electronWithdrawingNeighbourCount answers the total of unsaturated
// and saturated electron-withdrawing neighbours of this atom.
func (a *Atom) electronWithdrawingNeighbourCount() int {
	return a.unsatEwNbrCount + a.satEwNbrCount
}

// enolicHydrogenCount answers the number of attached hydrogen atoms,
// if this atom is enolic.
//
// For the attached hydrogen atoms to be enolic, this atom must be a
// carbon, must be saturated, must have at least one
// electron-withdrawing neighbour, and must not be a bridgehead.
func (a *Atom) enolicHydrogenCount() int {
	if a.atNum != 6 {
		return 0
	}
	if a.unsaturation > cmn.UnsaturationNone {
		return 0
	}
	if a.electronWithdrawingNeighbourCount() == 0 {
		return 0
	}
	if a.isBridgeHead {
		return 0
	}

	return int(a.hCount)
}

// isAtomicLeavingGroup answers if this atom can act as a leaving
// group.
func (a *Atom) isAtomicLeavingGroup() bool {
	return a.bonds.Count() == 1 && a.atNum != 6
}

// isCH2 answers if this atom is a carbon with exactly two hydrogen
// atoms bound to it.
func (a *Atom) isCH2() bool {
	return a.atNum == 6 && a.hCount == 2
}

// isCH3 answers if this atom is a carbon with exactly three hydrogen
// atoms bound to it.
func (a *Atom) isCH3() bool {
	return a.atNum == 6 && a.hCount == 3
}

// isCarbonylC answers if this atom is a carbon with exactly one
// double bond with an oxygen atom.
func (a *Atom) isCarbonylC() bool {
	if a.atNum != 6 {
		return false
	}
	if a.unsaturation != cmn.UnsaturationDoubleBondW {
		return false
	}

	mol := a.mol
	for bid, ok := a.bonds.NextSet(0); ok; bid, ok = a.bonds.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid - 1))
		if b.bType == cmn.BondTypeDouble {
			oaid := b.otherAtomIid(a.iId)
			oa := mol.atomWithIid(oaid)
			if oa.atNum == 8 {
				return true
			}
		}
	}

	return false
}

// isHydroxyl answers if the current atom is an oxygen with exactly
// one hydrogen atom bound to it.
func (a *Atom) isHydroxyl() bool {
	return a.atNum == 8 && a.hCount == 1
}
