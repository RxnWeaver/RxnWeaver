package data

import (
	"fmt"

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
type Ring struct {
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
func newRing(mol *Molecule) {
	r := new(Ring)
	r.mol = mol

	r.atoms = make([]uint16, 0, cmn.ListSizeSmall)
	r.bonds = make([]uint16, 0, cmn.ListSizeSmall)
	r.nbrs = make([]uint8, 0, cmn.ListSizeSmall)
}

// size answers the size of this ring.  It is equivalently the number
// of atoms or the number of bonds participating in this ring.
func (r *Ring) size() int {
	return len(r.atoms)
}

// hasAtom answers `true` if this ring includes the given atom.
// Answers `false` otherwise.
func (r *Ring) hasAtom(aid uint16) bool {
	return r.atomBitSet.Test(uint(aid))
}

// atomIndex answers the index of the given atom in this ring, if it
// is found.  Answers `-1` otherwise.
//
// Note that the answer may or may not be idempotent, depending on the
// normalisation status of the ring.  Should the ring be normalised
// between two invocations of this method, the answers could vary.
func (r *Ring) atomIndex(aid uint16) int {
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

// hasBond answers `true` if this ring includes the given bond.
// Answers `false` otherwise.
func (r *Ring) hasBond(bid uint16) bool {
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
func (r *Ring) addAtom(aid uint16) error {
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
	return nil
}

// complete closes the link between the last atom in the ring and the
// first.  This operation effectively freezes the ring.
//
// This method is idempotent.
func (r *Ring) complete() error {
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

	for _, aid := range r.atoms {
		r.atomBitSet.Set(uint(aid))
	}
	for _, bid := range r.bonds {
		r.bondBitSet.Set(uint(bid))
	}

	r.isComplete = true
	return nil
}
