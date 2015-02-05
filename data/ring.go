package data

import (
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
