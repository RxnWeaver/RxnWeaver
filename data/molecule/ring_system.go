package molecule

import (
	"fmt"

	bits "github.com/willf/bitset"

	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// RingSystem represents a set of physically fused rings.  Fusion is
// of one of the following types.
//
//   - One bond in common between a given pair of rings.
//   - One atom in common between a given pair of rings (spiro).
//   - More than one bond in common between a given pair of rings.
//
// The last type involves bridge atoms and bridgehead atoms.
//
// A ring system is mutable, unlike a ring.  Its composition can
// change during the course of the life of its molecule.
type _RingSystem struct {
	mol *Molecule // Containing molecule of this ring system.
	id  uint8     // Unique ID of this ring system in its molecule.

	rings      []uint8      // List of rings comprising this system.
	atomBitSet *bits.BitSet // All atoms from all rings in this system.
	bondBitSet *bits.BitSet // All bonds from all rings in this system.

	isAro bool // Is this ring system aromatic as a whole?
}

// newRingSystem creates and initialises a ring system with the given
// molecule and unique ID.
func newRingSystem(mol *Molecule, id uint8) *_RingSystem {
	rs := new(_RingSystem)
	rs.mol = mol
	rs.id = id

	rs.rings = make([]uint8, 0, cmn.MaxRings)
	rs.atomBitSet = bits.New(cmn.ListSizeSmall)
	rs.bondBitSet = bits.New(cmn.ListSizeSmall)

	return rs
}

// size answers the number of rings in this system.
func (rs *_RingSystem) size() int {
	return len(rs.rings)
}

// hasRing answers if this system includes the given ring.
func (rs *_RingSystem) hasRing(rid uint8) (bool, int) {
	for i, id := range rs.rings {
		if id == rid {
			return true, i
		}
	}

	return false, -1
}

// ringAt answers the ring at the given index.
//
// This method does not perform any index boundary checks!
func (rs *_RingSystem) ringAt(idx int) uint8 {
	return rs.rings[idx]
}

// addRing adds the given ring to this system, if it does not already
// belong to it.
//
// It also updates the internal bitsets of atoms and bonds
// appropriately.
func (rs *_RingSystem) addRing(r *_Ring) error {
	return rs.addRingAt(rs.size(), r)
}

// addRingAt add the given ring to this system, at the given index, if
// it does not already belong to it.
//
// This method does not perform any index boundary checks!
//
// It also updates the internal bitsets of atoms and bonds
// appropriately.
//
//  This method is idempotent.
func (rs *_RingSystem) addRingAt(idx int, r *_Ring) error {
	for _, rid := range rs.rings {
		if rid == r.id {
			return nil
		}
	}

	// Given ring should have at least one bond or atom in common with
	// others in this system.
	if rs.bondBitSet.Count() > 0 {
		if rs.bondBitSet.IntersectionCardinality(r.bondBitSet) == 0 {
			if rs.atomBitSet.IntersectionCardinality(r.atomBitSet) == 0 {
				return fmt.Errorf("Ring %d has no bonds or atoms in common with any others in this ring system")
			}
		}
	}

	t := rs.rings[idx:]
	rs.rings = append(rs.rings[:idx], r.id)
	rs.rings = append(rs.rings, t...)

	rs.atomBitSet.InPlaceUnion(r.atomBitSet)
	rs.bondBitSet.InPlaceUnion(r.bondBitSet)

	return nil
}

// removeRing removes the given ring from this system.
//
// It also updates the internal bitsets of atoms and bonds
// appropriately.
//
//  This method is idempotent.
func (rs *_RingSystem) removeRing(r *_Ring) error {
	idx := -1
	for i, rid := range rs.rings {
		if rid == r.id {
			idx = i
			break
		}
	}
	if idx == -1 {
		return nil
	}

	return rs.removeRingAt(idx)
}

// removeRingAt removes the ring at the given index from this system.
//
// It also updates the internal bitsets of atoms and bonds
// appropriately.
//
//  This method is idempotent.
func (rs *_RingSystem) removeRingAt(idx int) error {
	rs.rings = append(rs.rings[:idx], rs.rings[idx+1:]...)

	// Rebuild the bitsets.
	rs.atomBitSet.ClearAll()
	rs.bondBitSet.ClearAll()

	mol := rs.mol
	for _, rid := range rs.rings {
		r := mol.ringWithId(rid)
		rs.atomBitSet.InPlaceUnion(r.atomBitSet)
		rs.bondBitSet.InPlaceUnion(r.bondBitSet)
	}

	return nil
}

// piElectronCount answers the total number of delocalised pi
// electrons in this ring system.
func (rs *_RingSystem) piElectronCount() (int, bool) {
	n := 0
	mol := rs.mol
	abs := rs.atomBitSet
	for aiid, ok := abs.NextSet(0); ok; aiid, ok = abs.NextSet(aiid + 1) {
		a := mol.atomWithIid(uint16(aiid))
		if c, ok := a.piElectronCount(); ok {
			n += c
		} else {
			return 0, false
		}
	}
	return n, true
}

// determineAromaticity answers if this ring system, when considered
// as a whole, behaves like an aromatic ring.
//
// If the system is aromatic, its constituent rings are not tested
// individually for aromaticity.  This could change in future,
// depending on exceptions.
func (rs *_RingSystem) determineAromaticity() {
	err := false

	n, ok := rs.piElectronCount()
	if !ok { // Some condition preventing this ring from becoming aromatic.
		err = true
	}

	mol := rs.mol

	// First, we apply Huckel's rule.
	if (n-2)%4 != 0 {
		// TODO(js): Take exceptions into account.
		err = true
	}

	abs := rs.atomBitSet
	for aiid, ok := abs.NextSet(0); ok; aiid, ok = abs.NextSet(aiid + 1) {
		a := mol.atomWithIid(uint16(aiid))
		if a.atNum == 6 {
			if a.unsaturation == cmn.UnsaturationNone {
				err = true // There should not be any sp3 C atoms.
				break
			}
		}
	}

	// TODO(js): Take exceptions into account.

	if !err {
		// If we have come this far, this is an aromatic ring system.
		rs.isAro = true
		rs.markAtomsBondsAromatic()
	} else {
		for _, rid := range rs.rings {
			r := mol.ringWithId(rid)
			r.determineAromaticity()
		}
	}
}

// markAtomsBondsAromatic marks all participating atoms and bonds are
// being aromatic.
func (rs *_RingSystem) markAtomsBondsAromatic() {
	mol := rs.mol

	abs := rs.atomBitSet
	for aiid, ok := abs.NextSet(0); ok; aiid, ok = abs.NextSet(aiid + 1) {
		a := mol.atomWithIid(uint16(aiid))
		a.isInAroRing = true
	}

	bbs := rs.bondBitSet
	for bid, ok := bbs.NextSet(0); ok; bid, ok = bbs.NextSet(bid + 1) {
		b := mol.bondWithId(uint16(bid))
		b.isAro = true
	}
}
