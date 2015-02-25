package molecule

import (
	"fmt"
	"math"

	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// Bond represents a chemical bond.  This is strictly between two
// atoms, and does not cater to multi-bond requirements.
//
// Bonds always relate atoms by their input IDs, NOT by their
// normalised IDs.  This is because they are constructed while reading
// the input molecule.  This is also useful when debugging, since we
// can directly correlate the bonds to the original input structure.
type _Bond struct {
	mol *Molecule // Containing molecule of this bond.
	id  uint16    // Unique ID of this bond.

	a1      uint16         // iId of the first atom in the bond.
	a2      uint16         // iId of the second atom in the bond.
	bType   cmn.BondType   // Is this bond single, double or triple?
	bStereo cmn.BondStereo // See the enum definitions for details.

	isAro  bool   // Is this bond aromatic?
	isLink bool   // Is this bond part of a linking chain?
	hash   uint32 // For fast comparisons.

	rings []uint8 // The rings this bond participates in.
}

// newBond constructs and initialises a new bond between the two given
// atoms, with the provided configuration.
func newBond(mol *Molecule, id, a1, a2 uint16, bType cmn.BondType, bStereo cmn.BondStereo) *_Bond {
	bond := new(_Bond)
	bond.mol = mol
	bond.id = id
	bond.a1 = a1
	bond.a2 = a2

	bond.bType = bType
	bond.bStereo = bStereo

	bond.rings = make([]uint8, 0, cmn.MaxRings)

	return bond
}

// otherAtomIid answers the atom other than the given one that
// participates in this bond.  Answers `0` if the given atom does not
// participate in this bond at all.
func (b *_Bond) otherAtomIid(aid uint16) uint16 {
	if b.a1 == aid {
		return b.a2
	}

	if b.a2 == aid {
		return b.a1
	}

	return 0
}

// isCyclic answers if this bond participates in at least one ring.
func (b *_Bond) isCyclic() bool {
	return len(b.rings) > 0
}

// addRing adds the given ring to the list of rings in which this bond
// participates.
func (b *_Bond) addRing(rid uint8) {
	for _, id := range b.rings {
		if id == rid {
			return
		}
	}

	b.rings = append(b.rings, rid)
}

// removeRing removes the given ring from the list of rings in which
// this bond participates, if it does participate in the given ring.
func (b *_Bond) removeRing(rid uint8) {
	idx := -1
	for i, id := range b.rings {
		if id == rid {
			idx = i
			break
		}
	}

	if idx > -1 {
		b.rings = append(b.rings[:idx], b.rings[idx+1:]...)
	}
}

// isInRing answers if this bond participates in the given ring.
func (b *_Bond) isInRing(rid uint8) bool {
	for _, id := range b.rings {
		if id == rid {
			return true
		}
	}

	return false
}

// isInRingOfSize answers if this bond participates in at least one
// ring of the given size.
func (b *_Bond) isInRingOfSize(n int) bool {
	mol := b.mol
	for _, rid := range b.rings {
		r := mol.ringWithId(rid)
		if r.size() == n {
			return true
		}
	}

	return false
}

// smallestRing answers the smallest unique ring in which this atom
// participates.  If no such unique ring exists, an error is answered.
func (b *_Bond) smallestRing() (uint8, error) {
	if !b.isCyclic() {
		return 0, fmt.Errorf("Bond is not cyclic.")
	}

	min := int(math.MaxUint8)
	c := 0
	var ret uint8

	mol := b.mol
	for _, rid := range b.rings {
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
