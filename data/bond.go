package data

import (
	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// Bond represents a chemical bond.  This is strictly between two
// atoms, and does not cater to multi-bond requirements.
//
// Bonds always relate atoms by their input IDs, NOT by their
// normalised IDs.  This is because they are constructed while reading
// the input molecule.  This is also useful when debugging, since we
// can directly correlate the bonds to the original input structure.
type Bond struct {
	id uint16 // Unique ID of this bond.

	a1      uint16         // iId of the first atom in the bond.
	a2      uint16         // iId of the second atom in the bond.
	bType   cmn.BondType   // Is this bond single, double or triple?
	bStereo cmn.BondStereo // See the enum definitions for details.

	isAro  bool   // Is this bond aromatic?
	isLink bool   // Is this bond part of a linking chain?
	hash   uint32 // For fast comparisons.

	rings [cmn.MaxRings]uint8 // The rings this bond participates in.
}

// otherAtom answers the atom other than the given one that
// participates in this bond.  Answers `0` if the given atom does not
// participate in this bond at all.
func (b *Bond) otherAtom(a uint16) uint16 {
	if b.a1 == a {
		return b.a2
	}

	if b.a2 == a {
		return b.a1
	}

	return 0
}
