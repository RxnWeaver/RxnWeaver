package data

import (
	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// Bond represents a chemical bond.  This is strictly between two
// atoms, and does not cater to multi-bond requirements.
type Bond struct {
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
// participates in this bond.  Answers (0, false) if the given atom
// does not participate in this bond at all.
func (b *Bond) otherAtom(a uint16) (uint16, bool) {
	if b.a1 == a {
		return b.a2, true
	}

	if b.a2 == a {
		return b.a1, true
	}

	return 0, false
}
