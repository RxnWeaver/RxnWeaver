package molecule

import (
	"fmt"

	cmn "github.com/RxnWeaver/RxnWeaver/common"
)

// BondBuilder builds a bond, in typical builder fashion, one property
// at a time.
//
// This type is public.  This is the only way to create a new bond
// from outside this package.  The created and constructed bond can
// NOT be extracted or used directly.  The only way to make use of it
// is to send it to a molecule for inclusion in it.
//
// An instance of builder can create and build any number of bonds.
type BondBuilder struct {
	mol *Molecule // Molecule whose bonds this builder constructs.
	b   *_Bond    // Bond being built by this builder.
}

// New creates a new bond instance in this builder.
func (bb *BondBuilder) New(id int) (*BondBuilder, error) {
	if uint16(id) != bb.mol.nextBondId {
		return nil, fmt.Errorf("Possible out-of-sequence parsing.  Expected bond ID : %d, given : %d", bb.mol.nextBondId, id)
	}

	// The molecule, in which this bond gets eventually included,
	// should set itself as the containing molecule.
	bb.b = newBond(bb.mol, id)
	return bb, nil
}

// Atoms sets the input IDs of the two atoms being bound together by
// this bond.
func (bb *BondBuilder) Atoms(aiid1, aiid2 int) (*BondBuilder, error) {
	mol := bb.mol

	a1 := mol.atomWithIid(uint16(aiid1))
	a2 := mol.atomWithIid(uint16(aiid2))
	if a1 == nil {
		return nil, fmt.Errorf("Unknown atom input ID given : %d", aiid1)
	}
	if a2 == nil {
		return nil, fmt.Errorf("Unknown atom input ID given : %d", aiid2)
	}

	// We do not add bonds to hydrogen atoms.
	if a1.atNum == 1 {
		a2.hCount++
		bb.b = nil
		return bb, fmt.Errorf("Bond involves a hydrogen atom.")
	}
	if a2.atNum == 1 {
		a1.hCount++
		bb.b = nil
		return bb, fmt.Errorf("Bond involves a hydrogen atom.")
	}

	bb.b.a1 = uint16(aiid1)
	bb.b.a2 = uint16(aiid2)
	return bb, nil
}

// BondType sets the bond order of this bond.
func (bb *BondBuilder) BondType(bType cmn.BondType) (*BondBuilder, error) {
	if bType == cmn.BondTypeNone || bType == cmn.BondTypeAltern {
		return nil, fmt.Errorf("Unhandled bond type : %v", bType)
	}

	bb.b.bType = bType
	return bb, nil
}

// BondStereo sets the stereo type of this bond.
func (bb *BondBuilder) BondStereo(bStereo cmn.BondStereo) *BondBuilder {
	bb.b.bStereo = bStereo
	return bb
}
