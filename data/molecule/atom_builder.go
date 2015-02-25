package molecule

import (
	"fmt"

	cmn "github.com/RxnWeaver/rxnweaver/common"
)

// AtomBuilder builds an atom, in typical builder fashion, one
// property at a time.
//
// This type is public.  This is the only way to create a new atom
// from outside this package.  The created and constructed atom can
// NOT be extracted or used directly.  The only way to make use of it
// is to send it to a molecule for inclusion in it.
//
// An instance of builder can create and build only one atom.
type AtomBuilder struct {
	a *_Atom // Atom being built by this builder.
}

// NewAtomBuilder answers a new atom builder.
func NewAtomBuilder() *AtomBuilder {
	return &AtomBuilder{}
}

// New creates a new atom instance in this builder.
func (ab *AtomBuilder) New(iId int, sym string) (*AtomBuilder, error) {
	if iId <= 0 {
		return nil, fmt.Errorf("Input ID of an atom must be positive; given : %d", iId)
	}

	el, ok := cmn.PeriodicTable[sym]
	if !ok {
		return nil, fmt.Errorf("Unknown symbol : %s", sym)
	}

	// The molecule, in which this atom gets eventually included,
	// should set itself as the containing molecule.
	ab.a = newAtom(nil, el.Number, iId)
	return ab, nil
}

// Coordinates sets the given coordinates as the X-, Y- and
// Z-coordinates of this atom.
func (ab *AtomBuilder) Coordinates(x, y, z float32) *AtomBuilder {
	a := ab.a
	a.X = x
	a.Y = y
	a.Z = z
	return ab
}

// Charge sets the residual charge on this atom.
func (ab *AtomBuilder) Charge(ch int) *AtomBuilder {
	switch ch {
	case 1:
		ab.a.charge = int8(3)
	case 2:
		ab.a.charge = int8(2)
	case 3:
		ab.a.charge = int8(1)
	case 4:
		ab.a.radical = cmn.RadicalDoublet
	case 5:
		ab.a.charge = int8(-1)
	case 6:
		ab.a.charge = int8(-2)
	case 7:
		ab.a.charge = int8(-3)
	default:
		ab.a.charge = int8(0)
	}

	return ab
}

// Valence sets the current valence configuration of this atom.
func (ab *AtomBuilder) Valence(v int) *AtomBuilder {
	if v > 0 && v < 15 {
		ab.a.valence = int8(v)
	}

	return ab
}
