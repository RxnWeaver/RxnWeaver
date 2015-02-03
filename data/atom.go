package data

import (
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

	numH    uint8 // Number of implicit + explicit H atoms attached to this atom.
	charge  int8  // Residual net charge of this atom.
	valence int8  // Current valence configuration of this atom.

	pHash uint64 // A pseudo-hash of this atom, using some attributes.
	sHash uint64 // A pseudo-hash of this atom, using some attributes.

	bonds          [cmn.MaxBonds]uint8 // Expanded list of bonds of this atom.
	nbrs           [cmn.MaxBonds]uint8 // List of distinct neighbours of this atom.
	numSingleBonds uint8               // Number of single bonds this atom has.
	numDoubleBonds uint8               // Number of double bonds this atom has.
	numTripleBonds uint8               // Number of triple bonds this atom has.

	rings [cmn.MaxRings]uint8 // List of rings this atom participates in.
	// Does this atom participate in at least one aromatic ring?
	isInAroRing bool
	// Is this atom a bridgehead of a bicyclic system of rings?
	isBridgeHead bool
	// Is this atom the sole common atom of all of its rings?
	isSpiro bool

	// The functional groups substituted on this atom.  They are listed in
	// descending order of importance.  The first is the primary feature.
	features [cmn.MaxFeatures]uint8
}

// newAtom constructs and initialises a new atom of the given element
// type, and belonging to the given molecule.
func newAtom(mol *Molecule, atNum uint8) *Atom {
	atom := new(Atom)
	atom.mol = mol
	atom.atNum = atNum
	atom.valence = cmn.PeriodicTable[cmn.ElemSyms[atNum]].Valence

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

// numPiElectrons answers the number of delocalised pi electrons
// contributed by this atom.  This number is important for calculating
// the aromaticity of the rings this atom participates in.
func (a *Atom) numPiElectrons() int {
	wtSum := 100*int8(a.numDoubleBonds) + 10*int8(a.numSingleBonds) + a.charge

	switch a.atNum {
	case 6:
		switch wtSum {
		case 19:
			return 2
		case 110:
			return 1
		case 120:
			var b *Bond
			for _, idx := range a.bonds {
				b = a.mol.bonds[idx-1]
				if b.bType == cmn.BondTypeDouble {
					break
				}
			}
			if len(b.rings) > 0 {
				return 1
			} else {
				return 0
			}
		default:
			return 0
		}
	}

	return 0
}
