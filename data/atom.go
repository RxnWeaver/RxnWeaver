package data

import (
	"github.com/RxnWeaver/rxnweaver/common"
)

const (
	MaxBonds    = 20
	MaxRings    = common.ListSizeSmall
	MaxFeatures = common.ListSizeSmall
)

// Atom represents a chemical atom.
//
// It has various physical and chemical properties, apart from the
// transient information attached to it during its participation in
// reactions.
type Atom struct {
	mol   *Molecule // Containing molecule of this atom.
	atNum uint8     // Atomic number of this atom's element.
	iId   uint16    // Serial input ID of this atom.
	nId   uint16    // Normalised ID of this atom.

	X float32 // X-coordinate of this atom.
	Y float32 // Y-coordinate of this atom.
	Z float32 // Z-coordinate of this atom.

	numH    uint8 // Total of implicit and explicit hydrogen atoms attached to this atom.
	charge  int8  // Residual net charge of this atom.
	valence int8  // Current valence configuration of this atom.

	pHash uint64 // A pseudo-hash of this atom, using some attributes.
	sHash uint64 // A pseudo-hash of this atom, using some attributes.

	bonds [MaxBonds]uint8 // Expanded list of bonds of this atom.
	nbrs  [MaxBonds]uint8 // List of distinct neighbours of this atom.

	rings [MaxRings]uint8 // List of rings this atom participates in.
	// Does this atom participate in at least one aromatic ring?
	isInAroRing bool
	// Is this atom a bridgehead of a bicyclic system of rings?
	isBridgeHead bool
	// Is this atom the sole common atom of all of its rings?
	isSpiro bool

	// The functional groups substituted on this atom.  They are listed in
	// descending order of importance.  The first is the primary feature.
	features [MaxFeatures]uint8
}
