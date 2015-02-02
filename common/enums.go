package common

// Radical represents possible radical configurations of an atom.
type Radical uint8

const (
	RadicalNone Radical = 0
	RadicalSinglet
	RadicalDoublet
	RadicalTriplet
)

// BondType defines the possible types of bonds between a pair of
// atoms.
type BondType uint8

const (
	BondTypeNone BondType = 0
	BondTypeSingle
	BondTypeDouble
	BondTypeTriple
	BondTypeAltern // InChI says 'avoid by all means'!
)
