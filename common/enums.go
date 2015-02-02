package common

// The following `enum` definitions are in line with the corresponding
// ones in InChI 1.04 software.  A notable difference is that we DO
// NOT provide for specifying bond stereo with respect to the second
// atom in the pair.

// Radical represents possible radical configurations of an atom.
type Radical uint8

const (
	RadicalNone Radical = iota
	RadicalSinglet
	RadicalDoublet
	RadicalTriplet
)

// BondType defines the possible types of bonds between a pair of
// atoms.
type BondType uint8

const (
	BondTypeNone BondType = iota
	BondTypeSingle
	BondTypeDouble
	BondTypeTriple
	BondTypeAltern // InChI says 'avoid by all means'!
)

// BondStereo defines the possible stereo orientations of a given
// bond, when 2-D coordinates are given.
type BondStereo uint8

const (
	BondStereoNone         BondStereo = 0
	BondStereoUp           BondStereo = 1
	BondStereoEither       BondStereo = 4
	BondStereoDown         BondStereo = 6
	BondStereoDoubleEither BondStereo = 3
)

// StereoType specifies the nature of the origin of the stereo
// behaviour.
type StereoType uint8

const (
	StereoTypeNone StereoType = iota
	StereoTypeDoubleBond
	StereoTypeTetrahedral
	StereoTypeAllene
)

// StereoParity defines the possible stereo configurations, given a
// particular stereo centre (atom or bond).
type StereoParity uint8

const (
	StereoParityNone StereoParity = iota
	StereoParityOdd
	StereoParityEven
	StereoParityUnknown
	StereoParityUndefined
)
