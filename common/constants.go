package common

// Program-wide configuration constants.

const (
	IoBufferSize = 64 * 1024 // For streams.

	ListSizeTiny   = 2  // For reactant and product lists, etc.
	ListSizeSmall  = 10 // For functional group lists, etc.
	ListSizeMedium = 20 // For neighbour lists, etc.
	ListSizeLarge  = 64 // For atom and bond lists, etc.

	MaxBonds    = 20            // Maximum number of bonds an atom can have.
	MaxRings    = ListSizeSmall // Maximum number of rings an atom can be a part of.
	MaxFeatures = ListSizeSmall // Maximum number functional groups on an atom.
)
