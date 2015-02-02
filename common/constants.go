package common

// Program-wide configuration constants.

const (
	IoBufferSize = 64 * 1024 // For streams.

	ListSizeTiny   = 2  // For reactant and product lists, etc.
	ListSizeSmall  = 10 // For functional group lists, etc.
	ListSizeMedium = 20 // For neighbour lists, etc.
	ListSizeLarge  = 64 // For atom and bond lists, etc.
)
