package data

import (
	"sync"
)

// nextMolIdHolder is a synchronised struct used to assign a
// globally-unique ID to each molecule.
type nextMolIdHolder struct {
	mu     sync.Mutex
	nextId uint64
}

var nextMolId nextMolIdHolder

func nextMoleculeId() uint64 {
	nextMolId.mu.Lock()
	defer nextMolId.mu.Unlock()

	nextMolId.nextId++
	return nextMolId.nextId
}

// Molecule represents a chemical molecule.
//
// It holds information concerning its atom, bonds, rings, etc.  Note
// that a molecule is expected to be a single connected component.
type Molecule struct {
	id uint32 // The globally-unique ID of this molecule.

	atoms []*Atom // List of atoms in this molecule.
	bonds []*Bond // List of bonds in this molecule.
}
