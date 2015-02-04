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

// atomWithIid answers the atom for the given input ID, if found.
// Answers `nil` otherwise.
func (m *Molecule) atomWithIid(id uint16) *Atom {
	for _, a := range m.atoms {
		if a.iId == id {
			return a
		}
	}

	return nil
}

// bondBetween answers the bond between the two given atoms, if one
// such exists.  Answers `nil` otherwise.
//
// Note that the two given atoms are represented by their input IDs,
// NOT normalised IDs.
func (m *Molecule) bondBetween(a1, a2 uint16) *Bond {
	for _, b := range m.bonds {
		if (b.a1 == a1 && b.a2 == a2) || (b.a2 == a1 && b.a1 == a2) {
			return b
		}
	}

	return nil
}
