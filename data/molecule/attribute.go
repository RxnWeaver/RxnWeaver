package molecule

// Attribute represents a (key, value) pair that annotates this
// molecule.
//
// A given molecule can have zero or more such attributes.
type Attribute struct {
	Name  string
	Value string
}
