package common

import (
	"fmt"
)

// Element holds the essential chemical information of a given natural
// element.
type Element struct {
	Number   uint8   // Atomic number
	Symbol   string  // Chemical symbol
	Name     string  // Element's name
	Weight   float64 // Atomic weight of the most abundant isotope
	Valence  int8    // Default valence
	OxStates []int8  // Other oxidation states
}

// String answers a representation of the element that is easily
// readable.
func (e *Element) String() string {
	return fmt.Sprintf("%s : {Number: %d, Symbol: %s, Weight: %.4f, Valence: %d, Oxidation States: %v}\n",
		e.Name, e.Number, e.Symbol, e.Weight, e.Valence, e.OxStates)
}
