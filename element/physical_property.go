package element

import (
	"fmt"
)

// getMolarMass computes the molar mass of a compound from its constituent
// elements and stores it in compound.MolarMass.
func (compound *Compound) getMolarMass() error {
	if len(compound.Elements) == 0 {
		return fmt.Errorf("no elements passed")
	}
	for _, em := range compound.Elements {
		compound.MolarMass = compound.MolarMass.Add(em.Element.AtomicWeight.Mul(em.Moles))
	}
	return nil
}
