package element

import (
	"fmt"

	"chemhelper/units"

	"github.com/shopspring/decimal"
)

func (c *Compound) getMoles(mass decimal.Decimal) error {
	if mass.Equal(decimal.Zero) {
		return fmt.Errorf("no mass passed")
	}
	if c.MolarMass.Equal(decimal.Zero) {
		if err := c.getMolarMass(); err != nil {
			return err
		}
	}
	c.Moles = mass.Div(c.MolarMass)
	return nil
}

func (element *ElementMoles) getMoles(mass units.Mass) error {
	moles, err := mass.GetMoles(element.Element.AtomicWeight)
	if err != nil {
		return err
	}
	element.Moles = moles
	return nil
}

func (compound *Compound) getMolesFromMass(mass units.Mass) error {
	if compound.MolarMass.Equal(decimal.Zero) {
		if err := compound.getMolarMass(); err != nil {
			return err
		}
	}
	moles, err := mass.GetMoles(compound.MolarMass)
	if err != nil {
		return err
	}
	compound.Moles = moles
	return nil
}
