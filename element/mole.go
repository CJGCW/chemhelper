package element

import (
	"fmt"

	"github.com/shopspring/decimal"
)

func getMoles(p Property, value decimal.Decimal) (decimal.Decimal, error){
	return p.getMoles(value)
}

func (m Mass) getMoles(molarMass decimal.Decimal) (decimal.Decimal, error) {
	if molarMass.Equal(decimal.Zero){
		return decimal.Zero, fmt.Errorf("no molar mass passed")
	}
	standardMass, err := m.convertToStandard()
	if err != nil {
		return decimal.Zero, err
	}
	return standardMass.Div(molarMass), nil
}

func (c *Compound) getMoles(mass decimal.Decimal) ( error) {
	if mass.Equal(decimal.Zero) {
		return fmt.Errorf("no mass passed")
	}
	if c.MolarMass.Equal(decimal.Zero) {
		err := c.getMolarMass()
		if (err != nil) {
			return err
		}	
	}
	c.Moles = mass.Div(c.MolarMass)
	return nil
}

func (v Volume) getMoles(molarity decimal.Decimal) (decimal.Decimal, error) {
	zed := decimal.NewFromInt(0)
	if molarity.LessThanOrEqual(zed) {
		return zed, fmt.Errorf("molarity must be a nonzero, positive value, got %v ",molarity)
	}
	standardVol, err := v.convertToStandard()
	if err != nil {
		return zed, err
	}
	return standardVol.Mul(molarity), err
}

func (element *ElementMoles) getMoles(mass Mass) error {
	
	moles, err := mass.getMoles(element.Element.AtomicWeight)
	if err != nil {
		return err
	} 
	element.Moles = moles
	return nil
}

func (compound *Compound) getMolesFromMass(mass Mass) error {
	if compound.MolarMass.Equal(decimal.Zero) {
		err := compound.getMolarMass()
		if (err != nil) {
			return err
		}	
	}
	moles, err := mass.getMoles(compound.MolarMass)
	if (err != nil) {
		return err
	}	
	compound.Moles = moles
	return nil
}