package units

import (
	"fmt"
	"log"

	"github.com/shopspring/decimal"
)

type Prefix float64

const (
	None  Prefix = 1
	Kilo  Prefix = 1000
	Hecto Prefix = 100
	Deca  Prefix = 10
	Deci  Prefix = 0.1
	Centi Prefix = 0.01
	Milli Prefix = 0.001
	Micro Prefix = 0.000001
)

type MassUnit float64

const (
	UnknownMass MassUnit = -1
	Gram        MassUnit = 1
	Ounce       MassUnit = 28.349
	Pound       MassUnit = 453.592
)

type Mass struct {
	value  decimal.Decimal
	unit   MassUnit
	prefix Prefix
}

// Value returns the raw numeric value of the mass.
func (m Mass) Value() decimal.Decimal {
	return m.value
}

func NewMass(value decimal.Decimal, options ...interface{}) (Mass, error) {
	if value.Equal(decimal.Zero) {
		return Mass{}, fmt.Errorf("no mass value passed")
	}
	mass := Mass{
		value:  value,
		unit:   Gram,
		prefix: None,
	}
	for _, opt := range options {
		switch v := opt.(type) {
		case MassUnit:
			mass.unit = v
		case Prefix:
			mass.prefix = v
		default:
			log.Printf("%v is unexpected", v)
		}
	}
	return mass, nil
}

func (m Mass) ConvertToStandard() (decimal.Decimal, error) {
	if m.unit == UnknownMass {
		return decimal.Zero, fmt.Errorf("unknown mass unit")
	}
	if m.value.Equal(decimal.Zero) {
		return decimal.Zero, fmt.Errorf("empty property passed")
	}
	return m.value.
		Mul(decimal.NewFromFloat(float64(m.unit))).
		Mul(decimal.NewFromFloat(float64(m.prefix))), nil
}

// GetMoles returns moles given a molar mass (g/mol).
func (m Mass) GetMoles(molarMass decimal.Decimal) (decimal.Decimal, error) {
	if molarMass.Equal(decimal.Zero) {
		return decimal.Zero, fmt.Errorf("no molar mass passed")
	}
	standardMass, err := m.ConvertToStandard()
	if err != nil {
		return decimal.Zero, err
	}
	return standardMass.Div(molarMass), nil
}

type VolumeUnit float64

const (
	Liter VolumeUnit = 1
	// Room to add FluidOunce, Gallon, etc. in the future.
)

type Volume struct {
	value  decimal.Decimal
	unit   VolumeUnit
	prefix Prefix
}

func NewVolume(value decimal.Decimal, options ...interface{}) (Volume, error) {
	if value.Equal(decimal.Zero) {
		return Volume{}, fmt.Errorf("no volume value passed")
	}
	vol := Volume{
		value:  value,
		unit:   Liter,
		prefix: None,
	}
	for _, opt := range options {
		switch v := opt.(type) {
		case VolumeUnit:
			vol.unit = v
		case Prefix:
			vol.prefix = v
		default:
			log.Printf("%v is unexpected", v)
		}
	}
	return vol, nil
}

func (v Volume) ConvertToStandard() (decimal.Decimal, error) {
	if v.value.Equal(decimal.Zero) {
		return decimal.Zero, fmt.Errorf("empty property passed")
	}
	return v.value.
		Mul(decimal.NewFromFloat(float64(v.unit))).
		Mul(decimal.NewFromFloat(float64(v.prefix))), nil
}

// GetMoles returns moles given a molarity (mol/L).
func (v Volume) GetMoles(molarity decimal.Decimal) (decimal.Decimal, error) {
	if molarity.LessThanOrEqual(decimal.Zero) {
		return decimal.Zero, fmt.Errorf("molarity must be a nonzero, positive value, got %v", molarity)
	}
	standardVol, err := v.ConvertToStandard()
	if err != nil {
		return decimal.Zero, err
	}
	return standardVol.Mul(molarity), nil
}

// Property is implemented by any measurement that can be converted to a
// standard unit and used to derive moles.
type Property interface {
	ConvertToStandard() (decimal.Decimal, error)
	GetMoles(decimal.Decimal) (decimal.Decimal, error)
}

// GetMoles dispatches to the correct GetMoles implementation for a Property.
// For Mass, value is the molar mass (g/mol).
// For Volume, value is the molarity (mol/L).
func GetMoles(p Property, value decimal.Decimal) (decimal.Decimal, error) {
	return p.GetMoles(value)
}

// ConvertToStandard is a convenience wrapper around the Property interface.
func ConvertToStandard(p Property) (decimal.Decimal, error) {
	return p.ConvertToStandard()
}
