package element

import (
	"fmt"
	"testing"

	"chemhelper/units"

	"github.com/shopspring/decimal"
)

func TestMolesByMass(t *testing.T) {
	for _, test := range TestCompounds {
		t.Run(fmt.Sprintf("Testing Compound:%s", test.compound.Symbol), func(t *testing.T) {
			expected := test.expectedMoles
			err := test.compound.getMolesFromMass(test.massForMoles)
			if !test.expectedError && err != nil {
				t.Errorf("unexpected error: %v", err)
			}
			actual := test.compound.Moles
			if !actual.Equal(expected) {
				t.Errorf("Expected %v moles, but got %v", expected, actual)
			}
			if test.expectedError && err == nil {
				t.Errorf("Expected error but got none")
			}
		})
	}
}

func TestMolesByVolume(t *testing.T) {
	for _, test := range TestCompounds {
		t.Run(fmt.Sprintf("Testing Compound:%s", test.compound.Symbol), func(t *testing.T) {
			expected := test.expectedMoles
			actualMoles, err := test.compound.Volume.GetMoles(test.molarity)
			if !test.expectedError && err != nil {
				t.Errorf("unexpected error: %v", err)
			}
			if !actualMoles.Equal(expected) {
				t.Errorf("Expected %v moles, but got %v", expected, actualMoles)
			}
			if test.expectedError && err == nil {
				t.Errorf("Expected error but got none")
			}
		})
	}
}

func TestMolesFromStandardMass(t *testing.T) {
	for _, test := range TestCompounds {
		t.Run(fmt.Sprintf("Testing Compound:%s", test.compound.Symbol), func(t *testing.T) {
			expected := test.expectedMoles
			test.compound.MolarMass = test.compound.Mass.Value()
			standardMass, _ := test.massForMoles.ConvertToStandard()
			err := test.compound.getMoles(standardMass)
			if !test.expectedError && err != nil {
				t.Errorf("unexpected error: %v", err)
			}
			actual := test.compound.Moles
			if !actual.Equal(expected) {
				t.Errorf("Expected %v moles, but got %v", expected, actual)
			}
			if test.expectedError && err == nil {
				t.Errorf("Expected error but got none")
			}
		})
	}
}

func TestMolesOfElements(t *testing.T) {
	var preciseOMoles, _ = decimal.NewFromString("17.7192324520282518")

	testElementMoles := []struct {
		element       ElementMoles
		expectedError bool
		massForMoles  units.Mass
		expectedMoles decimal.Decimal
	}{
		{
			element:       ElementMoles{Element: Element{AtomicNumber: 1, Symbol: "H", Name: "Hydrogen", AtomicWeight: decimal.NewFromFloat(1.008)}},
			massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(10.08)); return m }(),
			expectedMoles: decimal.NewFromFloat(10),
			expectedError: false,
		},
		{
			element:       ElementMoles{Element: Element{AtomicNumber: 6, Symbol: "C", Name: "Carbon", AtomicWeight: decimal.NewFromFloat(12.011)}},
			massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(7.2066), units.Kilo); return m }(),
			expectedMoles: decimal.NewFromFloat(600),
			expectedError: false,
		},
		{
			element:       ElementMoles{Element: Element{AtomicNumber: 8, Symbol: "O", Name: "Oxygen", AtomicWeight: decimal.NewFromFloat(15.999)}},
			massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(10), units.Ounce); return m }(),
			expectedMoles: preciseOMoles,
			expectedError: false,
		},
		{
			element:       ElementMoles{Element: Element{AtomicNumber: 211, Symbol: "XX", Name: "Baddium", AtomicWeight: decimal.NewFromFloat(453.592)}},
			massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(100), units.Pound); return m }(),
			expectedMoles: decimal.NewFromFloat(100),
			expectedError: false,
		},
		{
			element:       ElementMoles{Element: Element{AtomicNumber: 17, Symbol: "Cl", Name: "Chlorine", AtomicWeight: decimal.NewFromFloat(35.45)}},
			massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(3.545), units.Milli); return m }(),
			expectedMoles: decimal.NewFromFloat(0.0001),
			expectedError: false,
		},
	}
	for _, test := range testElementMoles {
		t.Run(fmt.Sprintf("Testing Element:%s", test.element.Element.Symbol), func(t *testing.T) {
			err := test.element.getMoles(test.massForMoles)
			if !test.expectedError && err != nil {
				t.Errorf("unexpected error: %v", err)
			}
			actual := test.element.Moles
			if !actual.Equal(test.expectedMoles) {
				t.Errorf("Expected %v moles, but got %v", test.expectedMoles, actual)
			}
			if test.expectedError && err == nil {
				t.Errorf("Expected error but got none")
			}
		})
	}
}
