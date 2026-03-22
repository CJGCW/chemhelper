package units

import (
	"testing"

	"github.com/shopspring/decimal"
)

func TestConvertMassToStandard(t *testing.T) {
	tests := []struct {
		name           string
		value          float64
		unit           MassUnit
		prefix         Prefix
		expectedResult float64
		expectedError  bool
	}{
		{name: "1 kilogram to grams", value: 1, unit: Gram, prefix: Kilo, expectedResult: 1000, expectedError: false},
		{name: "1 hectogram to grams", value: 1, unit: Gram, prefix: Hecto, expectedResult: 100, expectedError: false},
		{name: "2 pounds to grams", value: 2, unit: Pound, prefix: None, expectedResult: 907.184, expectedError: false},
		{name: "3 ounces to grams", value: 3, unit: Ounce, prefix: None, expectedResult: 85.047, expectedError: false},
		{name: "100 milligrams to grams", value: 100, unit: Gram, prefix: Milli, expectedResult: 0.1, expectedError: false},
		{name: "1 microgram to grams", value: 1, unit: Gram, prefix: Micro, expectedResult: .000001, expectedError: false},
		{name: "1000 pounds to grams", value: 1000, unit: Pound, prefix: None, expectedResult: 453592, expectedError: false},
		{name: "Invalid unit (unknown)", value: 2, unit: UnknownMass, prefix: None, expectedResult: 2, expectedError: true},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			m := Mass{value: decimal.NewFromFloat(test.value), unit: test.unit, prefix: test.prefix}
			result, err := ConvertToStandard(m)
			if !test.expectedError && err != nil {
				t.Errorf("Unexpected error for %s: %s", test.name, err)
			}
			if test.expectedError && result.Equal(decimal.NewFromFloat(test.expectedResult)) {
				t.Errorf("Expected error for %s but got %v", test.name, result)
			}
			if !test.expectedError && !result.Equal(decimal.NewFromFloat(test.expectedResult)) {
				t.Errorf("Test %s failed: expected %v, got %v", test.name, test.expectedResult, result)
			}
		})
	}
}

func TestConvertVolume(t *testing.T) {
	volumeTests := []struct {
		name           string
		prefix         Prefix
		value          float64
		expectedResult float64
		expectedError  bool
	}{
		{prefix: None, value: 2, expectedResult: 2, expectedError: false},
		{prefix: Kilo, value: 8.25, expectedResult: 8250, expectedError: false},
		{prefix: Hecto, value: 9, expectedResult: 900, expectedError: false},
		{prefix: Deca, value: 10, expectedResult: 100, expectedError: false},
		{prefix: Deci, value: 1005, expectedResult: 100.5, expectedError: false},
		{prefix: Centi, value: 888, expectedResult: 8.88, expectedError: false},
		{prefix: Milli, value: 1618, expectedResult: 1.618, expectedError: false},
		{prefix: Micro, value: 1, expectedResult: 0.000001, expectedError: false},
	}
	for _, test := range volumeTests {
		t.Run(test.name, func(t *testing.T) {
			v := Volume{value: decimal.NewFromFloat(test.value), unit: Liter, prefix: test.prefix}
			result, err := ConvertToStandard(v)
			if !test.expectedError && err != nil {
				t.Errorf("Unexpected error for %s: %s", test.name, err)
			}
			if test.expectedError && result.Equal(decimal.NewFromFloat(test.expectedResult)) {
				t.Errorf("Expected error for %s but got %v", test.name, result)
			}
			if !test.expectedError && !result.Equal(decimal.NewFromFloat(test.expectedResult)) {
				t.Errorf("Test %s failed: expected %v, got %v", test.name, test.expectedResult, result)
			}
		})
	}
}

func TestNewMass(t *testing.T) {
	var gram MassUnit = Gram
	var pound MassUnit = Pound
	var kilo Prefix = Kilo

	tests := []struct {
		name           string
		value          float64
		unit           *MassUnit
		prefix         *Prefix
		expectedResult float64
		expectedError  bool
	}{
		{name: "1 kilogram to grams", value: 1, unit: &gram, prefix: &kilo, expectedResult: 1000, expectedError: false},
		{name: ".5 kilogram with no unit to grams", value: .5, unit: nil, prefix: &kilo, expectedResult: 500, expectedError: false},
		{name: "1 pounds to grams no prefix", value: 1, unit: &pound, prefix: nil, expectedResult: 453.592, expectedError: false},
		{name: "10 grams only passing value", value: 10, unit: nil, prefix: nil, expectedResult: 10, expectedError: false},
		{name: "Invalid value", value: 0, unit: &gram, prefix: &kilo, expectedResult: 1000, expectedError: true},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			var m Mass
			var e error
			value := decimal.NewFromFloat(test.value)
			if test.unit != nil && test.prefix != nil {
				m, e = NewMass(value, *test.unit, *test.prefix)
			} else if test.unit == nil && test.prefix != nil {
				m, e = NewMass(value, *test.prefix)
			} else if test.unit != nil && test.prefix == nil {
				m, e = NewMass(value, *test.unit)
			} else {
				m, e = NewMass(value)
			}

			result, err := ConvertToStandard(m)
			if !test.expectedError && (err != nil || e != nil) {
				t.Errorf("Unexpected error for %s: %s", test.name, err)
			}
			if test.expectedError && (err == nil && e == nil) {
				t.Errorf("Expected error for %s but got %v", test.name, result)
			}
			if !test.expectedError && !result.Equal(decimal.NewFromFloat(test.expectedResult)) {
				t.Errorf("Test %s failed: expected %v, got %v", test.name, test.expectedResult, result)
			}
		})
	}
}

func TestMolesOfProperty(t *testing.T) {
	testProperties := []struct {
		name          string
		property      Property
		value         decimal.Decimal
		expectedMoles decimal.Decimal
		expectedError bool
	}{
		{
			name:          "1 kg of 1g/mol",
			property:      Mass{value: decimal.NewFromInt(1), unit: Gram, prefix: Kilo},
			value:         decimal.NewFromInt(1),
			expectedMoles: decimal.NewFromInt(1000),
			expectedError: false,
		},
		{
			name:          "1 kg of 10g/mol",
			property:      Mass{value: decimal.NewFromInt(1), unit: Gram, prefix: Kilo},
			value:         decimal.NewFromInt(10),
			expectedMoles: decimal.NewFromInt(100),
			expectedError: false,
		},
		{
			name:          "1 lb of 1g/mol",
			property:      Mass{value: decimal.NewFromInt(1), unit: Pound, prefix: None},
			value:         decimal.NewFromInt(1),
			expectedMoles: decimal.NewFromFloat(453.592),
			expectedError: false,
		},
		{
			name:          "1 L of 1mol/L",
			property:      Volume{value: decimal.NewFromFloat(1), unit: Liter, prefix: None},
			value:         decimal.NewFromInt(1),
			expectedMoles: decimal.NewFromFloat(1),
			expectedError: false,
		},
		{
			name:          "1 mL of 10mol/L",
			property:      Volume{value: decimal.NewFromFloat(1), unit: Liter, prefix: Milli},
			value:         decimal.NewFromInt(10),
			expectedMoles: decimal.NewFromFloat(.01),
			expectedError: false,
		},
		{
			name:          "10 L of 0.05mol/L",
			property:      Volume{value: decimal.NewFromFloat(10), unit: Liter, prefix: None},
			value:         decimal.NewFromFloat(0.05),
			expectedMoles: decimal.NewFromFloat(0.5),
			expectedError: false,
		},
		{
			name:          "1 μL of 1mol/L",
			property:      Volume{value: decimal.NewFromFloat(1), unit: Liter, prefix: Micro},
			value:         decimal.NewFromInt(1),
			expectedMoles: decimal.NewFromFloat(0.000001),
			expectedError: false,
		},
		{
			name:          "0 molar mass throws error",
			property:      Mass{value: decimal.NewFromInt(1), unit: Gram, prefix: Kilo},
			value:         decimal.NewFromInt(0),
			expectedMoles: decimal.NewFromInt(0),
			expectedError: true,
		},
		{
			name:          "0 molarity throws an error",
			property:      Volume{value: decimal.NewFromFloat(100), unit: Liter, prefix: Kilo},
			value:         decimal.NewFromInt(0),
			expectedMoles: decimal.NewFromFloat(0),
			expectedError: true,
		},
		{
			name:          "0 mass throws error",
			property:      Mass{value: decimal.NewFromInt(0), unit: Gram, prefix: Kilo},
			value:         decimal.NewFromInt(1),
			expectedMoles: decimal.NewFromInt(0),
			expectedError: true,
		},
		{
			name:          "0 volume throws an error",
			property:      Volume{value: decimal.NewFromFloat(0), unit: Liter, prefix: Kilo},
			value:         decimal.NewFromInt(1),
			expectedMoles: decimal.NewFromFloat(0),
			expectedError: true,
		},
	}
	for _, testProperty := range testProperties {
		t.Run(testProperty.name, func(t *testing.T) {
			actualMoles, err := GetMoles(testProperty.property, testProperty.value)
			if !testProperty.expectedError && err != nil {
				t.Errorf("unexpected error: %v", err)
			}
			if !actualMoles.Equal(testProperty.expectedMoles) {
				t.Errorf("Expected %v moles, but got %v", testProperty.expectedMoles, actualMoles)
			}
			if testProperty.expectedError && err == nil {
				t.Errorf("Expected error but got none")
			}
		})
	}
}
