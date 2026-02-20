package element

import (
	"fmt"
	"testing"

	"github.com/shopspring/decimal"
)

func TestConvertMassToStandard(t *testing.T) {	
	tests:= []struct{
		name		string
		value		float64
		unit		MassUnit
		prefix		Prefix
		expectedResult float64
		expectedError bool
	}{

		{
			name:           "1 kilogram to grams",
			value:          1,
			unit:           gram,
			prefix:         kilo,
			expectedResult: 1000,
			expectedError:  false,
		},
		{
			name:           "1 hectogram to grams",
			value:          1,
			unit:           gram,
			prefix:         hecto,
			expectedResult: 100,
			expectedError:  false,
		},
		{
			name:           "2 pounds to grams",
			value:          2,
			unit:           pound,
			prefix:         none,
			expectedResult: 907.184, // 2 * 453.592
			expectedError:  false,
		},
		{
			name:           "3 ounces to grams",
			value:          3,
			unit:           ounce,
			prefix:         none,
			expectedResult: 85.047, // 3 * 28.349
			expectedError:  false,
		},
		{
			name:           "100 milligrams to grams",
			value:          100,
			unit:           gram,
			prefix:         milli,
			expectedResult: 0.1,
			expectedError:  false,
		},
		{
			name:           "1 micrograms to grams",
			value:          1,
			unit:           gram,
			prefix:         micro,
			expectedResult: .000001,
			expectedError:  false,
		},
		{
			name:           "1000 pounds to grams",
			value:          1000,
			unit:           pound,
			prefix:         none,
			expectedResult: 453592,
			expectedError:  false,
		},
		{
			name:           "Invalid unit (unknown)",
			value:          2,
			unit:           unknownMass,
			prefix:         none,
			expectedResult: 2,
			expectedError:  true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			m := Mass{value: decimal.NewFromFloat(test.value), unit: test.unit, prefix: test.prefix}
			result, err := convertToStandardValue(m)
			if !test.expectedError && err !=nil  {
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
	volumeTests:= []struct{
		name	string
		prefix		Prefix
		value		float64
		expectedResult float64
		expectedError bool
	}{
		{prefix: none, value: 2, expectedResult: 2, expectedError:  false}, 
		{prefix: kilo, value: 8.25, expectedResult:  8250, expectedError:  false},
		{prefix: hecto, value: 9, expectedResult:  900, expectedError:  false},
		{prefix: deca, value: 10, expectedResult:  100, expectedError:  false},
		{prefix: deci, value: 1005, expectedResult:  100.5, expectedError:  false},
		{prefix: centi, value: 888, expectedResult:  8.88, expectedError:  false},
		{prefix: milli, value: 1618, expectedResult:  1.618, expectedError:  false},
		{prefix: micro, value: 1, expectedResult:  0.000001, expectedError:  false},
	}
	for _, test := range volumeTests {
		t.Run(test.name, func(t *testing.T) {
			v := Volume{value: decimal.NewFromFloat(test.value), unit: test.prefix}
			result, err := convertToStandardValue(v)
			if !test.expectedError && err !=nil  {
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

func TestMolarMass(t *testing.T) {
	for _, test := range TestCompounds {
		t.Run(fmt.Sprintf("Testing Compound:%s", test.compound.Symbol), func(t *testing.T) {
			expected := test.compound.Mass.value
			err := test.compound.getMolarMass()
			if !test.expectedError && err!=nil {
				t.Errorf("unexpected error")
			}
			actual := test.compound.MolarMass
			if !actual.Equal(expected){
				t.Errorf("Expected %v mass, but got %v", expected, actual)
			}
			if (test.expectedError && err == nil) {
				t.Errorf("Expected error but got none")
			}
		})
	}
}

func TestNewMass(t *testing.T){
	var gram MassUnit = gram
	var pound MassUnit = pound
	var kilo Prefix = kilo
	
	tests := []struct {
		name		string
		value float64
		unit *MassUnit
		prefix *Prefix
		expectedResult float64
		expectedError bool
	}{
		{
			name:           "1 kilogram to grams",
			value:          1,
			unit:           &gram,
			prefix:         &kilo,
			expectedResult: 1000,
			expectedError:  false,
		},
		{
			name:           ".5 kilogram with no unit to grams",
			value:          .5,
			unit:           nil,
			prefix:         &kilo,
			expectedResult: 500,
			expectedError:  false,
		},
		{
			name:           "1 pounds to grams no prefix",
			value:          1,
			unit:           &pound,
			prefix:         nil,
			expectedResult: 453.592,
			expectedError:  false,
		},
		{
			name:           "10 grams only passing value",
			value:          10,
			unit:           nil,
			prefix:         nil,
			expectedResult: 10,
			expectedError:  false,
		},
		{
			name:           "Invalid value",
			value:          0,
			unit:           &gram,
			prefix:         &kilo,
			expectedResult: 1000,
			expectedError:  true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			var m Mass
			var e error
			value := decimal.NewFromFloat(test.value)
			if test.unit != nil && test.prefix !=nil {
				m,e =  NewMass(value, *test.unit ,*test.prefix)
			} else if test.unit == nil && test.prefix !=nil {
				m,e =  NewMass(value, *test.prefix)
			} else if test.unit != nil && test.prefix ==nil {
				m,e =  NewMass(value, *test.unit)
			} else if test.unit == nil && test.prefix ==nil {
				m,e =  NewMass(value)
			} else {
				t.Error("Unexpected testing condition")
			}

			result, err := convertToStandardValue(m)
			if !test.expectedError && (err !=nil || e !=nil)  {
				t.Errorf("Unexpected error for %s: %s", test.name, err)
			}
			if test.expectedError && (err == nil && e == nil){
				t.Errorf("Expected error for %s but got %v", test.name, result)
			}
			if !test.expectedError && !result.Equal(decimal.NewFromFloat(test.expectedResult)) {
				t.Errorf("Test %s failed: expected %v, got %v", test.name, test.expectedResult, result)
			}
		})
	}
}
