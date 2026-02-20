package element

import "testing"

func TestSigFigs(t *testing.T) {
	var testMasses = []struct {
		mass string
		expectedSigFigs int
		expectedError bool
	
	}{
		{
			mass: "1",
			expectedSigFigs: 1,
			expectedError: false,
		},
		{
			mass: "10000000",
			expectedSigFigs: 1,
			expectedError: false,
		},
		{
			mass: "100.00",
			expectedSigFigs: 5,
			expectedError: false,
		},{
			mass: "1 pound",
			expectedSigFigs: 1,
			expectedError: false,
		},{
			mass: "0.0",
			expectedSigFigs: 0,
			expectedError: false,
		},
		{
			mass: "fish",
			expectedSigFigs: 0,
			expectedError: true,
		},
	}
	for _, test := range testMasses {
		result, err := GetSignificantFigures(test.mass)
		if !test.expectedError && (err !=nil)  {
			t.Errorf("Unexpected error for %s: %s", test.mass, err)
		}
		if test.expectedError && err == nil{
			t.Errorf("Expected error for %s but got none", test.mass)
		}
		if !test.expectedError && result != test.expectedSigFigs {
			t.Errorf("Test %s failed: expected %d, got %d", test.mass, test.expectedSigFigs, result)
		}
	}
}
func TestLowestSigFigs(t *testing.T){
	var testMassGroup = []struct {
		masses []string
		expectedSigFigs int
		expectedError bool
	
	}{
		{
			masses: []string{"10","1000","10.0"},
			expectedSigFigs: 1,
			expectedError: false,
		},
		{
			masses: []string{"10.0","1000.0","100.0"},
			expectedSigFigs: 3,
			expectedError: false,
		},
		{
			masses: []string{"12345","999999999","3.14159265"},
			expectedSigFigs: 5,
			expectedError: false,
		},{
			masses: []string{"10.0"},
			expectedSigFigs: 3,
			expectedError: false,
		},{
			masses: []string{},
			expectedSigFigs: 0,
			expectedError: true,
		},
		{
			masses: []string{"10","1000","chips"},
			expectedSigFigs: 0,
			expectedError: true,
		},
	}
	for _, test := range testMassGroup {
		result, err := GetLowestSignificantFigures(test.masses)
		if !test.expectedError && (err !=nil)  {
			t.Errorf("Unexpected error for %s: %s", test.masses, err)
		}
		if test.expectedError && err == nil{
			t.Errorf("Expected error for %s but got none", test.masses)
		}
		if !test.expectedError && result != test.expectedSigFigs {
			t.Errorf("Test %s failed: expected %d, got %d", test.masses, test.expectedSigFigs, result)
		}
	}
}
func TestSetToSigFigs(t *testing.T){
	var tests = []struct {
		name string
		value float64
		sigfigs int32
		expectedResult float64
		expectedError bool

	}{
		{
			name: "general rounding",
			value: 123.45,
			sigfigs: 3,
			expectedResult: 123,
			expectedError: false,
		},
		{
			name: "Round to make even, rounding down",
			value: 123.45,
			sigfigs: 4,
			expectedResult: 123.4, // round to make even
			expectedError: false,
		},
		{
			name: "Round to make even, rounding up",
			value: 123.35,
			sigfigs: 4,
			expectedResult: 123.4,
			expectedError: false,
		},
		{
			name: "Rounding ints",
			value: 12345,
			sigfigs: 2,
			expectedResult: 12000,
			expectedError: false,
		},
		{
			name: "Significant zeros",
			value: 10.02,
			sigfigs: 3,
			expectedResult: 10.0,
			expectedError: false,
		},
		{
			name: "Less precision than desired sig figs",
			value: 10.0213,
			sigfigs: 30,
			expectedResult: 10.0213,
			expectedError: false,
		},
	}
	for _, test := range tests {
		result, err := SetToSigFigs(test.value, test.sigfigs)
		if !test.expectedError && (err !=nil)  {
			t.Errorf("Unexpected error for %s: %s", test.name, err)
		}
		if test.expectedError && err == nil{
			t.Errorf("Expected error for %v but got none", test.name)
		}
		if !test.expectedError && result != test.expectedResult {
			t.Errorf("Test %s failed: expected %v, got %v", test.name, test.expectedResult, result)
		}
	}
}