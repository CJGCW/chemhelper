package element

import (
	"fmt"
	"testing"
)

func TestMolarMass(t *testing.T) {
	for _, test := range TestCompounds {
		t.Run(fmt.Sprintf("Testing Compound:%s", test.compound.Symbol), func(t *testing.T) {
			expected := test.compound.Mass.Value()
			err := test.compound.getMolarMass()
			if !test.expectedError && err != nil {
				t.Errorf("unexpected error")
			}
			actual := test.compound.MolarMass
			if !actual.Equal(expected) {
				t.Errorf("Expected %v mass, but got %v", expected, actual)
			}
			if test.expectedError && err == nil {
				t.Errorf("Expected error but got none")
			}
		})
	}
}
