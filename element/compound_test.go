package element

import (
	"fmt"
	"testing"

	"github.com/shopspring/decimal"
)

// Simplified Periodic Table with only a few elements for testing.
func NewTestPeriodicTable() *PeriodicTable {
	return &PeriodicTable{
		Elements: []Element{
			{AtomicNumber: 1, Symbol: "H", Name: "Hydrogen", AtomicWeight: decimal.NewFromFloat(1.008)},
			{AtomicNumber: 6, Symbol: "C", Name: "Carbon", AtomicWeight: decimal.NewFromFloat(12.011)},
			{AtomicNumber: 8, Symbol: "O", Name: "Oxygen", AtomicWeight: decimal.NewFromFloat(15.999)},
			{AtomicNumber: 11, Symbol: "Na", Name: "Sodium", AtomicWeight: decimal.NewFromFloat(22.990)},
			{AtomicNumber: 17, Symbol: "Cl", Name: "Chlorine", AtomicWeight: decimal.NewFromFloat(35.45)},
		},
	}
}

func TestParseCompound(t *testing.T) {
	pt := NewTestPeriodicTable()

	for _, test := range TestCompounds { // test compounds were generated in property_test
		t.Run(fmt.Sprintf("Testing Compound:%s", test.compound.Symbol), func(t *testing.T) {
			result, err := ParseCompoundElements(test.compound.Symbol, pt)
			if (err != nil) != test.expectedError {
				t.Errorf("Expected error: %v, but got: %v", test.expectedError, err)
			}

			// If no error is expected, compare the result to the expected value
			if !test.expectedError {
				if len(result) != len(test.compound.Elements) {
					t.Errorf("Expected %d elements, but got %d", len(test.compound.Elements), len(result))
				}

				// Sort both the result and the expected elements to ensure order-independence
				sortElementMoles(result)
				sortElementMoles(test.compound.Elements)

				for i, elem := range result {
					if elem.Element.Symbol != test.compound.Elements[i].Element.Symbol || !elem.Moles.Equal(test.compound.Elements[i].Moles) {
						t.Errorf("Expected ElementMoles {%s, %v}, but got {%s, %v}",
							test.compound.Elements[i].Element.Symbol, test.compound.Elements[i].Moles,
							elem.Element.Symbol, elem.Moles)
					}
				}
			}
		})
	}
}
