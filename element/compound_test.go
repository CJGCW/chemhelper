package element

import (
	"fmt"
	"testing"

	"github.com/shopspring/decimal"
)

// Simplified Periodic Table with only a few elements for testing.
func NewTestPeriodicTable() *PeriodicTable {
	return newPeriodicTableFromElements([]Element{
		{AtomicNumber: 1, Symbol: "H", Name: "Hydrogen", AtomicWeight: decimal.NewFromFloat(1.008)},
		{AtomicNumber: 6, Symbol: "C", Name: "Carbon", AtomicWeight: decimal.NewFromFloat(12.011)},
		{AtomicNumber: 8, Symbol: "O", Name: "Oxygen", AtomicWeight: decimal.NewFromFloat(15.999)},
		{AtomicNumber: 11, Symbol: "Na", Name: "Sodium", AtomicWeight: decimal.NewFromFloat(22.990)},
		{AtomicNumber: 17, Symbol: "Cl", Name: "Chlorine", AtomicWeight: decimal.NewFromFloat(35.45)},
	})
}

func TestParseCompoundParentheses(t *testing.T) {
	pt := NewTestPeriodicTable()

	tests := []struct {
		input         string
		expectedError bool
		expected      map[string]int64 // symbol -> count
	}{
		{
			input:    "Ca(OH)2",
			expected: map[string]int64{"O": 2, "H": 2},
			// Ca not in test table, expect error
			expectedError: true,
		},
		{
			input:         "H2O",
			expected:      map[string]int64{"H": 2, "O": 1},
			expectedError: false,
		},
		{
			input:         "HOH",
			expected:      map[string]int64{"H": 2, "O": 1},
			expectedError: false,
		},
		{
			input:         "NaCl",
			expected:      map[string]int64{"Na": 1, "Cl": 1},
			expectedError: false,
		},
		{
			input:         "C6H12O6",
			expected:      map[string]int64{"C": 6, "H": 12, "O": 6},
			expectedError: false,
		},
		{
			input:         "(OH)2",
			expected:      map[string]int64{"O": 2, "H": 2},
			expectedError: false,
		},
		{
			input:         "H(OH)2",
			expected:      map[string]int64{"H": 3, "O": 2},
			expectedError: false,
		},
		{
			input:         "",
			expectedError: true,
		},
		{
			input:         "H2O1X",
			expectedError: true,
		},
		{
			input:         "H(OH",
			expectedError: true,
		},
	}

	for _, test := range tests {
		t.Run(fmt.Sprintf("Parsing:%q", test.input), func(t *testing.T) {
			result, err := ParseCompoundElements(test.input, pt)
			if (err != nil) != test.expectedError {
				t.Fatalf("expected error=%v, got err=%v", test.expectedError, err)
			}
			if test.expectedError {
				return
			}
			if len(result) != len(test.expected) {
				t.Fatalf("expected %d elements, got %d", len(test.expected), len(result))
			}
			for _, em := range result {
				expectedCount, ok := test.expected[em.Element.Symbol]
				if !ok {
					t.Errorf("unexpected element %s in result", em.Element.Symbol)
					continue
				}
				if !em.Moles.Equal(decimal.NewFromInt(expectedCount)) {
					t.Errorf("element %s: expected %d moles, got %s", em.Element.Symbol, expectedCount, em.Moles)
				}
			}
		})
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
