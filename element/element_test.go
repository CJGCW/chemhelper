package element

import (
	"testing"

	"github.com/shopspring/decimal"
)
var elementSymbols = []string{
	"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", 
	"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
	"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", 
	"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", 
	"Sb", "I", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
	"Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
	"Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
	"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", 
	"Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", 
	"Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", 
	"Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
}

func TestNewPeriodicTableSize(t *testing.T){
	pd := NewPeriodicTable()
	expected := 118
	if len(pd.Elements) != expected {
		t.Errorf("Expected table length of %v, but got %v", expected, len(pd.Elements))
	}
}

//We test the main organic elements since they are the most used, and most important
func TestNewPeriodicTableElements(t *testing.T){
	pd := NewPeriodicTable()
	var expectedElements = []struct {
		index   int
		element Element
	}{
		{0, newTestElement(1, "H", "Hydrogen", 1.008, 2.20, 120.0, 1, 1)},
		{5, newTestElement(6, "C", "Carbon", 12.011, 2.55, 170.0, 14, 2)},
		{6, newTestElement(7, "N", "Nitrogen", 14.007, 3.04, 155.0, 15, 2)},
		{7, newTestElement(8, "O", "Oxygen", 15.999, 3.44, 152.0, 16, 2)},
	}

	// Iterate over expected elements and compare
	for _, e := range expectedElements {
		actual := pd.Elements[e.index]
		checkElement(t, e.element, actual)
	}

}
func TestFindElementBySymbol(t *testing.T){
	pd := NewPeriodicTable()
	expected := newTestElement(46, "Pd", "Palladium", 106.42, 2.20, 163.0, 10, 5)
	actual, found := pd.FindElementBySymbol(expected.Symbol)
	if actual == nil || !found {
		t.Errorf("Expected to find %s in the periodic table ", expected.Name)
	}
	checkElement(t, expected, *actual)
	expected = newTestElement(27, "Co", "Cobalt", 58.933194, 1.88, 152.0, 9, 4)
	actual, found = pd.FindElementBySymbol(expected.Symbol)
	if actual == nil || !found {
		t.Errorf("Expected to find %s in the periodic table ", expected.Name)
	}
	checkElement(t, expected, *actual)
	for _, expectedSymbol := range elementSymbols {
		_, found = pd.FindElementBySymbol(expectedSymbol)
		if !found {
			t.Errorf("Expected to find %s in the periodic table ", expectedSymbol)
		}
	}
	unexpectedSymbol := "Xx"
	_, found = pd.FindElementBySymbol(unexpectedSymbol)
		if found {
			t.Errorf("Expected to find %s in the periodic table ", unexpectedSymbol)
		}
}

func TestElementGroups(t *testing.T) {
	var expectedElements = []struct {
		index   int
		element Element
		expectedGroup string
		expectedError bool
	}{
		{0, newTestElement(1, "H", "Hydrogen", 1.008, 2.20, 120.0, 1, 1), "None", false},
		{2, newTestElement(11, "Na", "Sodium",22.98976928, 0.93, 180.0, 1, 3), "Alkali Metals", false},
		{3, newTestElement(12, "Mg", "Magnesium", 24.305, 1.31, 173.0, 2, 3), "Alkaline Earth Metals", false},
		{5, newTestElement(6, "C", "Carbon", 12.011, 2.55, 170.0, 14, 2),"Carbon", false },
		{6, newTestElement(7, "N", "Nitrogen", 14.007, 3.04, 155.0, 15, 2), "Pnictogens", false},
		{7, newTestElement(8, "O", "Oxygen", 15.999, 3.44, 152.0, 16, 2), "Chalcogens", false},
		{35, newTestElement(35, "Br", "Bromine", 79.904, 2.96, 185.0, 17, 4), "Halogens", false},
		{10, newTestElement(10, "Ne", "Neon", 20.1797, 0.0, 154.0, 18, 2), "Noble Gases", false},
		{47, newTestElement(47, "Ag", "Silver", 107.8682, 1.93, 172.0, 11, 5), "Metals", false},
		{190, newTestElement(190, "XX", "ELEMENT X", 107.8682, 1.93, 172.0, 20, 5), "Unknown", true},
	}
	for _, e := range expectedElements {
		actual, err := e.element.GetGroup()
		if (err != nil && !e.expectedError){
			t.Errorf("Unexpected error")
		} else if (err == nil && e.expectedError){
			t.Errorf("Expected error but got none")
		} 
		if actual != e.expectedGroup {
			t.Errorf("Expected group to be %s, but got %s", e.expectedGroup, actual)
		}
		
	}
}
func TestIsTransitionMetal(t *testing.T) {
	cases := []struct {
		symbol string
		want   bool
	}{
		{"Sc", true},  // group 3, lower boundary
		{"Fe", true},  // group 8, mid-range
		{"Zn", true},  // group 12, upper boundary
		{"Ca", false}, // group 2, alkaline earth
		{"Al", false}, // group 13, post-transition
		{"O", false},  // group 16, nonmetal
		{"Ag", true},  // group 11
		{"Cu", true},  // group 11
	}
	pt := NewPeriodicTable()
	for _, c := range cases {
		elem, found := pt.FindElementBySymbol(c.symbol)
		if !found {
			t.Fatalf("element %s not found in periodic table", c.symbol)
		}
		if got := elem.IsTransitionMetal(); got != c.want {
			t.Errorf("IsTransitionMetal(%s) = %v, want %v", c.symbol, got, c.want)
		}
	}
}

func TestIsMetal(t *testing.T) {
	cases := []struct {
		symbol string
		want   bool
	}{
		{"H", false},  // group 1 by convention but is a nonmetal
		{"Na", true},  // group 1, alkali metal
		{"Ca", true},  // group 2, alkaline earth
		{"Fe", true},  // group 8, transition metal
		{"Al", true},  // group 13, post-transition (explicit)
		{"Ga", true},  // group 13, post-transition (explicit)
		{"In", true},  // group 13, post-transition (explicit)
		{"Tl", true},  // group 13, post-transition (explicit)
		{"Sn", true},  // group 14, post-transition (explicit)
		{"Pb", true},  // group 14, post-transition (explicit)
		{"Bi", true},  // group 15, post-transition (explicit)
		{"C", false},  // group 14, nonmetal
		{"O", false},  // group 16, nonmetal
		{"Cl", false}, // group 17, halogen
		{"He", false}, // group 18, noble gas
		{"Si", false}, // group 14, metalloid
		{"As", false}, // group 15, metalloid
	}
	pt := NewPeriodicTable()
	for _, c := range cases {
		elem, found := pt.FindElementBySymbol(c.symbol)
		if !found {
			t.Fatalf("element %s not found in periodic table", c.symbol)
		}
		if got := elem.IsMetal(); got != c.want {
			t.Errorf("IsMetal(%s) = %v, want %v", c.symbol, got, c.want)
		}
	}
}

func TestMonatomicIonCharge(t *testing.T) {
	cases := []struct {
		symbol string
		want   int
	}{
		{"N", -3},  // group 15
		{"P", -3},  // group 15
		{"O", -2},  // group 16
		{"S", -2},  // group 16
		{"F", -1},  // group 17
		{"Cl", -1}, // group 17
		{"Br", -1}, // group 17
		{"Na", 0},  // group 1, metal returns 0
		{"Ca", 0},  // group 2, metal returns 0
		{"Fe", 0},  // transition metal returns 0
		{"C", 0},   // group 14, returns 0
		{"He", 0},  // noble gas returns 0
	}
	pt := NewPeriodicTable()
	for _, c := range cases {
		elem, found := pt.FindElementBySymbol(c.symbol)
		if !found {
			t.Fatalf("element %s not found in periodic table", c.symbol)
		}
		if got := elem.MonatomicIonCharge(); got != c.want {
			t.Errorf("MonatomicIonCharge(%s) = %d, want %d", c.symbol, got, c.want)
		}
	}
}

func TestValence(t *testing.T) {
	cases := []struct {
		symbol string
		want   int
	}{
		{"H", 1},  // group 1
		{"Li", 1}, // group 1
		{"Ca", 2}, // group 2
		{"Mg", 2}, // group 2
		{"Fe", 8}, // group 8, TM: returns group number
		{"Cu", 11}, // group 11, TM
		{"Zn", 12}, // group 12, TM
		{"Al", 3}, // group 13
		{"C", 4},  // group 14
		{"N", 5},  // group 15
		{"O", 6},  // group 16
		{"Cl", 7}, // group 17
		{"He", 2}, // group 18, special case Z=2
		{"Ne", 8}, // group 18
		{"Ar", 8}, // group 18
	}
	pt := NewPeriodicTable()
	for _, c := range cases {
		elem, found := pt.FindElementBySymbol(c.symbol)
		if !found {
			t.Fatalf("element %s not found in periodic table", c.symbol)
		}
		if got := elem.Valence(); got != c.want {
			t.Errorf("Valence(%s) = %d, want %d", c.symbol, got, c.want)
		}
	}
}

func newTestElement(number int, symbol string, name string, weight float64, en float64, radius float64, group int, period int, ) Element{
	return Element{AtomicNumber: number, Symbol: symbol, Name: name, AtomicWeight: decimal.NewFromFloat(weight), Electronegativity: en, VanDerWaalsRadius: radius, Group: group, Period: period}
}

func checkElement(t *testing.T, expected Element, actual Element,) {
	if actual.AtomicNumber != expected.AtomicNumber && !actual.AtomicWeight.Equal(expected.AtomicWeight) && actual.Electronegativity != expected.Electronegativity && actual.Group != expected.Group && actual.Name != expected.Name && actual.Period != expected.Period && actual.Symbol != expected.Symbol && actual.VanDerWaalsRadius != expected.VanDerWaalsRadius{
		t.Errorf("Expected %s to equal %v, but got %v", expected.Name, expected, actual)
	}
}
