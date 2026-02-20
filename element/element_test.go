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
	"Rg", "Cn", "Fl", "Mc", "Lv", "Ts", "Og",
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
func newTestElement(number int, symbol string, name string, weight float64, en float64, radius float64, group int, period int, ) Element{
	return Element{AtomicNumber: number, Symbol: symbol, Name: name, AtomicWeight: decimal.NewFromFloat(weight), Electronegativity: en, VanDerWaalsRadius: radius, Group: group, Period: period}
}

func checkElement(t *testing.T, expected Element, actual Element,) {
	if actual.AtomicNumber != expected.AtomicNumber && !actual.AtomicWeight.Equal(expected.AtomicWeight) && actual.Electronegativity != expected.Electronegativity && actual.Group != expected.Group && actual.Name != expected.Name && actual.Period != expected.Period && actual.Symbol != expected.Symbol && actual.VanDerWaalsRadius != expected.VanDerWaalsRadius{
		t.Errorf("Expected %s to equal %v, but got %v", expected.Name, expected, actual)
	}
}
