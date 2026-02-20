package element

import (
	"github.com/shopspring/decimal"
)

// Shared testing objects
var preciseHOHMoles, _ =decimal.NewFromString("1258.9286705523175132") // floats are not precise enough
var TestCompounds = []struct {
	compound     Compound
	expectedError bool
	massForMoles Mass
	molarity decimal.Decimal
	expectedMoles decimal.Decimal
}{
	{
		compound: Compound{
			Symbol: "H2O",
			Elements: []ElementMoles{
				{Element: Element{Symbol: "H", Name: "Hydrogen", AtomicNumber: 1, AtomicWeight: decimal.NewFromFloat(1.008)}, Moles: decimal.NewFromFloat(2)},
				{Element: Element{Symbol: "O", Name: "Oxygen", AtomicNumber: 8, AtomicWeight: decimal.NewFromFloat(15.999)}, Moles: decimal.NewFromFloat(1)},
			},
			Mass: Mass{value: decimal.NewFromFloat(18.015)},
			Volume: Volume{value: decimal.NewFromInt(1), unit: none },
			},
			molarity: decimal.NewFromInt(1),
			massForMoles: Mass{value: decimal.NewFromFloat(18.015), unit: gram, prefix: none},
			expectedMoles: decimal.NewFromInt(1),
			expectedError: false,

	},
	{
		compound: Compound{
			Symbol: "NaCl",
			Elements: []ElementMoles{
				{Element: Element{Symbol: "Na", Name: "Sodium", AtomicNumber: 11, AtomicWeight: decimal.NewFromFloat(22.990)}, Moles: decimal.NewFromFloat(1)},
				{Element: Element{Symbol: "Cl", Name: "Chlorine", AtomicNumber: 17, AtomicWeight: decimal.NewFromFloat(35.45)}, Moles: decimal.NewFromFloat(1)},
			},
			Mass: Mass{value: decimal.NewFromFloat(58.44)},
			Volume: Volume{value: decimal.NewFromInt(500), unit: none}},
			molarity: decimal.NewFromInt(2),
			massForMoles: Mass{value: decimal.NewFromFloat(58.44), prefix: kilo, unit: gram},
			expectedMoles: decimal.NewFromFloat(1000),
			expectedError: false,
	},
	{
		compound: Compound{Symbol:  "HOH", // This is an alternative way to symbolize water
			Elements: []ElementMoles{
			{Element: Element{Symbol: "H", Name: "Hydrogen", AtomicNumber: 1, AtomicWeight: decimal.NewFromFloat(1.008)}, Moles: decimal.NewFromFloat(2)},
			{Element: Element{Symbol: "O", Name: "Oxygen", AtomicNumber: 8, AtomicWeight: decimal.NewFromFloat(15.999)}, Moles: decimal.NewFromFloat(1)},
		},
		Mass: Mass{value: decimal.NewFromFloat(18.015)},
		Volume: Volume{value: decimal.NewFromFloat(1), unit: micro},},
		molarity: preciseHOHMoles.Mul(decimal.NewFromInt(1000000)),
		massForMoles: Mass{value: decimal.NewFromFloat(50), prefix: none, unit: pound},
		expectedMoles: preciseHOHMoles,
		expectedError: false,
	},
	{
		compound: Compound{
			Symbol:  "C6H12O6",
			Elements: []ElementMoles{
				{Element: Element{Symbol: "C", Name: "Carbon", AtomicNumber: 6, AtomicWeight: decimal.NewFromFloat(12.011)}, Moles: decimal.NewFromFloat(6)},
				{Element: Element{Symbol: "H", Name: "Hydrogen", AtomicNumber: 1, AtomicWeight: decimal.NewFromFloat(1.008)}, Moles: decimal.NewFromFloat(12)},
				{Element: Element{Symbol: "O", Name: "Oxygen", AtomicNumber: 8, AtomicWeight: decimal.NewFromFloat(15.999)}, Moles: decimal.NewFromFloat(6)},
			},
		Mass: Mass{value: decimal.NewFromFloat(180.156)},
		Volume: Volume{value: decimal.NewFromFloat(2.50), unit: none}},
		molarity: decimal.NewFromFloat(0.4),
		massForMoles: Mass{value: decimal.NewFromFloat(18.0156),prefix: deca, unit: gram},
		expectedMoles: decimal.NewFromFloat(1),
		expectedError: false,
	},
	{
		compound: Compound{
			Symbol: "XYZ",
			Elements: nil},
		massForMoles: Mass{value: decimal.NewFromFloat(5),prefix: kilo, unit: gram},
		expectedMoles: decimal.Zero,
		expectedError: true,
	},
	{
		compound: Compound{Symbol: "H2O1X",
		Elements: nil},
		massForMoles: Mass{value: decimal.NewFromFloat(15), prefix: kilo, unit: gram},
		expectedMoles: decimal.Zero,
		expectedError: true,
	},
	{
		compound: Compound{Symbol: "",
		Elements: nil},
		massForMoles: Mass{value: decimal.NewFromFloat(1), prefix: kilo, unit: gram},
		expectedMoles: decimal.Zero,
		expectedError: true,
	},
	{
		compound: Compound{
			Symbol: "HHO",
			Elements: []ElementMoles{
				{Element: Element{Symbol: "H", Name: "Hydrogen", AtomicNumber: 1, AtomicWeight: decimal.NewFromFloat(1.008)}, Moles: decimal.NewFromFloat(2)},
				{Element: Element{Symbol: "O", Name: "Oxygen", AtomicNumber: 8, AtomicWeight: decimal.NewFromFloat(15.999)}, Moles: decimal.NewFromFloat(1)},
			},
			Mass: Mass{value: decimal.NewFromFloat(18.015)},
			Volume: Volume{value: decimal.NewFromInt(1), unit: micro}},
			molarity: decimal.NewFromInt(1),
			massForMoles: Mass{value: decimal.NewFromFloat(18.015), unit: gram, prefix: micro},
			expectedMoles: decimal.NewFromFloat(.000001),
			expectedError: false,

	},
}
