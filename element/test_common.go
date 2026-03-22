package element

import (
	"chemhelper/units"

	"github.com/shopspring/decimal"
)

// preciseHOHMoles is too precise for float64 literals.
var preciseHOHMoles, _ = decimal.NewFromString("1258.9286705523175132")

// Package-level volume helpers — NewVolume calls can't inline into struct
// composite literals, so we define them here.
var (
	vol1L, _    = units.NewVolume(decimal.NewFromInt(1))
	vol500L, _  = units.NewVolume(decimal.NewFromInt(500))
	vol1uL, _   = units.NewVolume(decimal.NewFromFloat(1), units.Micro)
	vol2pt5L, _ = units.NewVolume(decimal.NewFromFloat(2.50))
)

// Package-level mass helpers for Compound.Mass (stores molar mass reference).
var (
	massH2O, _  = units.NewMass(decimal.NewFromFloat(18.015))
	massNaCl, _ = units.NewMass(decimal.NewFromFloat(58.44))
	massGluc, _ = units.NewMass(decimal.NewFromFloat(180.156))
)

var TestCompounds = []struct {
	compound      Compound
	expectedError bool
	massForMoles  units.Mass
	molarity      decimal.Decimal
	expectedMoles decimal.Decimal
}{
	{
		compound: Compound{
			Symbol: "H2O",
			Elements: []ElementMoles{
				{Element: Element{Symbol: "H", Name: "Hydrogen", AtomicNumber: 1, AtomicWeight: decimal.NewFromFloat(1.008)}, Moles: decimal.NewFromFloat(2)},
				{Element: Element{Symbol: "O", Name: "Oxygen", AtomicNumber: 8, AtomicWeight: decimal.NewFromFloat(15.999)}, Moles: decimal.NewFromFloat(1)},
			},
			Mass:   massH2O,
			Volume: vol1L,
		},
		molarity:      decimal.NewFromInt(1),
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(18.015)); return m }(),
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
			Mass:   massNaCl,
			Volume: vol500L,
		},
		molarity:      decimal.NewFromInt(2),
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(58.44), units.Kilo); return m }(),
		expectedMoles: decimal.NewFromFloat(1000),
		expectedError: false,
	},
	{
		compound: Compound{
			Symbol: "HOH",
			Elements: []ElementMoles{
				{Element: Element{Symbol: "H", Name: "Hydrogen", AtomicNumber: 1, AtomicWeight: decimal.NewFromFloat(1.008)}, Moles: decimal.NewFromFloat(2)},
				{Element: Element{Symbol: "O", Name: "Oxygen", AtomicNumber: 8, AtomicWeight: decimal.NewFromFloat(15.999)}, Moles: decimal.NewFromFloat(1)},
			},
			Mass:   massH2O,
			Volume: vol1uL,
		},
		molarity:      preciseHOHMoles.Mul(decimal.NewFromInt(1000000)),
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(50), units.Pound); return m }(),
		expectedMoles: preciseHOHMoles,
		expectedError: false,
	},
	{
		compound: Compound{
			Symbol: "C6H12O6",
			Elements: []ElementMoles{
				{Element: Element{Symbol: "C", Name: "Carbon", AtomicNumber: 6, AtomicWeight: decimal.NewFromFloat(12.011)}, Moles: decimal.NewFromFloat(6)},
				{Element: Element{Symbol: "H", Name: "Hydrogen", AtomicNumber: 1, AtomicWeight: decimal.NewFromFloat(1.008)}, Moles: decimal.NewFromFloat(12)},
				{Element: Element{Symbol: "O", Name: "Oxygen", AtomicNumber: 8, AtomicWeight: decimal.NewFromFloat(15.999)}, Moles: decimal.NewFromFloat(6)},
			},
			Mass:   massGluc,
			Volume: vol2pt5L,
		},
		molarity:      decimal.NewFromFloat(0.4),
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(18.0156), units.Deca); return m }(),
		expectedMoles: decimal.NewFromFloat(1),
		expectedError: false,
	},
	{
		compound:      Compound{Symbol: "XYZ", Elements: nil},
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(5), units.Kilo); return m }(),
		expectedMoles: decimal.Zero,
		expectedError: true,
	},
	{
		compound:      Compound{Symbol: "H2O1X", Elements: nil},
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(15), units.Kilo); return m }(),
		expectedMoles: decimal.Zero,
		expectedError: true,
	},
	{
		compound:      Compound{Symbol: "", Elements: nil},
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(1), units.Kilo); return m }(),
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
			Mass:   massH2O,
			Volume: vol1uL,
		},
		molarity:      decimal.NewFromInt(1),
		massForMoles:  func() units.Mass { m, _ := units.NewMass(decimal.NewFromFloat(18.015), units.Micro); return m }(),
		expectedMoles: decimal.NewFromFloat(.000001),
		expectedError: false,
	},
}
