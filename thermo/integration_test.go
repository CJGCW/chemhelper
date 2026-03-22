// Integration tests for colligative properties. These tests mirror real
// word problems by chaining units → solution → thermo the same way
// application code would.
//
// Inputs are expressed in standard scientific notation with explicit precision
// (e.g. 0.500 kg rather than 500 g) so significant figures are unambiguous.
// Results are rounded to the number of sig figs of the least precise input
// before comparison.
package thermo_test

import (
	"testing"

	"chemhelper/solution"
	"chemhelper/thermo"
	"chemhelper/units"

	"github.com/shopspring/decimal"
)

// sigFigsFrom returns the lowest significant figure count across a set of
// input strings, determining the correct precision for the result.
func sigFigsFrom(t *testing.T, inputs ...string) int32 {
	t.Helper()
	n, err := units.GetLowestSignificantFigures(inputs)
	if err != nil {
		t.Fatalf("GetLowestSignificantFigures: %v", err)
	}
	return int32(n)
}

// toSigFigs rounds val to n significant figures using banker's rounding.
func toSigFigs(t *testing.T, val decimal.Decimal, n int32) decimal.Decimal {
	t.Helper()
	rounded, err := units.SetToSigFigs(val, n)
	if err != nil {
		t.Fatalf("SetToSigFigs: %v", err)
	}
	return rounded
}

// mustMoles derives moles from a units.Mass and a molar mass in g/mol.
func mustMoles(t *testing.T, mass units.Mass, molarMassGPerMol decimal.Decimal) decimal.Decimal {
	t.Helper()
	moles, err := mass.GetMoles(molarMassGPerMol)
	if err != nil {
		t.Fatalf("GetMoles: %v", err)
	}
	return moles
}

// mustMolality derives molality from moles and a units.Mass of solvent.
func mustMolality(t *testing.T, moles decimal.Decimal, solventMass units.Mass) decimal.Decimal {
	t.Helper()
	result, err := solution.FindMolality{
		Moles:       moles,
		SolventMass: solventMass,
	}.Calculate()
	if err != nil {
		t.Fatalf("FindMolality: %v", err)
	}
	return result.Value
}

// TestBPE_NaClInWater mirrors the problem:
// "0.500 kg of NaCl is dissolved in 1.250 kg of water. What is the boiling
// point elevation, and what is the new boiling point?"
//
// NaCl molar mass: 58.44 g/mol (4 sig figs)
// Water Kb: 0.512 degrees C * kg/mol (3 sig figs) <- limiting
// NaCl i = 2 (exact)
// Result precision: 3 sig figs
func TestBPE_NaClInWater(t *testing.T) {
	sf := sigFigsFrom(t, "0.500", "58.44", "1.250", "0.512")

	naclMass, _ := units.NewMass(decimal.NewFromFloat(0.500), units.Kilo)
	moles := mustMoles(t, naclMass, decimal.NewFromFloat(58.44))

	waterMass, _ := units.NewMass(decimal.NewFromFloat(1.250), units.Kilo)
	molality := mustMolality(t, moles, waterMass)

	bpe := thermo.BoilingPointElevation{
		Solvent:        thermo.Water,
		Molality:       molality,
		VantHoffFactor: decimal.NewFromInt(2),
	}

	deltaResult, err := bpe.Calculate()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if deltaResult.Unit != "°C" {
		t.Errorf("expected unit °C, got %s", deltaResult.Unit)
	}
	if len(deltaResult.Steps) == 0 {
		t.Errorf("expected steps but got none")
	}
	expectedDelta := toSigFigs(t, decimal.NewFromFloat(7.01), sf)
	actualDelta := toSigFigs(t, deltaResult.Value, sf)
	if !actualDelta.Equal(expectedDelta) {
		t.Errorf("dTb: expected %v, got %v (%d sig figs)", expectedDelta, actualDelta, sf)
	}

	bpResult, err := bpe.NewBoilingPoint()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expectedBP := toSigFigs(t, decimal.NewFromInt(107), sf)
	actualBP := toSigFigs(t, bpResult.Value, sf)
	if !actualBP.Equal(expectedBP) {
		t.Errorf("new BP: expected %v, got %v (%d sig figs)", expectedBP, actualBP, sf)
	}
}

// TestFPD_GlucoseInWater mirrors the problem:
// "0.200 kg of glucose (C6H12O6, molar mass 180.156 g/mol) is dissolved in
// 0.500 kg of water. What is the freezing point depression and the new
// freezing point?"
//
// Water Kf: 1.86 degrees C * kg/mol (3 sig figs) <- limiting
// Glucose i = 1 (non-electrolyte, exact) - VantHoffFactor left unset to
// verify it defaults to 1
// Result precision: 3 sig figs
func TestFPD_GlucoseInWater(t *testing.T) {
	sf := sigFigsFrom(t, "0.200", "180.156", "0.500", "1.86")

	glucoseMass, _ := units.NewMass(decimal.NewFromFloat(0.200), units.Kilo)
	moles := mustMoles(t, glucoseMass, decimal.NewFromFloat(180.156))

	waterMass, _ := units.NewMass(decimal.NewFromFloat(0.500), units.Kilo)
	molality := mustMolality(t, moles, waterMass)

	fpd := thermo.FreezingPointDepression{
		Solvent:  thermo.Water,
		Molality: molality,
	}

	deltaResult, err := fpd.Calculate()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expectedDelta := toSigFigs(t, decimal.NewFromFloat(4.13), sf)
	actualDelta := toSigFigs(t, deltaResult.Value, sf)
	if !actualDelta.Equal(expectedDelta) {
		t.Errorf("dTf: expected %v, got %v (%d sig figs)", expectedDelta, actualDelta, sf)
	}

	fpResult, err := fpd.NewFreezingPoint()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expectedFP := toSigFigs(t, decimal.NewFromFloat(-4.13), sf)
	actualFP := toSigFigs(t, fpResult.Value, sf)
	if !actualFP.Equal(expectedFP) {
		t.Errorf("new FP: expected %v, got %v (%d sig figs)", expectedFP, actualFP, sf)
	}
}

// TestFPD_CaCl2InWater mirrors the problem:
// "0.100 kg of CaCl2 (molar mass 110.98 g/mol) is dissolved in 0.800 kg of
// water. What is the freezing point depression and the new freezing point?"
//
// CaCl2 -> Ca2+ + 2 Cl-, i = 3 (exact)
// Water Kf: 1.86 degrees C * kg/mol (3 sig figs) <- limiting
// Result precision: 3 sig figs
func TestFPD_CaCl2InWater(t *testing.T) {
	sf := sigFigsFrom(t, "0.100", "110.98", "0.800", "1.86")

	cacl2Mass, _ := units.NewMass(decimal.NewFromFloat(0.100), units.Kilo)
	moles := mustMoles(t, cacl2Mass, decimal.NewFromFloat(110.98))

	waterMass, _ := units.NewMass(decimal.NewFromFloat(0.800), units.Kilo)
	molality := mustMolality(t, moles, waterMass)

	fpd := thermo.FreezingPointDepression{
		Solvent:        thermo.Water,
		Molality:       molality,
		VantHoffFactor: decimal.NewFromInt(3),
	}

	deltaResult, err := fpd.Calculate()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expectedDelta := toSigFigs(t, decimal.NewFromFloat(6.28), sf)
	actualDelta := toSigFigs(t, deltaResult.Value, sf)
	if !actualDelta.Equal(expectedDelta) {
		t.Errorf("dTf: expected %v, got %v (%d sig figs)", expectedDelta, actualDelta, sf)
	}

	fpResult, err := fpd.NewFreezingPoint()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expectedFP := toSigFigs(t, decimal.NewFromFloat(-6.28), sf)
	actualFP := toSigFigs(t, fpResult.Value, sf)
	if !actualFP.Equal(expectedFP) {
		t.Errorf("new FP: expected %v, got %v (%d sig figs)", expectedFP, actualFP, sf)
	}
}

// TestBPE_EthanolSolvent mirrors the problem:
// "0.0500 kg of a non-electrolyte solute (molar mass 250.0 g/mol) is dissolved
// in 0.200 kg of ethanol. What is the boiling point elevation?"
//
// Ethanol Kb: 1.22 degrees C * kg/mol (3 sig figs) <- limiting
// moles = 50/250 = 0.200, molality = 0.200/0.200 = 1.00
// dTb = 1.22 * 1.00 * 1 = 1.22
// Result precision: 3 sig figs
func TestBPE_EthanolSolvent(t *testing.T) {
	sf := sigFigsFrom(t, "0.0500", "250.0", "0.200", "1.22")

	soluteMass, _ := units.NewMass(decimal.NewFromFloat(0.0500), units.Kilo)
	moles := mustMoles(t, soluteMass, decimal.NewFromFloat(250.0))

	solventMass, _ := units.NewMass(decimal.NewFromFloat(0.200), units.Kilo)
	molality := mustMolality(t, moles, solventMass)

	bpe := thermo.BoilingPointElevation{
		Solvent:  thermo.Ethanol,
		Molality: molality,
	}

	deltaResult, err := bpe.Calculate()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expectedDelta := toSigFigs(t, decimal.NewFromFloat(1.22), sf)
	actualDelta := toSigFigs(t, deltaResult.Value, sf)
	if !actualDelta.Equal(expectedDelta) {
		t.Errorf("dTb: expected %v, got %v (%d sig figs)", expectedDelta, actualDelta, sf)
	}
}

// TestFPD_BackCalculateMolality mirrors the problem:
// "The freezing point of a water solution is measured at -3.72 degrees C.
// The solute is NaCl (i=2). What is the molality?"
//
// Water Kf: 1.86 degrees C * kg/mol (3 sig figs) <- limiting
// m = 3.72 / (1.86 * 2) = 1.00 mol/kg
// Result precision: 3 sig figs
func TestFPD_BackCalculateMolality(t *testing.T) {
	sf := sigFigsFrom(t, "3.72", "1.86")

	f := thermo.FindMolalityFromFPD{
		Solvent:        thermo.Water,
		DeltaTf:        decimal.NewFromFloat(3.72),
		VantHoffFactor: decimal.NewFromInt(2),
	}

	result, err := f.Calculate()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if result.Unit != "mol/kg" {
		t.Errorf("expected unit mol/kg, got %s", result.Unit)
	}
	expected := toSigFigs(t, decimal.NewFromInt(1), sf)
	actual := toSigFigs(t, result.Value, sf)
	if !actual.Equal(expected) {
		t.Errorf("molality: expected %v, got %v (%d sig figs)", expected, actual, sf)
	}
}

// TestBPE_BackCalculateMolality mirrors the problem:
// "The boiling point of a benzene solution is measured at 83.91 degrees C.
// The solute is a non-electrolyte. What is the molality?"
//
// Benzene BP: 80.1 degrees C (3 sig figs), Kb: 2.53 degrees C * kg/mol (3 sig figs) <- limiting
// dTb = 83.91 - 80.1 = 3.81 degrees C
// m = 3.81 / 2.53 = 1.506... -> 1.51 at 3 sig figs
// Result precision: 3 sig figs
func TestBPE_BackCalculateMolality(t *testing.T) {
	sf := sigFigsFrom(t, "83.91", "80.1", "2.53")

	deltaTb := decimal.NewFromFloat(83.91).Sub(thermo.Benzene.BoilingPoint)

	f := thermo.FindMolalityFromBPE{
		Solvent: thermo.Benzene,
		DeltaTb: deltaTb,
	}

	result, err := f.Calculate()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expected := toSigFigs(t, decimal.NewFromFloat(3.81).Div(decimal.NewFromFloat(2.53)), sf)
	actual := toSigFigs(t, result.Value, sf)
	if !actual.Equal(expected) {
		t.Errorf("molality: expected %v, got %v (%d sig figs)", expected, actual, sf)
	}
}
