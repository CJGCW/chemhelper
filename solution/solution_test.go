package solution

import (
	"testing"

	"chemhelper/units"

	"github.com/shopspring/decimal"
)

// newVol is a test helper to construct a Volume without boilerplate error checking.
func newVol(val float64, opts ...interface{}) units.Volume {
	v, err := units.NewVolume(decimal.NewFromFloat(val), opts...)
	if err != nil {
		panic(err)
	}
	return v
}

// newMass is a test helper to construct a Mass without boilerplate error checking.
func newMass(val float64, opts ...interface{}) units.Mass {
	m, err := units.NewMass(decimal.NewFromFloat(val), opts...)
	if err != nil {
		panic(err)
	}
	return m
}

func TestFindMolarity(t *testing.T) {
	tests := []struct {
		name          string
		moles         float64
		volume        units.Volume
		expectedValue float64
		expectedError bool
	}{
		{
			name:          "1 mol in 1 L",
			moles:         1,
			volume:        newVol(1),
			expectedValue: 1,
			expectedError: false,
		},
		{
			name:          "2 mol in 500 mL",
			moles:         2,
			volume:        newVol(500, units.Milli),
			expectedValue: 4,
			expectedError: false,
		},
		{
			name:          "0.5 mol in 250 mL",
			moles:         0.5,
			volume:        newVol(250, units.Milli),
			expectedValue: 2,
			expectedError: false,
		},
		{
			name:          "zero moles returns error",
			moles:         0,
			volume:        newVol(1),
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FindMolarity{
				Moles:  decimal.NewFromFloat(test.moles),
				Volume: test.volume,
			}
			result, err := f.Calculate()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(decimal.NewFromFloat(test.expectedValue)) {
				t.Errorf("expected %v mol/L, got %v", test.expectedValue, result.Value)
			}
			if result.Unit != "mol/L" {
				t.Errorf("expected unit mol/L, got %s", result.Unit)
			}
			if len(result.Steps) == 0 {
				t.Errorf("expected steps but got none")
			}
		})
	}
}

func TestFindMolarityFromMass(t *testing.T) {
	tests := []struct {
		name          string
		mass          units.Mass
		molarMass     float64
		volume        units.Volume
		expectedValue float64
		expectedError bool
	}{
		{
			// NaCl: 58.44 g/mol; 58.44 g in 1 L → 1 mol/L
			name:          "58.44g NaCl in 1L",
			mass:          newMass(58.44),
			molarMass:     58.44,
			volume:        newVol(1),
			expectedValue: 1,
			expectedError: false,
		},
		{
			// 180.156g glucose (C6H12O6) in 500 mL → 2 mol/L
			name:          "180.156g glucose in 500 mL",
			mass:          newMass(180.156),
			molarMass:     180.156,
			volume:        newVol(500, units.Milli),
			expectedValue: 2,
			expectedError: false,
		},
		{
			name:          "zero molar mass returns error",
			mass:          newMass(10),
			molarMass:     0,
			volume:        newVol(1),
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FindMolarityFromMass{
				Mass:      test.mass,
				MolarMass: decimal.NewFromFloat(test.molarMass),
				Volume:    test.volume,
			}
			result, err := f.Calculate()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(decimal.NewFromFloat(test.expectedValue)) {
				t.Errorf("expected %v mol/L, got %v", test.expectedValue, result.Value)
			}
		})
	}
}

func TestFindMolesFromMolarity(t *testing.T) {
	tests := []struct {
		name          string
		molarity      float64
		volume        units.Volume
		expectedValue float64
		expectedError bool
	}{
		{
			name:          "1 mol/L × 1 L",
			molarity:      1,
			volume:        newVol(1),
			expectedValue: 1,
			expectedError: false,
		},
		{
			name:          "2 mol/L × 250 mL",
			molarity:      2,
			volume:        newVol(250, units.Milli),
			expectedValue: 0.5,
			expectedError: false,
		},
		{
			name:          "0.1 mol/L × 100 mL",
			molarity:      0.1,
			volume:        newVol(100, units.Milli),
			expectedValue: 0.01,
			expectedError: false,
		},
		{
			name:          "zero molarity returns error",
			molarity:      0,
			volume:        newVol(1),
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FindMolesFromMolarity{
				Molarity: decimal.NewFromFloat(test.molarity),
				Volume:   test.volume,
			}
			result, err := f.Calculate()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(decimal.NewFromFloat(test.expectedValue)) {
				t.Errorf("expected %v mol, got %v", test.expectedValue, result.Value)
			}
			if result.Unit != "mol" {
				t.Errorf("expected unit mol, got %s", result.Unit)
			}
		})
	}
}

func TestFindMolality(t *testing.T) {
	tests := []struct {
		name          string
		moles         float64
		solventMass   units.Mass
		expectedValue float64
		expectedError bool
	}{
		{
			// 1 mol in 1 kg → 1 mol/kg
			name:          "1 mol in 1 kg water",
			moles:         1,
			solventMass:   newMass(1000),
			expectedValue: 1,
			expectedError: false,
		},
		{
			// 0.5 mol in 500 g → 1 mol/kg
			name:          "0.5 mol in 500 g solvent",
			moles:         0.5,
			solventMass:   newMass(500),
			expectedValue: 1,
			expectedError: false,
		},
		{
			// 2 mol in 1 kg → 2 mol/kg
			name:          "2 mol in 1 kg solvent",
			moles:         2,
			solventMass:   newMass(1, units.Kilo),
			expectedValue: 2,
			expectedError: false,
		},
		{
			name:          "zero moles returns error",
			moles:         0,
			solventMass:   newMass(1000),
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FindMolality{
				Moles:       decimal.NewFromFloat(test.moles),
				SolventMass: test.solventMass,
			}
			result, err := f.Calculate()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(decimal.NewFromFloat(test.expectedValue)) {
				t.Errorf("expected %v mol/kg, got %v", test.expectedValue, result.Value)
			}
			if result.Unit != "mol/kg" {
				t.Errorf("expected unit mol/kg, got %s", result.Unit)
			}
		})
	}
}

func TestDilutionFindFinalConcentration(t *testing.T) {
	tests := []struct {
		name                 string
		initialConcentration float64
		initialVolume        units.Volume
		finalVolume          units.Volume
		expectedValue        float64
		expectedError        bool
	}{
		{
			// 1 mol/L, 100 mL → dilute to 1 L → 0.1 mol/L
			name:                 "dilute 100 mL of 1M to 1 L",
			initialConcentration: 1,
			initialVolume:        newVol(100, units.Milli),
			finalVolume:          newVol(1),
			expectedValue:        0.1,
			expectedError:        false,
		},
		{
			// 6 mol/L, 50 mL → dilute to 300 mL → 1 mol/L
			name:                 "dilute 50 mL of 6M to 300 mL",
			initialConcentration: 6,
			initialVolume:        newVol(50, units.Milli),
			finalVolume:          newVol(300, units.Milli),
			expectedValue:        1,
			expectedError:        false,
		},
		{
			name:                 "final volume less than initial returns error",
			initialConcentration: 1,
			initialVolume:        newVol(1),
			finalVolume:          newVol(500, units.Milli),
			expectedError:        true,
		},
		{
			name:                 "zero concentration returns error",
			initialConcentration: 0,
			initialVolume:        newVol(1),
			finalVolume:          newVol(2),
			expectedError:        true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			d := DilutionFindFinalConcentration{
				InitialConcentration: decimal.NewFromFloat(test.initialConcentration),
				InitialVolume:        test.initialVolume,
				FinalVolume:          test.finalVolume,
			}
			result, err := d.Calculate()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(decimal.NewFromFloat(test.expectedValue)) {
				t.Errorf("expected %v mol/L, got %v", test.expectedValue, result.Value)
			}
		})
	}
}

func TestDilutionFindFinalVolume(t *testing.T) {
	tests := []struct {
		name                 string
		initialConcentration float64
		initialVolume        units.Volume
		finalConcentration   float64
		expectedValue        float64
		expectedError        bool
	}{
		{
			// 1 mol/L, 100 mL → dilute to 0.1 mol/L → 1 L
			name:                 "dilute 100 mL of 1M to 0.1M",
			initialConcentration: 1,
			initialVolume:        newVol(100, units.Milli),
			finalConcentration:   0.1,
			expectedValue:        1,
			expectedError:        false,
		},
		{
			// 6 mol/L, 50 mL → dilute to 1 mol/L → 300 mL = 0.3 L
			name:                 "dilute 50 mL of 6M to 1M",
			initialConcentration: 6,
			initialVolume:        newVol(50, units.Milli),
			finalConcentration:   1,
			expectedValue:        0.3,
			expectedError:        false,
		},
		{
			name:                 "final concentration >= initial returns error",
			initialConcentration: 1,
			initialVolume:        newVol(1),
			finalConcentration:   2,
			expectedError:        true,
		},
		{
			name:                 "zero final concentration returns error",
			initialConcentration: 1,
			initialVolume:        newVol(1),
			finalConcentration:   0,
			expectedError:        true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			d := DilutionFindFinalVolume{
				InitialConcentration: decimal.NewFromFloat(test.initialConcentration),
				InitialVolume:        test.initialVolume,
				FinalConcentration:   decimal.NewFromFloat(test.finalConcentration),
			}
			result, err := d.Calculate()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(decimal.NewFromFloat(test.expectedValue)) {
				t.Errorf("expected %v L, got %v", test.expectedValue, result.Value)
			}
			if result.Unit != "L" {
				t.Errorf("expected unit L, got %s", result.Unit)
			}
		})
	}
}
