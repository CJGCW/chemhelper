package thermo

import (
	"testing"

	"github.com/shopspring/decimal"
)

// dec is a test helper for cleaner decimal literals.
func dec(f float64) decimal.Decimal {
	return decimal.NewFromFloat(f)
}

func TestBoilingPointElevation(t *testing.T) {
	tests := []struct {
		name           string
		solvent        Solvent
		molality       float64
		vantHoff       float64 // 0 means unset, defaults to 1
		expectedDelta  float64
		expectedError  bool
	}{
		{
			// Water, 1 mol/kg, non-electrolyte: ΔTb = 0.512 × 1 × 1 = 0.512
			name:          "water, 1 mol/kg non-electrolyte",
			solvent:       Water,
			molality:      1,
			expectedDelta: 0.512,
			expectedError: false,
		},
		{
			// Water, 1 mol/kg NaCl (i=2): ΔTb = 0.512 × 1 × 2 = 1.024
			name:          "water, 1 mol/kg NaCl (i=2)",
			solvent:       Water,
			molality:      1,
			vantHoff:      2,
			expectedDelta: 1.024,
			expectedError: false,
		},
		{
			// Benzene, 0.5 mol/kg: ΔTb = 2.53 × 0.5 × 1 = 1.265
			name:          "benzene, 0.5 mol/kg non-electrolyte",
			solvent:       Benzene,
			molality:      0.5,
			expectedDelta: 1.265,
			expectedError: false,
		},
		{
			// Water, 2 mol/kg CaCl2 (i=3): ΔTb = 0.512 × 2 × 3 = 3.072
			name:          "water, 2 mol/kg CaCl2 (i=3)",
			solvent:       Water,
			molality:      2,
			vantHoff:      3,
			expectedDelta: 3.072,
			expectedError: false,
		},
		{
			name:          "zero molality returns error",
			solvent:       Water,
			molality:      0,
			expectedError: true,
		},
		{
			name:          "negative van't Hoff factor returns error",
			solvent:       Water,
			molality:      1,
			vantHoff:      -1,
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			b := BoilingPointElevation{
				Solvent:        test.solvent,
				Molality:       dec(test.molality),
				VantHoffFactor: dec(test.vantHoff),
			}
			result, err := b.Calculate()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(dec(test.expectedDelta)) {
				t.Errorf("expected ΔTb = %v, got %v", test.expectedDelta, result.Value)
			}
			if result.Unit != "°C" {
				t.Errorf("expected unit °C, got %s", result.Unit)
			}
			if len(result.Steps) == 0 {
				t.Errorf("expected steps but got none")
			}
		})
	}
}

func TestNewBoilingPoint(t *testing.T) {
	tests := []struct {
		name          string
		solvent       Solvent
		molality      float64
		vantHoff      float64
		expectedBP    float64
		expectedError bool
	}{
		{
			// Water, 1 mol/kg: new bp = 100 + 0.512 = 100.512
			name:       "water, 1 mol/kg",
			solvent:    Water,
			molality:   1,
			expectedBP: 100.512,
		},
		{
			// Benzene, 1 mol/kg: new bp = 80.1 + 2.53 = 82.63
			name:       "benzene, 1 mol/kg",
			solvent:    Benzene,
			molality:   1,
			expectedBP: 82.63,
		},
		{
			// Water, 1 mol/kg NaCl (i=2): new bp = 100 + 1.024 = 101.024
			name:       "water, 1 mol/kg NaCl (i=2)",
			solvent:    Water,
			molality:   1,
			vantHoff:   2,
			expectedBP: 101.024,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			b := BoilingPointElevation{
				Solvent:        test.solvent,
				Molality:       dec(test.molality),
				VantHoffFactor: dec(test.vantHoff),
			}
			result, err := b.NewBoilingPoint()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(dec(test.expectedBP)) {
				t.Errorf("expected new bp = %v, got %v", test.expectedBP, result.Value)
			}
		})
	}
}

func TestFreezingPointDepression(t *testing.T) {
	tests := []struct {
		name          string
		solvent       Solvent
		molality      float64
		vantHoff      float64
		expectedDelta float64
		expectedError bool
	}{
		{
			// Water, 1 mol/kg, non-electrolyte: ΔTf = 1.86 × 1 × 1 = 1.86
			name:          "water, 1 mol/kg non-electrolyte",
			solvent:       Water,
			molality:      1,
			expectedDelta: 1.86,
			expectedError: false,
		},
		{
			// Water, 1 mol/kg NaCl (i=2): ΔTf = 1.86 × 1 × 2 = 3.72
			name:          "water, 1 mol/kg NaCl (i=2)",
			solvent:       Water,
			molality:      1,
			vantHoff:      2,
			expectedDelta: 3.72,
			expectedError: false,
		},
		{
			// Cyclohexane, 0.5 mol/kg: ΔTf = 20.2 × 0.5 × 1 = 10.1
			name:          "cyclohexane, 0.5 mol/kg",
			solvent:       Cyclohexane,
			molality:      0.5,
			expectedDelta: 10.1,
			expectedError: false,
		},
		{
			// Water, 2 mol/kg CaCl2 (i=3): ΔTf = 1.86 × 2 × 3 = 11.16
			name:          "water, 2 mol/kg CaCl2 (i=3)",
			solvent:       Water,
			molality:      2,
			vantHoff:      3,
			expectedDelta: 11.16,
			expectedError: false,
		},
		{
			name:          "zero molality returns error",
			solvent:       Water,
			molality:      0,
			expectedError: true,
		},
		{
			name:          "negative van't Hoff factor returns error",
			solvent:       Water,
			molality:      1,
			vantHoff:      -1,
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FreezingPointDepression{
				Solvent:        test.solvent,
				Molality:       dec(test.molality),
				VantHoffFactor: dec(test.vantHoff),
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
			if !result.Value.Equal(dec(test.expectedDelta)) {
				t.Errorf("expected ΔTf = %v, got %v", test.expectedDelta, result.Value)
			}
			if result.Unit != "°C" {
				t.Errorf("expected unit °C, got %s", result.Unit)
			}
			if len(result.Steps) == 0 {
				t.Errorf("expected steps but got none")
			}
		})
	}
}

func TestNewFreezingPoint(t *testing.T) {
	tests := []struct {
		name          string
		solvent       Solvent
		molality      float64
		vantHoff      float64
		expectedFP    float64
		expectedError bool
	}{
		{
			// Water, 1 mol/kg: new fp = 0 - 1.86 = -1.86
			name:       "water, 1 mol/kg",
			solvent:    Water,
			molality:   1,
			expectedFP: -1.86,
		},
		{
			// Water, 1 mol/kg NaCl (i=2): new fp = 0 - 3.72 = -3.72
			name:       "water, 1 mol/kg NaCl (i=2)",
			solvent:    Water,
			molality:   1,
			vantHoff:   2,
			expectedFP: -3.72,
		},
		{
			// Benzene, 1 mol/kg: new fp = 5.5 - 5.12 = 0.38
			name:       "benzene, 1 mol/kg",
			solvent:    Benzene,
			molality:   1,
			expectedFP: 0.38,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FreezingPointDepression{
				Solvent:        test.solvent,
				Molality:       dec(test.molality),
				VantHoffFactor: dec(test.vantHoff),
			}
			result, err := f.NewFreezingPoint()
			if test.expectedError {
				if err == nil {
					t.Errorf("expected error but got none")
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !result.Value.Equal(dec(test.expectedFP)) {
				t.Errorf("expected new fp = %v, got %v", test.expectedFP, result.Value)
			}
		})
	}
}

func TestFindMolalityFromBPE(t *testing.T) {
	tests := []struct {
		name             string
		solvent          Solvent
		deltaTb          float64
		vantHoff         float64
		expectedMolality float64
		expectedError    bool
	}{
		{
			// Water: m = 0.512 / (0.512 × 1) = 1
			name:             "water, ΔTb=0.512, i=1",
			solvent:          Water,
			deltaTb:          0.512,
			expectedMolality: 1,
			expectedError:    false,
		},
		{
			// Water: m = 1.024 / (0.512 × 2) = 1
			name:             "water, ΔTb=1.024, i=2 (NaCl)",
			solvent:          Water,
			deltaTb:          1.024,
			vantHoff:         2,
			expectedMolality: 1,
			expectedError:    false,
		},
		{
			// Benzene: m = 2.53 / (2.53 × 1) = 1
			name:             "benzene, ΔTb=2.53, i=1",
			solvent:          Benzene,
			deltaTb:          2.53,
			expectedMolality: 1,
			expectedError:    false,
		},
		{
			name:          "zero ΔTb returns error",
			solvent:       Water,
			deltaTb:       0,
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FindMolalityFromBPE{
				Solvent:        test.solvent,
				DeltaTb:        dec(test.deltaTb),
				VantHoffFactor: dec(test.vantHoff),
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
			if !result.Value.Equal(dec(test.expectedMolality)) {
				t.Errorf("expected molality = %v, got %v", test.expectedMolality, result.Value)
			}
			if result.Unit != "mol/kg" {
				t.Errorf("expected unit mol/kg, got %s", result.Unit)
			}
		})
	}
}

func TestFindMolalityFromFPD(t *testing.T) {
	tests := []struct {
		name             string
		solvent          Solvent
		deltaTf          float64
		vantHoff         float64
		expectedMolality float64
		expectedError    bool
	}{
		{
			// Water: m = 1.86 / (1.86 × 1) = 1
			name:             "water, ΔTf=1.86, i=1",
			solvent:          Water,
			deltaTf:          1.86,
			expectedMolality: 1,
			expectedError:    false,
		},
		{
			// Water: m = 3.72 / (1.86 × 2) = 1
			name:             "water, ΔTf=3.72, i=2 (NaCl)",
			solvent:          Water,
			deltaTf:          3.72,
			vantHoff:         2,
			expectedMolality: 1,
			expectedError:    false,
		},
		{
			// Cyclohexane: m = 20.2 / (20.2 × 1) = 1
			name:             "cyclohexane, ΔTf=20.2, i=1",
			solvent:          Cyclohexane,
			deltaTf:          20.2,
			expectedMolality: 1,
			expectedError:    false,
		},
		{
			name:          "zero ΔTf returns error",
			solvent:       Water,
			deltaTf:       0,
			expectedError: true,
		},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			f := FindMolalityFromFPD{
				Solvent:        test.solvent,
				DeltaTf:        dec(test.deltaTf),
				VantHoffFactor: dec(test.vantHoff),
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
			if !result.Value.Equal(dec(test.expectedMolality)) {
				t.Errorf("expected molality = %v, got %v", test.expectedMolality, result.Value)
			}
			if result.Unit != "mol/kg" {
				t.Errorf("expected unit mol/kg, got %s", result.Unit)
			}
		})
	}
}
