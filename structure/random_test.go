package structure

import (
	"fmt"
	"testing"
)

// expectedGeometry derives the molecular geometry from bond count and lone pairs
// using the same lookup table as determineGeometry, so the test is self-consistent
// with the generator's own rules.
func expectedGeometry(n, lpCenter int) string {
	return determineGeometry(n, lpCenter)
}

// lpCenterFor computes the number of lone pairs on the central atom for a given
// (center, terminal, n) combo using the same valence-electron arithmetic as
// buildValidCombos.
func lpCenterFor(center, terminal string, n int) int {
	vc  := elemValence(center)
	vt  := elemValence(terminal)
	tlp := 6
	if terminal == "H" {
		tlp = 0
	}
	totalVE  := vc + n*vt
	onCenter := totalVE - 2*n - n*tlp
	return onCenter / 2
}

// TestValidCombos_AllGenerate checks that every combo produced by buildValidCombos
// successfully generates a Lewis structure without error.
func TestValidCombos_AllGenerate(t *testing.T) {
	for _, c := range validCombos {
		c := c
		formula := fmt.Sprintf("%s%s%d", c.Center, c.Terminal, c.N)
		t.Run(formula, func(t *testing.T) {
			t.Parallel()
			ls, err := generateFromParsed(formula, 0)
			if err != nil {
				t.Fatalf("generateFromParsed(%q) error: %v", formula, err)
			}
			if ls.Geometry == "unknown" {
				t.Errorf("geometry is \"unknown\" for %s — determineGeometry table may be missing a case", formula)
			}
		})
	}
}

// TestValidCombos_GeometryMatchesPrediction checks that the geometry returned by
// the generator matches what the valence-electron arithmetic predicts for every
// combo in validCombos.
func TestValidCombos_GeometryMatchesPrediction(t *testing.T) {
	for _, c := range validCombos {
		c := c
		formula := fmt.Sprintf("%s%s%d", c.Center, c.Terminal, c.N)
		t.Run(formula, func(t *testing.T) {
			t.Parallel()
			ls, err := generateFromParsed(formula, 0)
			if err != nil {
				t.Fatalf("generateFromParsed(%q) error: %v", formula, err)
			}
			lp   := lpCenterFor(c.Center, c.Terminal, c.N)
			want := expectedGeometry(c.N, lp)
			if ls.Geometry != want {
				t.Errorf("%s: geometry = %q, want %q (n=%d, lpCenter=%d)",
					formula, ls.Geometry, want, c.N, lp)
			}
		})
	}
}

// TestLinearCombos verifies that every combo the electron-pair math predicts
// to be linear actually comes back with geometry "linear".
func TestLinearCombos(t *testing.T) {
	for _, c := range validCombos {
		lp := lpCenterFor(c.Center, c.Terminal, c.N)
		// Linear occurs at {n=2, lp=0} (e.g. BeF2) and {n=2, lp=3} (e.g. XeF2)
		if !(c.N == 2 && (lp == 0 || lp == 3)) {
			continue
		}
		c := c
		formula := fmt.Sprintf("%s%s%d", c.Center, c.Terminal, c.N)
		t.Run(formula, func(t *testing.T) {
			t.Parallel()
			ls, err := generateFromParsed(formula, 0)
			if err != nil {
				t.Fatalf("generateFromParsed(%q) error: %v", formula, err)
			}
			if ls.Geometry != "linear" {
				t.Errorf("%s: geometry = %q, want \"linear\" (n=%d, lpCenter=%d)",
					formula, ls.Geometry, c.N, lp)
			}
		})
	}
}

// TestGenerateRandom checks that GenerateRandom returns a valid non-nil structure
// with a recognised geometry across multiple calls.
func TestGenerateRandom(t *testing.T) {
	known := map[string]bool{
		"linear": true, "bent": true, "trigonal_planar": true,
		"trigonal_pyramidal": true, "tetrahedral": true,
		"seesaw": true, "t_shaped": true, "square_planar": true,
		"trigonal_bipyramidal": true, "square_pyramidal": true,
		"octahedral": true,
	}
	for i := range 30 {
		ls, err := GenerateRandom()
		if err != nil {
			t.Fatalf("call %d: GenerateRandom() error: %v", i, err)
		}
		if ls == nil {
			t.Fatalf("call %d: GenerateRandom() returned nil", i)
		}
		if !known[ls.Geometry] {
			t.Errorf("call %d (%s): unrecognised geometry %q", i, ls.Formula, ls.Geometry)
		}
	}
}
