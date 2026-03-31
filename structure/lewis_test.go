package structure

import (
	"strings"
	"testing"
)

// ── extractCharge ─────────────────────────────────────────────────────────────

func TestExtractCharge(t *testing.T) {
	// extractCharge reads trailing digits as the charge magnitude, so inputs
	// whose formula itself ends in a digit (e.g. "NH4+") are ambiguous — the 4
	// would be consumed as magnitude. Only test unambiguous inputs here.
	cases := []struct {
		input       string
		wantFormula string
		wantCharge  int
		wantErr     bool
	}{
		{"H2O", "H2O", 0, false},
		{"OH-", "OH", -1, false},
		{"Al3+", "Al", 3, false},
		{"Fe2+", "Fe", 2, false},
		{"Fe3+", "Fe", 3, false},
		{"Cl-", "Cl", -1, false},
		{"", "", 0, true},
		{"-", "", 0, true},
		{"+", "", 0, true},
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			formula, charge, err := extractCharge(c.input)
			if c.wantErr {
				if err == nil {
					t.Errorf("expected error, got formula=%q charge=%d", formula, charge)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if formula != c.wantFormula {
				t.Errorf("formula: got %q, want %q", formula, c.wantFormula)
			}
			if charge != c.wantCharge {
				t.Errorf("charge: got %d, want %d", charge, c.wantCharge)
			}
		})
	}
}

// ── parseFlat ─────────────────────────────────────────────────────────────────

func TestParseFlat(t *testing.T) {
	cases := []struct {
		input   string
		want    []atomCount
		wantErr bool
	}{
		{"H2O", []atomCount{{"H", 2}, {"O", 1}}, false},
		{"NaCl", []atomCount{{"Na", 1}, {"Cl", 1}}, false},
		{"C6H12O6", []atomCount{{"C", 6}, {"H", 12}, {"O", 6}}, false},
		{"Fe2O3", []atomCount{{"Fe", 2}, {"O", 3}}, false},
		{"H", []atomCount{{"H", 1}}, false},
		// invalid inputs
		{"1H2O", nil, true},
		{"H2x", nil, true},
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			got, err := parseFlat(c.input)
			if c.wantErr {
				if err == nil {
					t.Errorf("expected error for %q but got none", c.input)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if len(got) != len(c.want) {
				t.Fatalf("got %d atoms, want %d: %v", len(got), len(c.want), got)
			}
			for i, want := range c.want {
				if got[i].Symbol != want.Symbol || got[i].Count != want.Count {
					t.Errorf("[%d] got %v, want %v", i, got[i], want)
				}
			}
		})
	}
}

// ── parseGroups ───────────────────────────────────────────────────────────────

func TestParseGroups(t *testing.T) {
	cases := []struct {
		input   string
		wantLen int
		wantErr bool
	}{
		{"Fe(OH)2", 2, false},  // Fe + (OH)×2
		{"Ca(OH)2", 2, false},  // Ca + (OH)×2
		{"Mn(NO3)2", 2, false}, // Mn + (NO3)×2
		{"FeCl3", 2, false},    // Fe + Cl×3 (no parens, flat)
		{"(OH", 0, true},       // unmatched paren
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			got, err := parseGroups(c.input)
			if c.wantErr {
				if err == nil {
					t.Errorf("expected error for %q but got none", c.input)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if len(got) != c.wantLen {
				t.Errorf("got %d groups, want %d: %v", len(got), c.wantLen, got)
			}
		})
	}
}

// ── LookupLewisWithError — main-group molecules ───────────────────────────────

func TestLewis_Water(t *testing.T) {
	ls, err := LookupLewisWithError("H2O")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if ls.Geometry != "bent" {
		t.Errorf("geometry: got %q, want %q", ls.Geometry, "bent")
	}
	if len(ls.Atoms) != 3 {
		t.Errorf("atom count: got %d, want 3", len(ls.Atoms))
	}
	if len(ls.Bonds) != 2 {
		t.Errorf("bond count: got %d, want 2", len(ls.Bonds))
	}
	// Oxygen should have 2 lone pairs
	o := findAtom(ls, "O")
	if o == nil {
		t.Fatal("O atom not found")
	}
	if o.LonePairs != 2 {
		t.Errorf("O lone pairs: got %d, want 2", o.LonePairs)
	}
	if ls.TotalValenceElectrons != 8 {
		t.Errorf("total VE: got %d, want 8", ls.TotalValenceElectrons)
	}
}

func TestLewis_CO2(t *testing.T) {
	ls, err := LookupLewisWithError("CO2")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if ls.Geometry != "linear" {
		t.Errorf("geometry: got %q, want %q", ls.Geometry, "linear")
	}
	// Both C=O bonds should be double bonds
	for _, b := range ls.Bonds {
		if b.Order != 2 {
			t.Errorf("bond %s→%s: got order %d, want 2", b.From, b.To, b.Order)
		}
	}
	c := findAtom(ls, "C")
	if c == nil {
		t.Fatal("C atom not found")
	}
	if c.LonePairs != 0 {
		t.Errorf("C lone pairs: got %d, want 0", c.LonePairs)
	}
}

func TestLewis_NH3(t *testing.T) {
	ls, err := LookupLewisWithError("NH3")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if ls.Geometry != "trigonal_pyramidal" {
		t.Errorf("geometry: got %q, want %q", ls.Geometry, "trigonal_pyramidal")
	}
	n := findAtom(ls, "N")
	if n == nil {
		t.Fatal("N atom not found")
	}
	if n.LonePairs != 1 {
		t.Errorf("N lone pairs: got %d, want 1", n.LonePairs)
	}
}

func TestLewis_HCl(t *testing.T) {
	ls, err := LookupLewisWithError("HCl")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(ls.Atoms) != 2 {
		t.Errorf("atom count: got %d, want 2", len(ls.Atoms))
	}
	if len(ls.Bonds) != 1 {
		t.Errorf("bond count: got %d, want 1", len(ls.Bonds))
	}
	if ls.Bonds[0].Order != 1 {
		t.Errorf("bond order: got %d, want 1", ls.Bonds[0].Order)
	}
}

func TestLewis_H2(t *testing.T) {
	ls, err := LookupLewisWithError("H2")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(ls.Bonds) != 1 || ls.Bonds[0].Order != 1 {
		t.Errorf("expected single bond, got %v", ls.Bonds)
	}
}

func TestLewis_N2(t *testing.T) {
	ls, err := LookupLewisWithError("N2")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(ls.Bonds) != 1 || ls.Bonds[0].Order != 3 {
		t.Errorf("expected triple bond, got order %d", ls.Bonds[0].Order)
	}
}

func TestLewis_CH4(t *testing.T) {
	ls, err := LookupLewisWithError("CH4")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if ls.Geometry != "tetrahedral" {
		t.Errorf("geometry: got %q, want %q", ls.Geometry, "tetrahedral")
	}
	c := findAtom(ls, "C")
	if c == nil {
		t.Fatal("C atom not found")
	}
	if c.LonePairs != 0 {
		t.Errorf("C lone pairs: got %d, want 0", c.LonePairs)
	}
}

// ── LookupLewisWithError — ions ───────────────────────────────────────────────

func TestLewis_HydroxideIon(t *testing.T) {
	ls, err := LookupLewisWithError("OH-")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if ls.Charge != -1 {
		t.Errorf("charge: got %d, want -1", ls.Charge)
	}
	o := findAtom(ls, "O")
	if o == nil {
		t.Fatal("O atom not found")
	}
	if o.LonePairs != 3 {
		t.Errorf("O lone pairs in OH⁻: got %d, want 3", o.LonePairs)
	}
}

func TestLewis_FluorideIon(t *testing.T) {
	ls, err := LookupLewisWithError("F-")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if ls.Charge != -1 {
		t.Errorf("charge: got %d, want -1", ls.Charge)
	}
	f := findAtom(ls, "F")
	if f == nil {
		t.Fatal("F atom not found")
	}
	if f.LonePairs != 4 {
		t.Errorf("F lone pairs: got %d, want 4", f.LonePairs)
	}
}

// ── LookupLewisWithError — transition metal compounds ────────────────────────

func TestLewis_FeCl3(t *testing.T) {
	ls, err := LookupLewisWithError("FeCl3")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	fe := findAtom(ls, "Fe")
	if fe == nil {
		t.Fatal("Fe atom not found")
	}
	if fe.FormalCharge != 3 {
		t.Errorf("Fe formal charge: got %d, want +3", fe.FormalCharge)
	}
	if fe.LonePairs != 0 {
		t.Errorf("Fe lone pairs: got %d, want 0 (ionic model)", fe.LonePairs)
	}
	if ls.Geometry != "trigonal_planar" {
		t.Errorf("geometry: got %q, want trigonal_planar", ls.Geometry)
	}
}

func TestLewis_FeOH2(t *testing.T) {
	ls, err := LookupLewisWithError("Fe(OH)2")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	fe := findAtom(ls, "Fe")
	if fe == nil {
		t.Fatal("Fe atom not found")
	}
	if fe.FormalCharge != 2 {
		t.Errorf("Fe formal charge: got %d, want +2", fe.FormalCharge)
	}
	// Steps should mention d-electron configuration for TMs
	stepsJoined := strings.Join(ls.Steps, " ")
	if !strings.Contains(stepsJoined, "d6") {
		t.Errorf("expected d6 configuration in steps for Fe²⁺, got: %s", stepsJoined)
	}
}

func TestLewis_MnNO32(t *testing.T) {
	ls, err := LookupLewisWithError("Mn(NO3)2")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	mn := findAtom(ls, "Mn")
	if mn == nil {
		t.Fatal("Mn atom not found")
	}
	if mn.FormalCharge != 2 {
		t.Errorf("Mn formal charge: got %d, want +2", mn.FormalCharge)
	}
}

// ── LookupLewisWithError — main-group metal compounds ────────────────────────

func TestLewis_CaOH2(t *testing.T) {
	ls, err := LookupLewisWithError("Ca(OH)2")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	ca := findAtom(ls, "Ca")
	if ca == nil {
		t.Fatal("Ca atom not found")
	}
	if ca.FormalCharge != 2 {
		t.Errorf("Ca formal charge: got %d, want +2", ca.FormalCharge)
	}
	// Steps should NOT mention d-electron configuration for Ca
	stepsJoined := strings.Join(ls.Steps, " ")
	if strings.Contains(stepsJoined, "electron configuration") {
		t.Errorf("d-electron language should not appear for Ca, got: %s", stepsJoined)
	}
	// Step should say "Identify the metal" not "transition metal"
	if !strings.Contains(stepsJoined, "Identify the metal: Calcium") {
		t.Errorf("expected 'Identify the metal: Calcium' in steps, got: %s", stepsJoined)
	}
}

func TestLewis_NaCl(t *testing.T) {
	ls, err := LookupLewisWithError("NaCl")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	na := findAtom(ls, "Na")
	if na == nil {
		t.Fatal("Na atom not found")
	}
	if na.FormalCharge != 1 {
		t.Errorf("Na formal charge: got %d, want +1", na.FormalCharge)
	}
}

// ── LookupLewisWithError — error cases ───────────────────────────────────────

func TestLewis_Errors(t *testing.T) {
	cases := []struct {
		input       string
		errContains string
	}{
		{"", "empty input"},
		{"Xx", "unknown element"},
		{"H4", "formula contains only hydrogen"},
		{"1H2O", "unexpected character"},
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			_, err := LookupLewisWithError(c.input)
			if err == nil {
				t.Fatalf("expected error containing %q but got none", c.errContains)
			}
			if !strings.Contains(err.Error(), c.errContains) {
				t.Errorf("error %q does not contain %q", err.Error(), c.errContains)
			}
		})
	}
}

// ── helpers ───────────────────────────────────────────────────────────────────

// findAtom returns the first atom in ls with the given element symbol, or nil.
func findAtom(ls *LewisStructure, symbol string) *LewisAtom {
	for i := range ls.Atoms {
		if ls.Atoms[i].Element == symbol {
			return &ls.Atoms[i]
		}
	}
	return nil
}
