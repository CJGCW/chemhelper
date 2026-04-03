package smiles_test

import (
	"testing"

	"chemhelper/smiles"
)

// ── DetectInputType ───────────────────────────────────────────────────────────

func TestDetectInputType(t *testing.T) {
	cases := []struct {
		input string
		want  string
	}{
		// CID — plain integers
		{"962", "cid"},
		{"5234", "cid"},
		{"2244", "cid"},

		// InChI — InChI= prefix (any casing)
		{"InChI=1S/H2O/h1H2", "inchi"},
		{"INCHI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "inchi"},

		// InChIKey — 27-char uppercase with two hyphens
		{"XLYOFNOQVPJJNP-UHFFFAOYSA-N", "inchikey"},
		{"RYYVLZVUVIJVGH-UHFFFAOYSA-N", "inchikey"}, // caffeine

		// SMILES with special chars
		{"CC(=O)O", "smiles"},       // acetic acid
		{"O=C=O", "smiles"},         // CO2
		{"[Na+].[Cl-]", "smiles"},   // NaCl ionic
		{"C(=O)(=O)", "smiles"},     // CO2 branched
		{"C#N", "smiles"},           // hydrogen cyanide
		{"C/C=C/C", "smiles"},       // trans-2-butene

		// Plain organic SMILES (no special chars, no digits)
		{"CCO", "smiles"},  // ethanol
		{"CC", "smiles"},   // ethane
		{"C", "smiles"},    // methane
		{"N", "smiles"},    // ammonia
		{"O", "smiles"},    // water
		{"CCCl", "smiles"}, // 1-chloropropane
		{"CBr", "smiles"},  // bromomethane
		{"CS", "smiles"},   // methanethiol
		{"CCCO", "smiles"}, // 1-propanol

		// Aromatic SMILES with ring-closure digits
		{"c1ccccc1", "smiles"},   // benzene
		{"c1ccncc1", "smiles"},   // pyridine
		{"c1ccoc1", "smiles"},    // furan

		// Names / IUPAC names
		{"water", "name"},
		{"ethanol", "name"},
		{"caffeine", "name"},
		{"acetic acid", "name"},
		{"sodium chloride", "name"},

		// Molecular formulas with digits → name (PubChem name search handles formulas)
		{"H2O", "name"},
		{"C6H12O6", "name"},
		{"CH4", "name"},
		{"CO2", "name"},
		{"NaCl", "name"}, // no digits but contains 'a' — not organic SMILES chars

		// Edge cases
		{"", "name"},
		{"   ", "name"},
	}

	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			got := smiles.DetectInputType(c.input)
			if got != c.want {
				t.Errorf("DetectInputType(%q) = %q, want %q", c.input, got, c.want)
			}
		})
	}
}

// ── Lookup — live network tests ───────────────────────────────────────────────
// These make real PubChem requests.
// Run with: go test ./smiles/... -v -run TestLookup -timeout 60s

func TestLookup_ByName_Water(t *testing.T) {
	p, err := smiles.Lookup("water")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "H2O" {
		t.Errorf("formula: got %q, want H2O", p.MolecularFormula)
	}
	if p.SMILES == "" {
		t.Error("expected non-empty SMILES for water")
	}
	if p.InChIKey == "" {
		t.Error("expected non-empty InChIKey for water")
	}
	if p.IUPACName == "" {
		t.Error("expected non-empty IUPAC name for water")
	}
	if p.CID == 0 {
		t.Error("expected non-zero CID")
	}
	if p.InputType != "name" {
		t.Errorf("input_type: got %q, want name", p.InputType)
	}
}

func TestLookup_BySMILES_Ethanol(t *testing.T) {
	p, err := smiles.Lookup("CCO")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "C2H6O" {
		t.Errorf("formula: got %q, want C2H6O", p.MolecularFormula)
	}
	if p.CID != 702 {
		t.Errorf("CID: got %d, want 702", p.CID)
	}
	if p.InputType != "smiles" {
		t.Errorf("input_type: got %q, want smiles", p.InputType)
	}
}

func TestLookup_BySMILES_Caffeine(t *testing.T) {
	p, err := smiles.Lookup("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "C8H10N4O2" {
		t.Errorf("formula: got %q, want C8H10N4O2", p.MolecularFormula)
	}
	if p.SMILES == "" {
		t.Error("expected non-empty SMILES")
	}
}

func TestLookup_ByCID_NaCl(t *testing.T) {
	p, err := smiles.Lookup("5234")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "ClNa" && p.MolecularFormula != "NaCl" {
		t.Errorf("formula: got %q, want NaCl or ClNa", p.MolecularFormula)
	}
	if p.CID != 5234 {
		t.Errorf("CID: got %d, want 5234", p.CID)
	}
	if p.InputType != "cid" {
		t.Errorf("input_type: got %q, want cid", p.InputType)
	}
}

func TestLookup_ByInChIKey_Ethanol(t *testing.T) {
	p, err := smiles.Lookup("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "C2H6O" {
		t.Errorf("formula: got %q, want C2H6O", p.MolecularFormula)
	}
	if p.InputType != "inchikey" {
		t.Errorf("input_type: got %q, want inchikey", p.InputType)
	}
}

func TestLookup_AllFieldsPopulated(t *testing.T) {
	p, err := smiles.Lookup("ethanol")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	checks := []struct {
		name  string
		value string
	}{
		{"MolecularFormula", p.MolecularFormula},
		{"IUPACName", p.IUPACName},
		{"SMILES", p.SMILES},
		{"InChI", p.InChI},
		{"InChIKey", p.InChIKey},
	}
	for _, c := range checks {
		if c.value == "" {
			t.Errorf("field %s is empty", c.name)
		}
	}
	if p.MolecularWeight.IsZero() {
		t.Error("MolecularWeight is zero")
	}
	if p.CID == 0 {
		t.Error("CID is zero")
	}
}

func TestLookup_NotFound(t *testing.T) {
	_, err := smiles.Lookup("xyzzy_not_a_compound_123456789")
	if err == nil {
		t.Error("expected error for unknown compound")
	}
}
