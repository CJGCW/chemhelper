package smiles_test

import (
	"testing"

	"chemhelper/smiles"

	"github.com/shopspring/decimal"
)

// These tests make real network calls to PubChem.
// Run with: go test ./smiles/... -v -timeout 30s

func TestResolve_Water(t *testing.T) {
	p, err := smiles.Resolve("O")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "H2O" {
		t.Errorf("expected formula H2O, got %s", p.MolecularFormula)
	}
	// Water molar mass ~18.015 g/mol — check within 0.01
	expected := decimal.NewFromFloat(18.015)
	diff := p.MolecularWeight.Sub(expected).Abs()
	if diff.GreaterThan(decimal.NewFromFloat(0.01)) {
		t.Errorf("expected molar mass ~18.015, got %s", p.MolecularWeight)
	}
	if p.CID == 0 {
		t.Errorf("expected non-zero CID")
	}
}

func TestResolve_NaCl(t *testing.T) {
	p, err := smiles.Resolve("[Na+].[Cl-]")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expected := decimal.NewFromFloat(58.44)
	diff := p.MolecularWeight.Sub(expected).Abs()
	if diff.GreaterThan(decimal.NewFromFloat(0.01)) {
		t.Errorf("expected molar mass ~58.44, got %s", p.MolecularWeight)
	}
}

func TestResolve_Glucose(t *testing.T) {
	p, err := smiles.Resolve("C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "C6H12O6" {
		t.Errorf("expected formula C6H12O6, got %s", p.MolecularFormula)
	}
	expected := decimal.NewFromFloat(180.156)
	diff := p.MolecularWeight.Sub(expected).Abs()
	if diff.GreaterThan(decimal.NewFromFloat(0.01)) {
		t.Errorf("expected molar mass ~180.156, got %s", p.MolecularWeight)
	}
}

func TestResolve_Caffeine(t *testing.T) {
	// Caffeine: C8H10N4O2, molar mass 194.19 g/mol
	p, err := smiles.Resolve("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if p.MolecularFormula != "C8H10N4O2" {
		t.Errorf("expected formula C8H10N4O2, got %s", p.MolecularFormula)
	}
}

func TestResolve_CachesResult(t *testing.T) {
	// Two calls with the same SMILES — second should come from cache.
	// We can't directly observe the cache, but we can verify both return
	// the same value and neither errors.
	p1, err := smiles.Resolve("CCO") // ethanol
	if err != nil {
		t.Fatalf("first call error: %v", err)
	}
	p2, err := smiles.Resolve("CCO")
	if err != nil {
		t.Fatalf("second call error: %v", err)
	}
	if !p1.MolecularWeight.Equal(p2.MolecularWeight) {
		t.Errorf("cached result differs: %s vs %s", p1.MolecularWeight, p2.MolecularWeight)
	}
}

func TestResolve_InvalidSMILES(t *testing.T) {
	_, err := smiles.Resolve("not-a-smiles-string!!!")
	if err == nil {
		t.Error("expected error for invalid SMILES but got none")
	}
}

func TestResolveToMolarMass(t *testing.T) {
	mw, err := smiles.ResolveToMolarMass("O") // water
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	expected := decimal.NewFromFloat(18.015)
	diff := mw.Sub(expected).Abs()
	if diff.GreaterThan(decimal.NewFromFloat(0.01)) {
		t.Errorf("expected ~18.015, got %s", mw)
	}
}
