package render

import (
	"strings"
	"testing"
)

// Minimal SDF snippets that exercise different structural features.
// Format: 3 header lines + counts + atom block + bond block + M  END

const ethanolSDF = `Ethanol
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
`

// acetone: CC(=O)C  — has a double bond (C=O)
const acetoneSDF = `Acetone
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
M  END
$$$$
`

// benzene: planar ring (aromatic bond type 4)
const benzeneSDF = `Benzene
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  6  1  4  0
M  END
$$$$
`

// (R)-2-bromobutane: has a wedge bond (stereo 1)
const brombutaneSDF = `2-Bromobutane
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  1
M  END
$$$$
`

// minimal 1-atom SDF (single oxygen atom — edge case)
const singleAtomSDF = `Water
     RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
`

func TestParseSDF_Ethanol(t *testing.T) {
	mol, err := parseSDF(ethanolSDF)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(mol.Atoms) != 3 {
		t.Errorf("want 3 atoms, got %d", len(mol.Atoms))
	}
	if len(mol.Bonds) != 2 {
		t.Errorf("want 2 bonds, got %d", len(mol.Bonds))
	}
	if mol.Atoms[2].Symbol != "O" {
		t.Errorf("want 3rd atom O, got %q", mol.Atoms[2].Symbol)
	}
}

func TestParseSDF_Acetone(t *testing.T) {
	mol, err := parseSDF(acetoneSDF)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// Check the C=O bond has order 2
	found := false
	for _, b := range mol.Bonds {
		if (b.A1 == 1 && b.A2 == 2) || (b.A1 == 2 && b.A2 == 1) {
			if b.Order != 2 {
				t.Errorf("C=O bond: want order 2, got %d", b.Order)
			}
			found = true
		}
	}
	if !found {
		t.Error("did not find C=O bond")
	}
}

func TestParseSDF_Benzene(t *testing.T) {
	mol, err := parseSDF(benzeneSDF)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(mol.Atoms) != 6 {
		t.Errorf("want 6 atoms, got %d", len(mol.Atoms))
	}
	for _, b := range mol.Bonds {
		if b.Order != 4 {
			t.Errorf("benzene bond order: want 4, got %d", b.Order)
		}
	}
}

func TestParseSDF_WedgeBond(t *testing.T) {
	mol, err := parseSDF(brombutaneSDF)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	found := false
	for _, b := range mol.Bonds {
		if b.Stereo == 1 {
			found = true
		}
	}
	if !found {
		t.Error("did not find wedge bond (stereo=1)")
	}
}

func TestToSVG_ContainsElements(t *testing.T) {
	mol, _ := parseSDF(ethanolSDF)
	svg := toSVG(mol, 200, 150)

	// Must have SVG root element
	if !strings.Contains(svg, "<svg") {
		t.Error("SVG missing <svg element")
	}
	// Oxygen label must appear (skeletal style shows heteroatoms)
	if !strings.Contains(svg, ">O<") {
		t.Error("SVG missing O atom label")
	}
	// Carbon must NOT appear as explicit text
	if strings.Contains(svg, ">C<") {
		t.Error("SVG should not render explicit C label")
	}
}

func TestToSVG_DoubleBondHasParallelLines(t *testing.T) {
	mol, _ := parseSDF(acetoneSDF)
	svg := toSVG(mol, 200, 150)
	// A double bond produces 2 <line> elements for the C=O bond;
	// the total line count should be > the total bond count
	lineCount := strings.Count(svg, "<line")
	bondCount := len(mol.Bonds)
	if lineCount <= bondCount {
		t.Errorf("double bond should produce extra lines: got %d lines for %d bonds", lineCount, bondCount)
	}
}

func TestToSVG_WedgeBond(t *testing.T) {
	mol, _ := parseSDF(brombutaneSDF)
	svg := toSVG(mol, 200, 150)
	if !strings.Contains(svg, "<polygon") {
		t.Error("wedge bond should produce a <polygon>")
	}
}

func TestToSVG_SingleAtom(t *testing.T) {
	mol, _ := parseSDF(singleAtomSDF)
	svg := toSVG(mol, 200, 150)
	if !strings.Contains(svg, "<svg") {
		t.Error("single-atom SVG must still produce <svg>")
	}
}

func TestToSVG_CurrentColor(t *testing.T) {
	mol, _ := parseSDF(benzeneSDF)
	svg := toSVG(mol, 200, 150)
	if !strings.Contains(svg, "currentColor") {
		t.Error("SVG must use currentColor for dark-mode compatibility")
	}
}

// ── Local renderer tests ──────────────────────────────────────────────────────

func TestNeedsLocalRenderer(t *testing.T) {
	cases := []struct {
		smiles string
		want   bool
	}{
		{"[R]", true},
		{"[Nu]", true},
		{"[LG]", true},
		{"[X]", true},
		{"[E+]", true},
		{"[Nu-]", true},
		{"[R]C([R])([R])[LG]", true},
		{"O", false},
		{"CC(=O)O", false},
		{"[H3O+]", false},
		{"[OH-]", false},
		{"C1CO1", false},
	}
	for _, tc := range cases {
		got := needsLocalRenderer(tc.smiles)
		if got != tc.want {
			t.Errorf("needsLocalRenderer(%q) = %v, want %v", tc.smiles, got, tc.want)
		}
	}
}

func TestParseSMILESLocal_SingleAtom(t *testing.T) {
	mol := parseSMILESLocal("[Nu]")
	if mol == nil {
		t.Fatal("nil mol")
	}
	if len(mol.atoms) != 1 {
		t.Fatalf("want 1 atom, got %d", len(mol.atoms))
	}
	if mol.atoms[0].symbol != "Nu" {
		t.Errorf("want symbol Nu, got %q", mol.atoms[0].symbol)
	}
}

func TestParseSMILESLocal_TwoAtoms(t *testing.T) {
	mol := parseSMILESLocal("[R][LG]")
	if mol == nil {
		t.Fatal("nil mol")
	}
	if len(mol.atoms) != 2 {
		t.Fatalf("want 2 atoms, got %d", len(mol.atoms))
	}
	if len(mol.bonds) != 1 {
		t.Fatalf("want 1 bond, got %d", len(mol.bonds))
	}
}

func TestParseSMILESLocal_Tetrahedral(t *testing.T) {
	mol := parseSMILESLocal("[R]C([R])([R])[LG]")
	if mol == nil {
		t.Fatal("nil mol")
	}
	if len(mol.atoms) != 5 {
		t.Fatalf("want 5 atoms, got %d: %+v", len(mol.atoms), mol.atoms)
	}
	// C should have 4 bonds
	cIdx := -1
	for i, a := range mol.atoms {
		if a.symbol == "C" {
			cIdx = i
			break
		}
	}
	if cIdx < 0 {
		t.Fatal("no C atom found")
	}
	if len(mol.adj[cIdx]) != 4 {
		t.Errorf("C atom degree want 4, got %d", len(mol.adj[cIdx]))
	}
}

func TestParseBracket_H3OPlus(t *testing.T) {
	sym, chg, hcnt := parseBracket("H3O+")
	if sym != "O" {
		t.Errorf("want sym O, got %q", sym)
	}
	if chg != 1 {
		t.Errorf("want charge +1, got %d", chg)
	}
	if hcnt != 3 {
		t.Errorf("want hCount 3, got %d", hcnt)
	}
}

func TestParseBracket_OHMinus(t *testing.T) {
	sym, chg, hcnt := parseBracket("OH-")
	if sym != "O" {
		t.Errorf("want sym O, got %q", sym)
	}
	if chg != -1 {
		t.Errorf("want charge -1, got %d", chg)
	}
	if hcnt != 1 {
		t.Errorf("want hCount 1, got %d", hcnt)
	}
}

func TestParseBracket_CustomLabel(t *testing.T) {
	sym, chg, hcnt := parseBracket("R")
	if sym != "R" {
		t.Errorf("want sym R, got %q", sym)
	}
	if chg != 0 || hcnt != 0 {
		t.Errorf("want charge 0, hCount 0; got %d, %d", chg, hcnt)
	}

	sym, chg, _ = parseBracket("Nu-")
	if sym != "Nu" || chg != -1 {
		t.Errorf("Nu-: want sym Nu chg -1; got %q %d", sym, chg)
	}
}

func TestRenderLocalSMILES_SingleLabel(t *testing.T) {
	svg, err := renderLocalSMILES("[Nu]", 100, 80, false)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(svg, "<svg") {
		t.Error("missing <svg")
	}
	if !strings.Contains(svg, "Nu") {
		t.Error("missing Nu label in SVG")
	}
}

func TestRenderLocalSMILES_Tetrahedral(t *testing.T) {
	svg, err := renderLocalSMILES("[R]C([R])([R])[LG]", 200, 150, false)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(svg, "R") {
		t.Error("missing R label")
	}
	if !strings.Contains(svg, "LG") {
		t.Error("missing LG label")
	}
	// C should be visible (local renderer shows all atoms)
	if !strings.Contains(svg, "C") {
		t.Error("missing C label")
	}
}

func TestRenderLocalSMILES_LonePairs(t *testing.T) {
	svg, err := renderLocalSMILES("[R]O", 120, 90, true)
	if err != nil {
		t.Fatal(err)
	}
	// Should contain circle elements for lone pair dots
	if !strings.Contains(svg, "<circle") {
		t.Error("lone pairs should produce <circle> elements")
	}
}

func TestLonePairCount(t *testing.T) {
	cases := []struct {
		sym    string
		charge int
		bonds  int
		want   int
	}{
		{"O", 0, 2, 2},  // water: 2 LP
		{"O", 0, 1, 2},  // alcohol O: still 2 LP in simplified model (6 - 0 - 1 = 5, /2 = 2)
		{"O", -1, 1, 3}, // hydroxide: 6 - (-1) - 1 = 6... wait that's wrong
		{"N", 0, 3, 1},  // amine: 1 LP
		{"N", 1, 4, 0},  // ammonium: 0 LP
		{"F", 0, 1, 3},  // fluoride: 3 LP
		{"C", 0, 4, 0},  // sp3 C: 0 LP
	}
	for _, tc := range cases {
		got := lonePairCount(tc.sym, tc.charge, tc.bonds)
		if got != tc.want {
			t.Errorf("lonePairCount(%q, %d, %d) = %d, want %d",
				tc.sym, tc.charge, tc.bonds, got, tc.want)
		}
	}
}

func TestRender_CustomLabel(t *testing.T) {
	svg, err := Render("[R]C([R])([R])[LG]", 200, 150, RenderOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(svg, "<svg") {
		t.Error("missing svg element")
	}
	if !strings.Contains(svg, "LG") {
		t.Error("LG label not in SVG")
	}
}
