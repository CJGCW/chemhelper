package estimate

import (
	"math"
	"testing"
)

// hasPeakNear returns true if peaks contains a peak within tolerance of x.
func hasPeakNear(peaks []Peak, x float64, tolerance float64) bool {
	for _, p := range peaks {
		if math.Abs(p.X-x) <= tolerance {
			return true
		}
	}
	return false
}

// peakCount returns the number of peaks.
func peakCount(peaks []Peak) int { return len(peaks) }

// totalIntegration sums integration values across all ¹H peaks.
func totalIntegration(peaks []Peak) int {
	total := 0
	for _, p := range peaks {
		total += p.Integration
	}
	return total
}

// ── Core correctness tests ─────────────────────────────────────────────────────

func TestEstimateAcetone(t *testing.T) {
	result, err := Estimate("CC(=O)C")
	if err != nil {
		t.Fatal(err)
	}

	// IR: C=O ketone ~1715 ± 30
	if !hasPeakNear(result.IR, 1715, 30) {
		t.Error("acetone IR: missing C=O stretch near 1715")
	}
	// IR: sp3 C-H ~2920
	if !hasPeakNear(result.IR, 2920, 60) {
		t.Error("acetone IR: missing sp3 C-H stretch near 2920")
	}
	// IR: no O-H
	if hasPeakNear(result.IR, 3400, 200) {
		t.Error("acetone IR: spurious O-H peak")
	}

	// ¹H NMR: two equivalent methyls → 1 peak, integration=6
	if len(result.NMR1H) != 1 {
		t.Errorf("acetone ¹H NMR: expected 1 peak, got %d", len(result.NMR1H))
	} else {
		if result.NMR1H[0].Integration != 6 {
			t.Errorf("acetone ¹H NMR: expected integration 6, got %d", result.NMR1H[0].Integration)
		}
		if result.NMR1H[0].Splitting != "s" {
			t.Errorf("acetone ¹H NMR: expected singlet, got %q", result.NMR1H[0].Splitting)
		}
	}

	// ¹³C NMR: carbonyl ~205, methyl ~30
	if !hasPeakNear(result.NMR13C, 205, 15) {
		t.Error("acetone ¹³C NMR: missing carbonyl near 205")
	}
	if !hasPeakNear(result.NMR13C, 30, 15) {
		t.Error("acetone ¹³C NMR: missing methyl near 30")
	}
	// Two distinct carbons in acetone (C=O and CH3, where two CH3's are equivalent)
	if len(result.NMR13C) != 2 {
		t.Errorf("acetone ¹³C NMR: expected 2 peaks, got %d", len(result.NMR13C))
	}

	// MS: M+ at 58
	if !hasPeakNear(result.MS, 58, 1) {
		t.Errorf("acetone MS: missing M+ at 58 (MW=%f)", result.MolecularWeight)
	}

	// Molecular properties
	if result.MolecularFormula != "C3H6O" {
		t.Errorf("acetone formula: expected C3H6O, got %q", result.MolecularFormula)
	}
	if math.Abs(result.MolecularWeight-58.08) > 0.5 {
		t.Errorf("acetone MW: expected ~58.08, got %f", result.MolecularWeight)
	}
	if result.DegreesUnsaturation != 1 {
		t.Errorf("acetone DoU: expected 1, got %d", result.DegreesUnsaturation)
	}
}

func TestEstimateEthanol(t *testing.T) {
	result, err := Estimate("CCO")
	if err != nil {
		t.Fatal(err)
	}

	// IR: O-H stretch
	if !hasPeakNear(result.IR, 3400, 200) {
		t.Error("ethanol IR: missing O-H stretch")
	}
	// IR: sp3 C-H
	if !hasPeakNear(result.IR, 2920, 80) {
		t.Error("ethanol IR: missing sp3 C-H stretch")
	}
	// IR: C-O
	if !hasPeakNear(result.IR, 1050, 100) {
		t.Error("ethanol IR: missing C-O stretch")
	}

	// ¹H NMR: CH3 (~0.9–1.3 ppm, t), CH2 (~3.5 ppm), OH (~2.5 ppm) → 3 peaks
	if len(result.NMR1H) != 3 {
		t.Errorf("ethanol ¹H NMR: expected 3 peaks, got %d", len(result.NMR1H))
	}
	total := totalIntegration(result.NMR1H)
	if total != 6 {
		t.Errorf("ethanol ¹H NMR: expected total integration 6, got %d", total)
	}

	// MS: M+ at 46
	if !hasPeakNear(result.MS, 46, 1) {
		t.Error("ethanol MS: missing M+ at 46")
	}
}

func TestEstimateAceticAcid(t *testing.T) {
	result, err := Estimate("CC(=O)O")
	if err != nil {
		t.Fatal(err)
	}

	// IR: broad O-H ~2700 (acid), C=O ~1710
	if !hasPeakNear(result.IR, 2700, 500) {
		t.Error("acetic acid IR: missing broad O-H stretch")
	}
	if !hasPeakNear(result.IR, 1710, 30) {
		t.Error("acetic acid IR: missing C=O stretch")
	}

	// ¹H NMR: COOH at ~11.5, CH3 at ~2.3
	if !hasPeakNear(result.NMR1H, 11.5, 1) {
		t.Error("acetic acid ¹H NMR: missing COOH peak")
	}
	if !hasPeakNear(result.NMR1H, 2.3, 0.8) {
		t.Error("acetic acid ¹H NMR: missing CH3 peak")
	}

	// MS: M+ at 60
	if !hasPeakNear(result.MS, 60, 1) {
		t.Error("acetic acid MS: missing M+ at 60")
	}
}

func TestEstimateBenzene(t *testing.T) {
	result, err := Estimate("c1ccccc1")
	if err != nil {
		t.Fatal(err)
	}

	// IR: aromatic C=C, C-H sp2
	if !hasPeakNear(result.IR, 3050, 50) {
		t.Error("benzene IR: missing sp2 C-H stretch")
	}
	if !hasPeakNear(result.IR, 1600, 50) {
		t.Error("benzene IR: missing aromatic C=C stretch")
	}

	// ¹H NMR: 1 peak at ~7.3 ppm, integration=6 (all ArH equivalent)
	if len(result.NMR1H) != 1 {
		t.Errorf("benzene ¹H NMR: expected 1 peak, got %d", len(result.NMR1H))
	} else if result.NMR1H[0].Integration != 6 {
		t.Errorf("benzene ¹H NMR: expected integration 6, got %d", result.NMR1H[0].Integration)
	}

	// MS: M+ at 78
	if !hasPeakNear(result.MS, 78, 1) {
		t.Error("benzene MS: missing M+ at 78")
	}

	// Formula
	if result.MolecularFormula != "C6H6" {
		t.Errorf("benzene formula: expected C6H6, got %q", result.MolecularFormula)
	}
	if result.DegreesUnsaturation != 4 {
		t.Errorf("benzene DoU: expected 4, got %d", result.DegreesUnsaturation)
	}
}

func TestEstimateEthylAcetate(t *testing.T) {
	result, err := Estimate("CC(=O)OCC")
	if err != nil {
		t.Fatal(err)
	}

	// IR: ester C=O ~1740, C-O ~1200, sp3 C-H
	if !hasPeakNear(result.IR, 1740, 30) {
		t.Error("ethyl acetate IR: missing ester C=O near 1740")
	}
	if !hasPeakNear(result.IR, 1200, 100) {
		t.Error("ethyl acetate IR: missing ester C-O stretch")
	}
	if !hasPeakNear(result.IR, 2920, 80) {
		t.Error("ethyl acetate IR: missing sp3 C-H")
	}

	// MS: M+ at 88
	if !hasPeakNear(result.MS, 88, 1) {
		t.Error("ethyl acetate MS: missing M+ at 88")
	}
	// m/z 43: acetyl cation
	if !hasPeakNear(result.MS, 43, 1) {
		t.Error("ethyl acetate MS: missing m/z 43 acetyl cation")
	}
}

func TestEstimateButanal(t *testing.T) {
	result, err := Estimate("CCCC=O")
	if err != nil {
		t.Fatal(err)
	}

	// IR: aldehyde C=O ~1725, aldehyde C-H ~2820, sp3 C-H
	if !hasPeakNear(result.IR, 1725, 30) {
		t.Error("butanal IR: missing aldehyde C=O")
	}
	if !hasPeakNear(result.IR, 2820, 50) {
		t.Error("butanal IR: missing aldehyde C-H")
	}

	// ¹H NMR: aldehyde H at ~9.7 ppm
	if !hasPeakNear(result.NMR1H, 9.7, 0.5) {
		t.Error("butanal ¹H NMR: missing aldehyde H")
	}

	// MS: M+ at 72, M-29 (loss of CHO)
	if !hasPeakNear(result.MS, 72, 1) {
		t.Error("butanal MS: missing M+ at 72")
	}
	if !hasPeakNear(result.MS, 43, 1) {
		t.Error("butanal MS: missing M-29 at 43")
	}
}

func TestEstimateAcetonitrile(t *testing.T) {
	result, err := Estimate("CC#N")
	if err != nil {
		t.Fatal(err)
	}

	// IR: C≡N ~2240, sp3 C-H
	if !hasPeakNear(result.IR, 2240, 30) {
		t.Error("acetonitrile IR: missing C≡N stretch")
	}

	// MS: M+ at 41
	if !hasPeakNear(result.MS, 41, 1) {
		t.Error("acetonitrile MS: missing M+ at 41")
	}
}

func TestEstimateDiethylEther(t *testing.T) {
	result, err := Estimate("CCOCC")
	if err != nil {
		t.Fatal(err)
	}

	// IR: C-O ether ~1100, sp3 C-H; no C=O, no O-H
	if !hasPeakNear(result.IR, 1100, 100) {
		t.Error("diethyl ether IR: missing C-O stretch")
	}
	if hasPeakNear(result.IR, 1715, 50) {
		t.Error("diethyl ether IR: spurious C=O peak")
	}
	if hasPeakNear(result.IR, 3400, 200) {
		t.Error("diethyl ether IR: spurious O-H peak")
	}

	// MS: M+ at 74
	if !hasPeakNear(result.MS, 74, 1) {
		t.Error("diethyl ether MS: missing M+ at 74")
	}
}

func TestEstimateAniline(t *testing.T) {
	result, err := Estimate("Nc1ccccc1")
	if err != nil {
		t.Fatal(err)
	}

	// IR: N-H primary amine (two bands ~3380/3280), aromatic
	if !hasPeakNear(result.IR, 3380, 100) && !hasPeakNear(result.IR, 3280, 100) {
		t.Error("aniline IR: missing N-H stretch")
	}
	if !hasPeakNear(result.IR, 1600, 50) {
		t.Error("aniline IR: missing aromatic C=C")
	}

	// MS: M+ at 93
	if !hasPeakNear(result.MS, 93, 1) {
		t.Error("aniline MS: missing M+ at 93")
	}
	// m/z 77: phenyl
	if !hasPeakNear(result.MS, 77, 1) {
		t.Error("aniline MS: missing m/z 77 phenyl cation")
	}
}

func TestEstimatePropylamine(t *testing.T) {
	result, err := Estimate("CCCN")
	if err != nil {
		t.Fatal(err)
	}

	// IR: N-H primary amine
	if !hasPeakNear(result.IR, 3380, 100) && !hasPeakNear(result.IR, 3280, 100) {
		t.Error("propylamine IR: missing N-H stretch")
	}

	// MS: M+ at 59
	if !hasPeakNear(result.MS, 59, 1) {
		t.Error("propylamine MS: missing M+ at 59")
	}
}

func TestEstimateBromobenzene(t *testing.T) {
	result, err := Estimate("Brc1ccccc1")
	if err != nil {
		t.Fatal(err)
	}

	// MS: M+ at 156 (79Br), M+2 at 158 (81Br) — about equal intensity
	if !hasPeakNear(result.MS, 156, 1) {
		t.Error("bromobenzene MS: missing M+ at 156")
	}
	if !hasPeakNear(result.MS, 158, 1) {
		t.Error("bromobenzene MS: missing M+2 (81Br) at 158")
	}
}

func TestEstimateChloroform(t *testing.T) {
	result, err := Estimate("ClCCl")
	if err != nil {
		t.Fatal(err)
	}

	// MS: M+ at 84 (CH2Cl2), M+2 due to 37Cl
	// Actually ClCCl is CH2Cl2: MW = 12+2+35+35 = 84
	if !hasPeakNear(result.MS, 84, 1) {
		t.Error("CH2Cl2 MS: missing M+ at 84")
	}
	if !hasPeakNear(result.MS, 86, 1) {
		t.Error("CH2Cl2 MS: missing M+2 (37Cl)")
	}
}

func TestEstimatePropyne(t *testing.T) {
	result, err := Estimate("CC#C")
	if err != nil {
		t.Fatal(err)
	}

	// IR: ≡C-H ~3300, C≡C ~2150, sp3 C-H
	if !hasPeakNear(result.IR, 3300, 50) {
		t.Error("propyne IR: missing ≡C-H stretch")
	}
	if !hasPeakNear(result.IR, 2150, 50) {
		t.Error("propyne IR: missing C≡C stretch")
	}

	// ¹H NMR: terminal alkyne H at ~2.5, CH3 at ~0.9–2.3
	if !hasPeakNear(result.NMR1H, 2.5, 0.5) {
		t.Error("propyne ¹H NMR: missing terminal alkyne H near 2.5")
	}
}

func TestEstimateToluene(t *testing.T) {
	result, err := Estimate("Cc1ccccc1")
	if err != nil {
		t.Fatal(err)
	}

	// IR: sp3 C-H, sp2/ArH, aromatic C=C
	if !hasPeakNear(result.IR, 2920, 80) {
		t.Error("toluene IR: missing sp3 C-H (CH3)")
	}
	if !hasPeakNear(result.IR, 1600, 50) {
		t.Error("toluene IR: missing aromatic C=C")
	}

	// ¹H NMR: ArH (7.3 ppm, 5H) + CH3 (2.3 ppm, 3H)
	// The aromatic H's in monosubstituted benzene are 3 groups (o/m/p)
	// OR all merge — either way total ArH integration should be 5
	arHTotal := 0
	for _, p := range result.NMR1H {
		if p.X > 6.5 {
			arHTotal += p.Integration
		}
	}
	if arHTotal != 5 {
		t.Errorf("toluene ¹H NMR: expected 5 ArH, got %d", arHTotal)
	}

	// MS: M+ at 92, m/z 91 (tropylium)
	if !hasPeakNear(result.MS, 92, 1) {
		t.Error("toluene MS: missing M+ at 92")
	}
	if !hasPeakNear(result.MS, 91, 1) {
		t.Error("toluene MS: missing tropylium at 91")
	}
}

func TestEstimateAcetamide(t *testing.T) {
	result, err := Estimate("CC(=O)N")
	if err != nil {
		t.Fatal(err)
	}

	// IR: amide C=O ~1660, N-H primary amide
	if !hasPeakNear(result.IR, 1660, 40) {
		t.Error("acetamide IR: missing amide C=O")
	}
	if !hasPeakNear(result.IR, 3300, 150) {
		t.Error("acetamide IR: missing N-H stretch")
	}
}

func TestEstimateNitrobenzene(t *testing.T) {
	result, err := Estimate("[N+](=O)[O-]c1ccccc1")
	if err != nil {
		t.Fatal(err)
	}

	// IR: N=O stretches ~1530 and ~1350
	if !hasPeakNear(result.IR, 1530, 40) {
		t.Error("nitrobenzene IR: missing N=O asym stretch ~1530")
	}
	if !hasPeakNear(result.IR, 1350, 40) {
		t.Error("nitrobenzene IR: missing N=O sym stretch ~1350")
	}
}

func TestEstimateAceticAnhydride(t *testing.T) {
	result, err := Estimate("CC(=O)OC(=O)C")
	if err != nil {
		t.Fatal(err)
	}

	// IR: two C=O bands for anhydride ~1820 and ~1760
	if !hasPeakNear(result.IR, 1820, 40) {
		t.Error("acetic anhydride IR: missing high-frequency C=O")
	}
	if !hasPeakNear(result.IR, 1760, 40) {
		t.Error("acetic anhydride IR: missing low-frequency C=O")
	}
}

func TestEstimateAcetylChloride(t *testing.T) {
	result, err := Estimate("CC(=O)Cl")
	if err != nil {
		t.Fatal(err)
	}

	// IR: acyl chloride C=O ~1800
	if !hasPeakNear(result.IR, 1800, 30) {
		t.Error("acetyl chloride IR: missing C=O near 1800")
	}
}

func TestEstimate1Propene(t *testing.T) {
	result, err := Estimate("CC=C")
	if err != nil {
		t.Fatal(err)
	}

	// IR: C=C alkene ~1650, sp2 C-H ~3050
	if !hasPeakNear(result.IR, 1650, 40) {
		t.Error("propene IR: missing C=C stretch")
	}
	if !hasPeakNear(result.IR, 3050, 60) {
		t.Error("propene IR: missing sp2 C-H stretch")
	}

	// ¹H NMR: should have alkene H peaks (~5.5 ppm)
	if !hasPeakNear(result.NMR1H, 5.5, 0.5) {
		t.Error("propene ¹H NMR: missing alkene H")
	}
}

func TestEstimatePhenol(t *testing.T) {
	result, err := Estimate("Oc1ccccc1")
	if err != nil {
		t.Fatal(err)
	}

	// IR: O-H stretch, aromatic
	if !hasPeakNear(result.IR, 3400, 300) {
		t.Error("phenol IR: missing O-H stretch")
	}
	if !hasPeakNear(result.IR, 1600, 50) {
		t.Error("phenol IR: missing aromatic C=C")
	}

	// MS: M+ at 94
	if !hasPeakNear(result.MS, 94, 1) {
		t.Error("phenol MS: missing M+ at 94")
	}
}

// ── Parse correctness ─────────────────────────────────────────────────────────

func TestParseMolecularFormulas(t *testing.T) {
	cases := []struct {
		smi     string
		formula string
		mw      float64
	}{
		{"C", "CH4", 16.04},
		{"CC", "C2H6", 30.07},
		{"CCO", "C2H6O", 46.07},
		{"CC(=O)C", "C3H6O", 58.08},
		{"c1ccccc1", "C6H6", 78.11},
		{"CC#N", "C2H3N", 41.05},
	}

	for _, tc := range cases {
		result, err := Estimate(tc.smi)
		if err != nil {
			t.Errorf("%s: parse error: %v", tc.smi, err)
			continue
		}
		if result.MolecularFormula != tc.formula {
			t.Errorf("%s: formula expected %s, got %s", tc.smi, tc.formula, result.MolecularFormula)
		}
		if math.Abs(result.MolecularWeight-tc.mw) > 0.5 {
			t.Errorf("%s: MW expected %.2f, got %.2f", tc.smi, tc.mw, result.MolecularWeight)
		}
	}
}

func TestEstimateEmptySMILES(t *testing.T) {
	_, err := Estimate("")
	if err == nil {
		t.Error("expected error for empty SMILES")
	}
}

func TestEstimateInvalidSMILES(t *testing.T) {
	_, err := Estimate("invalid!@#$")
	if err == nil {
		t.Error("expected error for invalid SMILES")
	}
}

func TestAllPeaksValid(t *testing.T) {
	smiles := []string{
		"CC(=O)C", "CCO", "CC(=O)O", "c1ccccc1", "Cc1ccccc1",
		"CC(=O)OCC", "CCCC=O", "CC#N", "CCOCC", "CCCN",
		"CC=C", "CC#CC", "Oc1ccccc1", "CC(=O)N", "Nc1ccccc1",
	}
	for _, smi := range smiles {
		result, err := Estimate(smi)
		if err != nil {
			t.Errorf("%s: unexpected error: %v", smi, err)
			continue
		}
		for _, p := range result.IR {
			if p.X <= 0 || p.Y <= 0 || p.Label == "" {
				t.Errorf("%s IR: invalid peak %+v", smi, p)
			}
		}
		for _, p := range result.NMR1H {
			if p.X < 0 || p.X > 15 || p.Y <= 0 || p.Integration <= 0 {
				t.Errorf("%s ¹H NMR: invalid peak %+v", smi, p)
			}
		}
		for _, p := range result.NMR13C {
			if p.X < 0 || p.X > 230 || p.Y <= 0 {
				t.Errorf("%s ¹³C NMR: invalid peak %+v", smi, p)
			}
		}
		for _, p := range result.MS {
			if p.X <= 0 || p.Y <= 0 {
				t.Errorf("%s MS: invalid peak %+v", smi, p)
			}
		}
	}
}
