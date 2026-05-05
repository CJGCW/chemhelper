package estimate

import "fmt"

// Peak represents a single spectral peak for SpectrumViewer rendering.
type Peak struct {
	X           float64 `json:"x"`
	Y           float64 `json:"y"`
	Label       string  `json:"label"`
	Width       float64 `json:"width"`
	Splitting   string  `json:"splitting,omitempty"`
	Integration int     `json:"integration,omitempty"`
	Source      string  `json:"source,omitempty"`
}

// EstimateResult contains estimated spectra for all four spectrum types.
type EstimateResult struct {
	SMILES              string   `json:"smiles"`
	MolecularFormula    string   `json:"molecular_formula"`
	MolecularWeight     float64  `json:"molecular_weight"`
	DegreesUnsaturation int      `json:"degrees_unsaturation"`
	IR                  []Peak   `json:"ir"`
	NMR1H               []Peak   `json:"nmr_1h"`
	NMR13C              []Peak   `json:"nmr_13c"`
	MS                  []Peak   `json:"ms"`
	Warnings            []string `json:"warnings,omitempty"`
}

// Estimate returns estimated IR, ¹H NMR, ¹³C NMR, and MS spectra for a SMILES string.
// Returns an error if the SMILES string cannot be parsed.
func Estimate(smi string) (EstimateResult, error) {
	mol, err := parseSMILES(smi)
	if err != nil {
		return EstimateResult{}, fmt.Errorf("unparseable SMILES: %w", err)
	}

	return EstimateResult{
		SMILES:              smi,
		MolecularFormula:    mol.molecularFormula(),
		MolecularWeight:     mol.molecularWeight(),
		DegreesUnsaturation: mol.degreesUnsaturation(),
		IR:                  estimateIR(mol),
		NMR1H:               estimateHNMR(mol),
		NMR13C:              estimateCNMR(mol),
		MS:                  estimateMS(mol),
	}, nil
}
