// Package spectral provides types and stub implementations for IR, NMR, and MS analysis.
package spectral

// SpectrumType identifies which type of spectrum is represented.
type SpectrumType string

const (
	IR      SpectrumType = "ir"
	HNMR    SpectrumType = "1h_nmr"
	CNMR    SpectrumType = "13c_nmr"
	MassSpec SpectrumType = "mass_spec"
)

// Peak represents a single absorption or fragment peak.
type Peak struct {
	X     float64 `json:"x"`     // wavenumber cm⁻¹ (IR), ppm (NMR), m/z (MS)
	Y     float64 `json:"y"`     // intensity / %T / abundance
	Label string  `json:"label"` // functional group assignment
	Width float64 `json:"width"` // Gaussian half-width parameter
}

// NMRPeak extends Peak with splitting and integration data.
type NMRPeak struct {
	Peak
	Splitting   string  `json:"splitting"`   // singlet, doublet, triplet, quartet, multiplet
	Integration float64 `json:"integration"` // relative area (number of protons for ¹H)
}

// SpectrumPrediction is the result of a predicted IR or MS spectrum.
type SpectrumPrediction struct {
	Type  SpectrumType `json:"type"`
	Peaks []Peak       `json:"peaks"`
	Notes []string     `json:"notes"`
}

// NMRPrediction is the result of a predicted ¹H or ¹³C NMR spectrum.
type NMRPrediction struct {
	Type  SpectrumType `json:"type"`
	Peaks []NMRPeak    `json:"peaks"`
	Notes []string     `json:"notes"`
}

// Candidate is a compound proposed by the analyzer as consistent with the observed spectrum.
type Candidate struct {
	Name          string   `json:"name"`
	Formula       string   `json:"formula"`
	SMILES        string   `json:"smiles"`
	Confidence    float64  `json:"confidence"`     // 0–100
	SVG           string   `json:"svg"`             // structure SVG (may be empty)
	MatchedGroups []string `json:"matched_groups"`  // functional groups confirmed by spectrum
}

// AnalysisResult summarises the output of AnalyzeSpectrum.
type AnalysisResult struct {
	SpectrumType     SpectrumType `json:"spectrum_type"`
	Peaks            []Peak       `json:"peaks"`
	FunctionalGroups []string     `json:"functional_groups"`
	Candidates       []Candidate  `json:"candidates"`
}
