// Package structure provides Lewis structure data for common intro-chemistry molecules.
package structure

// LewisAtom describes a single atom in a Lewis structure.
type LewisAtom struct {
	ID           string `json:"id"`
	Element      string `json:"element"`
	LonePairs    int    `json:"lone_pairs"` // number of lone PAIRS (each = 2 electrons)
	FormalCharge int    `json:"formal_charge"`
}

// LewisBond describes a bond between two atoms.
type LewisBond struct {
	From  string `json:"from"`
	To    string `json:"to"`
	Order int    `json:"order"` // 1 = single, 2 = double, 3 = triple
}

// LewisStructure is the full Lewis structure representation returned by the API.
type LewisStructure struct {
	Name                  string      `json:"name"`
	Formula               string      `json:"formula"` // display formula with unicode subscripts
	Charge                int         `json:"charge"`  // overall ionic charge (0 for neutral)
	TotalValenceElectrons int         `json:"total_valence_electrons"`
	Geometry              string      `json:"geometry"` // VSEPR geometry hint for frontend layout
	Atoms                 []LewisAtom `json:"atoms"`
	Bonds                 []LewisBond `json:"bonds"`
	Steps                 []string    `json:"steps"`
	Notes                 string      `json:"notes,omitempty"`
}

// registry maps a normalised formula key → hand-verified static structure.
var registry map[string]*LewisStructure

// LookupLewis returns the Lewis structure for the given input string.
// It first checks the static registry (hand-verified entries), then falls back
// to the algorithmic generator.
func LookupLewis(input string) (*LewisStructure, bool) {
	ls, err := LookupLewisWithError(input)
	if err != nil {
		return nil, false
	}
	return ls, true
}

// LookupLewisWithError is like LookupLewis but returns the underlying error
// so the HTTP handler can report a meaningful message to the client.
func LookupLewisWithError(input string) (*LewisStructure, error) {
	return generate(input)
}

// LookupLewisWithCharge generates a Lewis structure for the given molecular
// formula and an explicit integer charge, bypassing the charge-suffix parser.
// The formula must be a bare formula with no trailing +/- notation.
// Parenthesised formulas (e.g. "Ca(OH)2") are not supported here; use
// LookupLewisWithError for those.
func LookupLewisWithCharge(formula string, charge int) (*LewisStructure, error) {
	return generateFromParsed(formula, charge)
}
