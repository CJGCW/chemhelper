package api

import (
	"fmt"
	"strings"

	"chemhelper/smiles"
	"chemhelper/units"

	"github.com/shopspring/decimal"
)

// --- Request types ---

// MassRequest is the JSON representation of a mass measurement.
// Unit: "gram" (default), "pound", "ounce"
// Prefix: "none" (default), "kilo", "hecto", "deca", "deci", "centi", "milli", "micro"
type MassRequest struct {
	Value  string `json:"value"`
	Unit   string `json:"unit,omitempty"`
	Prefix string `json:"prefix,omitempty"`
}

// VolumeRequest is the JSON representation of a volume measurement.
// The unit is always litres; Prefix selects the scale.
// Prefix: "none" (default), "kilo", "hecto", "deca", "deci", "centi", "milli", "micro"
type VolumeRequest struct {
	Value  string `json:"value"`
	Prefix string `json:"prefix,omitempty"`
}

// CompoundRequest identifies a compound either by SMILES string or by a
// pre-known molar mass. Exactly one of SMILES or MolarMass must be set.
//
// SMILES example:  {"smiles": "[Na+].[Cl-]"}
// MolarMass example: {"molar_mass": "58.44"}
type CompoundRequest struct {
	SMILES    string `json:"smiles,omitempty"`
	MolarMass string `json:"molar_mass,omitempty"`
}

// ResolveMolarMass returns the molar mass in g/mol. If SMILES is set it
// queries PubChem (with caching); otherwise it parses MolarMass directly.
func (c CompoundRequest) ResolveMolarMass() (decimal.Decimal, error) {
	if c.SMILES != "" && c.MolarMass != "" {
		return decimal.Zero, fmt.Errorf("provide either smiles or molar_mass, not both")
	}
	if c.SMILES != "" {
		return smiles.ResolveToMolarMass(c.SMILES)
	}
	if c.MolarMass != "" {
		v, err := decimal.NewFromString(c.MolarMass)
		if err != nil {
			return decimal.Zero, fmt.Errorf("invalid molar_mass %q: %w", c.MolarMass, err)
		}
		return v, nil
	}
	return decimal.Zero, fmt.Errorf("either smiles or molar_mass must be provided")
}

// ResolveProperties returns the full PubChem properties for a SMILES compound.
// Returns an error if MolarMass was provided instead of SMILES.
func (c CompoundRequest) ResolveProperties() (smiles.Properties, error) {
	if c.SMILES == "" {
		return smiles.Properties{}, fmt.Errorf("smiles is required for full property resolution")
	}
	return smiles.Resolve(c.SMILES)
}
type BackCalcTbRequest struct {
	Solvent        string `json:"solvent"`
	DeltaTb        string `json:"delta_tb"`
	VantHoffFactor string `json:"vant_hoff_factor,omitempty"`
}

// BackCalcTfRequest is the body for POST /api/thermo/molality-from-fpd.
type BackCalcTfRequest struct {
	Solvent        string `json:"solvent"`
	DeltaTf        string `json:"delta_tf"`
	VantHoffFactor string `json:"vant_hoff_factor,omitempty"`
}

// --- Response types ---

// CalcResponse is the standard JSON envelope for all calculation results.
type CalcResponse struct {
	Value   string   `json:"value"`
	Unit    string   `json:"unit"`
	SigFigs int      `json:"sig_figs,omitempty"`
	Steps   []string `json:"steps,omitempty"`
}

// ColligativeResponse is returned by BPE and FPD endpoints.
// It includes both the delta (elevation or depression) and the resulting
// new boiling or freezing point, so the client never needs two round-trips.
type ColligativeResponse struct {
	Delta    string   `json:"delta"`
	NewPoint string   `json:"new_point"`
	Unit     string   `json:"unit"`
	Steps    []string `json:"steps,omitempty"`
}

// ErrorResponse is returned on 4xx/5xx.
type ErrorResponse struct {
	Error string `json:"error"`
}

// SolventResponse is the JSON representation of a Solvent.
type SolventResponse struct {
	Name          string `json:"name"`
	BoilingPoint  string `json:"boiling_point_c"`
	FreezingPoint string `json:"freezing_point_c"`
	Kb            string `json:"kb"`
	Kf            string `json:"kf"`
}

// CompoundResponse is returned by POST /api/compound/resolve.
type CompoundResponse struct {
	CID              int    `json:"cid"`
	MolecularFormula string `json:"molecular_formula"`
	MolecularWeight  string `json:"molecular_weight"`
	IUPACName        string `json:"iupac_name"`
}

// CompoundLookupResponse is returned by POST /api/compound/lookup.
type CompoundLookupResponse struct {
	CID              int    `json:"cid"`
	MolecularFormula string `json:"molecular_formula"`
	MolecularWeight  string `json:"molecular_weight"`
	IUPACName        string `json:"iupac_name"`
	SMILES           string `json:"smiles"`
	InChI            string `json:"inchi"`
	InChIKey         string `json:"inchi_key"`
	InputType        string `json:"input_type"`
}

// ElementResponse is the JSON representation of a periodic table element.
type ElementResponse struct {
	AtomicNumber      int     `json:"atomic_number"`
	Symbol            string  `json:"symbol"`
	Name              string  `json:"name"`
	AtomicWeight      string  `json:"atomic_weight"`
	Electronegativity float64 `json:"electronegativity"`
	VanDerWaalsRadius float64 `json:"van_der_waals_radius_pm"`
	Group             int     `json:"group"`
	Period            int     `json:"period"`
	GroupName         string  `json:"group_name,omitempty"`
}

// CompoundElementResponse is one entry in a parsed compound result.
type CompoundElementResponse struct {
	Symbol string `json:"symbol"`
	Name   string `json:"name"`
	Moles  string `json:"moles"`
}

// --- Conversion helpers ---

func parsePrefix(s string) (units.Prefix, error) {
	switch strings.ToLower(s) {
	case "", "none":
		return units.None, nil
	case "kilo":
		return units.Kilo, nil
	case "hecto":
		return units.Hecto, nil
	case "deca":
		return units.Deca, nil
	case "deci":
		return units.Deci, nil
	case "centi":
		return units.Centi, nil
	case "milli":
		return units.Milli, nil
	case "micro":
		return units.Micro, nil
	default:
		return units.None, fmt.Errorf("unknown prefix %q", s)
	}
}

func parseMassUnit(s string) (units.MassUnit, error) {
	switch strings.ToLower(s) {
	case "", "gram":
		return units.Gram, nil
	case "pound":
		return units.Pound, nil
	case "ounce":
		return units.Ounce, nil
	default:
		return units.Gram, fmt.Errorf("unknown mass unit %q", s)
	}
}

// ToMass converts a MassRequest to a units.Mass.
func (r MassRequest) ToMass() (units.Mass, error) {
	val, err := decimal.NewFromString(r.Value)
	if err != nil {
		return units.Mass{}, fmt.Errorf("invalid mass value %q: %w", r.Value, err)
	}
	unit, err := parseMassUnit(r.Unit)
	if err != nil {
		return units.Mass{}, err
	}
	prefix, err := parsePrefix(r.Prefix)
	if err != nil {
		return units.Mass{}, err
	}
	return units.NewMass(val, unit, prefix)
}

// ToVolume converts a VolumeRequest to a units.Volume.
func (r VolumeRequest) ToVolume() (units.Volume, error) {
	val, err := decimal.NewFromString(r.Value)
	if err != nil {
		return units.Volume{}, fmt.Errorf("invalid volume value %q: %w", r.Value, err)
	}
	prefix, err := parsePrefix(r.Prefix)
	if err != nil {
		return units.Volume{}, err
	}
	return units.NewVolume(val, prefix)
}

// parseVantHoff parses an optional van't Hoff factor string.
// An empty string returns zero (which the domain layer treats as 1).
func parseVantHoff(s string) (decimal.Decimal, error) {
	if s == "" {
		return decimal.Zero, nil
	}
	v, err := decimal.NewFromString(s)
	if err != nil {
		return decimal.Zero, fmt.Errorf("invalid vant_hoff_factor %q: %w", s, err)
	}
	return v, nil
}
