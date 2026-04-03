package smiles

import (
	"encoding/json"
	"fmt"
	"net/http"
	"net/url"
	"regexp"
	"strconv"
	"strings"

	"github.com/shopspring/decimal"
)

const propList = "MolecularFormula,MolecularWeight,IUPACName,SMILES,InChI,InChIKey"

// ExtendedProperties holds all retrievable compound data from PubChem.
type ExtendedProperties struct {
	CID              int
	MolecularFormula string
	MolecularWeight  decimal.Decimal
	IUPACName        string
	SMILES           string
	InChI            string
	InChIKey         string
	InputType        string // "cid" | "inchi" | "inchikey" | "smiles" | "name"
}

var inchiKeyRe = regexp.MustCompile(`^[A-Z]{14}-[A-Z]{10}-[A-Z]$`)

// smilesOnlyRe matches strings composed exclusively of organic-subset SMILES
// atoms (no special chars, no digits). Such strings can't be molecular formulas
// (formulas use digits for counts) or common names (which contain non-SMILES
// letters like 'a', 'e', 't'). Examples: "CCO", "CC", "c1ccccc1" minus digits.
// Atoms covered: B C F H I N O P S and two-char Cl Br, plus aromatic bcnops.
var smilesOnlyRe = regexp.MustCompile(`^([BCFHINOPSbcnops]|[Cc]l|[Bb]r)+$`)

// aromaticDigitSmilesRe matches SMILES with ring-closure digits alongside
// lowercase aromatic atoms — e.g. "c1ccccc1" for benzene.
var aromaticDigitSmilesRe = regexp.MustCompile(`^([BCFHINOPSbcnops]|[Cc]l|[Bb]r|\d)+$`)

// DetectInputType identifies what kind of chemical identifier the input is.
// Exported so the API handler can include it in the response.
func DetectInputType(s string) string {
	s = strings.TrimSpace(s)
	if s == "" {
		return "name"
	}
	if _, err := strconv.Atoi(s); err == nil {
		return "cid"
	}
	if strings.HasPrefix(strings.ToUpper(s), "INCHI=") {
		return "inchi"
	}
	if inchiKeyRe.MatchString(s) {
		return "inchikey"
	}
	// SMILES-specific characters not found in IUPAC names or molecular formulas.
	if strings.ContainsAny(s, "=#@\\/[]") {
		return "smiles"
	}
	// Parentheses alone appear in IUPAC names, so only treat as SMILES when
	// combined with lowercase element symbols (common in SMILES notation).
	if strings.Contains(s, "(") && strings.ContainsAny(s, "cnosp") {
		return "smiles"
	}
	// Plain organic SMILES (e.g. "CCO", "CC", "Cl") have no digits or special
	// chars, only organic-subset atom letters. Molecular formulas always use
	// digits for element counts (H2O, CH4), so digit-free all-atom strings are
	// almost certainly SMILES, not formulas or names.
	if smilesOnlyRe.MatchString(s) {
		return "smiles"
	}
	// Aromatic SMILES with ring-closure digits (e.g. "c1ccccc1") contain
	// lowercase aromatic atom symbols alongside digits.
	if aromaticDigitSmilesRe.MatchString(s) && strings.ContainsAny(s, "bcnops") {
		return "smiles"
	}
	return "name"
}

// Lookup auto-detects the input type and fetches extended compound properties
// from PubChem. Returns an error if PubChem cannot resolve the input.
func Lookup(input string) (ExtendedProperties, error) {
	input = strings.TrimSpace(input)
	inputType := DetectInputType(input)

	resp, err := fetchExtended(input, inputType)
	if err != nil {
		return ExtendedProperties{}, err
	}
	resp.InputType = inputType
	return resp, nil
}

func fetchExtended(input, inputType string) (ExtendedProperties, error) {
	const base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

	var (
		httpResp *http.Response
		err      error
	)

	switch inputType {
	case "inchi":
		// InChI strings contain slashes, so they must be sent via POST form.
		form := url.Values{"inchi": {input}}
		httpResp, err = client.PostForm(
			fmt.Sprintf("%s/inchi/property/%s/JSON", base, propList),
			form,
		)
	default:
		var segment string
		switch inputType {
		case "cid":
			segment = "cid/" + url.PathEscape(input)
		case "inchikey":
			segment = "inchikey/" + url.PathEscape(input)
		case "smiles":
			segment = "smiles/" + url.PathEscape(input)
		default: // "name"
			segment = "name/" + url.PathEscape(input)
		}
		httpResp, err = client.Get(fmt.Sprintf("%s/%s/property/%s/JSON", base, segment, propList))
	}

	if err != nil {
		return ExtendedProperties{}, fmt.Errorf("PubChem request failed: %w", err)
	}
	defer httpResp.Body.Close()

	if httpResp.StatusCode == http.StatusNotFound {
		var fault pubchemFaultResponse
		if decErr := json.NewDecoder(httpResp.Body).Decode(&fault); decErr == nil && fault.Fault.Message != "" {
			return ExtendedProperties{}, fmt.Errorf("not found in PubChem: %s", fault.Fault.Message)
		}
		return ExtendedProperties{}, fmt.Errorf("%q not found in PubChem", input)
	}
	if httpResp.StatusCode != http.StatusOK {
		return ExtendedProperties{}, fmt.Errorf("PubChem returned status %d", httpResp.StatusCode)
	}

	var result extendedPubchemResponse
	if err := json.NewDecoder(httpResp.Body).Decode(&result); err != nil {
		return ExtendedProperties{}, fmt.Errorf("failed to decode PubChem response: %w", err)
	}
	if len(result.PropertyTable.Properties) == 0 {
		return ExtendedProperties{}, fmt.Errorf("PubChem returned no properties for %q", input)
	}

	p := result.PropertyTable.Properties[0]
	return ExtendedProperties{
		CID:              p.CID,
		MolecularFormula: p.MolecularFormula,
		MolecularWeight:  decimal.NewFromFloat(p.MolecularWeight),
		IUPACName:        p.IUPACName,
		SMILES:           p.SMILES,
		InChI:            p.InChI,
		InChIKey:         p.InChIKey,
	}, nil
}

type extendedPubchemResponse struct {
	PropertyTable struct {
		Properties []struct {
			CID              int     `json:"CID"`
			MolecularFormula string  `json:"MolecularFormula"`
			MolecularWeight  float64 `json:"MolecularWeight,string"`
			IUPACName        string  `json:"IUPACName"`
			SMILES           string  `json:"SMILES"`
			InChI            string  `json:"InChI"`
			InChIKey         string  `json:"InChIKey"`
		} `json:"Properties"`
	} `json:"PropertyTable"`
}
