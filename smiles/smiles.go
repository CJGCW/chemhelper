// Package smiles resolves SMILES strings to compound properties via the
// PubChem REST API. Results are cached in memory for the lifetime of the
// process to avoid redundant network calls.
package smiles

import (
	"encoding/json"
	"fmt"
	"net/http"
	"net/url"
	"sync"
	"time"

	"github.com/shopspring/decimal"
)

const (
	pubchemBase = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles"
	timeout     = 10 * time.Second
)

var client = &http.Client{Timeout: timeout}

// Properties holds the resolved properties of a compound.
type Properties struct {
	MolecularFormula string
	MolecularWeight  decimal.Decimal // g/mol
	IUPACName        string
	CID              int // PubChem Compound ID
}

// cache stores resolved SMILES to avoid redundant PubChem calls.
var (
	cacheMu sync.RWMutex
	cache   = make(map[string]Properties)
)

// Resolve resolves a SMILES string to compound properties via PubChem.
// Results are cached — repeated calls with the same SMILES string are free.
func Resolve(smiles string) (Properties, error) {
	cacheMu.RLock()
	if p, ok := cache[smiles]; ok {
		cacheMu.RUnlock()
		return p, nil
	}
	cacheMu.RUnlock()

	p, err := fetchFromPubChem(smiles)
	if err != nil {
		return Properties{}, err
	}

	cacheMu.Lock()
	cache[smiles] = p
	cacheMu.Unlock()

	return p, nil
}

// ResolveToMolarMass is a convenience wrapper that returns only the molar mass.
func ResolveToMolarMass(smiles string) (decimal.Decimal, error) {
	p, err := Resolve(smiles)
	if err != nil {
		return decimal.Zero, err
	}
	return p.MolecularWeight, nil
}

// pubchemResponse mirrors the PubChem property table JSON envelope.
type pubchemResponse struct {
	PropertyTable struct {
		Properties []struct {
			CID              int     `json:"CID"`
			MolecularFormula string  `json:"MolecularFormula"`
			MolecularWeight  float64 `json:"MolecularWeight,string"`
			IUPACName        string  `json:"IUPACName"`
		} `json:"Properties"`
	} `json:"PropertyTable"`
}

type pubchemFaultResponse struct {
	Fault struct {
		Code    string `json:"Code"`
		Message string `json:"Message"`
		Details []string `json:"Details"`
	} `json:"Fault"`
}

func fetchFromPubChem(smiles string) (Properties, error) {
	encoded := url.PathEscape(smiles)
	reqURL := fmt.Sprintf(
		"%s/%s/property/MolecularFormula,MolecularWeight,IUPACName/JSON",
		pubchemBase, encoded,
	)

	resp, err := client.Get(reqURL)
	if err != nil {
		return Properties{}, fmt.Errorf("PubChem request failed: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode == http.StatusNotFound {
		// PubChem returns a structured fault on 404
		var fault pubchemFaultResponse
		if err := json.NewDecoder(resp.Body).Decode(&fault); err == nil && fault.Fault.Message != "" {
			return Properties{}, fmt.Errorf("SMILES not found in PubChem: %s", fault.Fault.Message)
		}
		return Properties{}, fmt.Errorf("SMILES %q not found in PubChem", smiles)
	}
	if resp.StatusCode != http.StatusOK {
		return Properties{}, fmt.Errorf("PubChem returned status %d", resp.StatusCode)
	}

	var result pubchemResponse
	if err := json.NewDecoder(resp.Body).Decode(&result); err != nil {
		return Properties{}, fmt.Errorf("failed to decode PubChem response: %w", err)
	}
	if len(result.PropertyTable.Properties) == 0 {
		return Properties{}, fmt.Errorf("PubChem returned no properties for SMILES %q", smiles)
	}

	p := result.PropertyTable.Properties[0]
	return Properties{
		CID:              p.CID,
		MolecularFormula: p.MolecularFormula,
		MolecularWeight:  decimal.NewFromFloat(p.MolecularWeight),
		IUPACName:        p.IUPACName,
	}, nil
}
