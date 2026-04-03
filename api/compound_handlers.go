package api

import (
	"net/http"

	"chemhelper/smiles"
)

// POST /api/compound/lookup
// Accepts any common compound identifier (name, SMILES, InChI, InChIKey, CID,
// molecular formula) and returns all available representations from PubChem.
// Body: { "input": "water" }
func HandleCompoundLookup(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Input string `json:"input"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.Input == "" {
		writeError(w, http.StatusBadRequest, "input is required")
		return
	}
	props, err := smiles.Lookup(req.Input)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, CompoundLookupResponse{
		CID:              props.CID,
		MolecularFormula: props.MolecularFormula,
		MolecularWeight:  props.MolecularWeight.String(),
		IUPACName:        props.IUPACName,
		SMILES:           props.SMILES,
		InChI:            props.InChI,
		InChIKey:         props.InChIKey,
		InputType:        props.InputType,
	})
}

// POST /api/compound/resolve
// Resolves a SMILES string to compound properties via PubChem.
// Body: { "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C" }
// Returns: molecular formula, molar weight, IUPAC name, and PubChem CID.
func HandleCompoundResolve(w http.ResponseWriter, r *http.Request) {
	var req struct {
		SMILES string `json:"smiles"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.SMILES == "" {
		writeError(w, http.StatusBadRequest, "smiles is required")
		return
	}
	props, err := smiles.Resolve(req.SMILES)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, CompoundResponse{
		CID:              props.CID,
		MolecularFormula: props.MolecularFormula,
		MolecularWeight:  props.MolecularWeight.String(),
		IUPACName:        props.IUPACName,
	})
}
