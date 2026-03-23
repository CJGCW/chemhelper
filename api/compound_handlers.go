package api

import (
	"net/http"

	"chemhelper/smiles"
)

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
