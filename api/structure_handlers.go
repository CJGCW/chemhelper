package api

import (
	"net/http"
	"strings"

	"chemhelper/smiles"
	"chemhelper/structure"
)

// isSMILES reports whether s contains characters that only appear in SMILES
// notation and never in plain molecular formulas.
func isSMILES(s string) bool {
	return strings.ContainsAny(s, "=#[]@/\\.")
}

// POST /api/structure/lewis
// Body: { "input": "H2O" }
//
// Accepts molecular formulas ("H2O", "CaCl2", "Fe(OH)2"), parenthesised groups
// ("Mn(NO3)2"), or SMILES strings ("O=C=O", "[Ca+2].[OH-].[OH-]").
// SMILES inputs are resolved to a molecular formula via PubChem before
// Lewis structure generation. Returns 422 with a descriptive error on failure.
func HandleLewis(w http.ResponseWriter, r *http.Request) {
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

	formula := req.Input
	if isSMILES(req.Input) {
		props, err := smiles.Resolve(req.Input)
		if err != nil {
			writeError(w, http.StatusUnprocessableEntity, "SMILES resolution failed: "+err.Error())
			return
		}
		formula = props.MolecularFormula
	}

	ls, err := structure.LookupLewisWithError(formula)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}

	writeJSON(w, http.StatusOK, ls)
}
