package api

import (
	"fmt"
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

// maxChargeMagnitude is the largest absolute ionic charge this endpoint will accept.
// Real compounds rarely exceed ±8; values beyond that almost certainly indicate a
// mis-formatted formula (e.g. "NO33-" parsed as N at charge −32).
const maxChargeMagnitude = 8

// GET /api/structure/random
// Returns a randomly generated symmetric AB_n Lewis structure.
func HandleRandomStructure(w http.ResponseWriter, r *http.Request) {
	ls, err := structure.GenerateRandom()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, ls)
}

// POST /api/structure/lewis
// Body: { "input": "H2O" }
//    or { "input": "NO3", "charge": -1 }
//
// `input` accepts molecular formulas ("H2O", "CaCl2"), parenthesised groups
// ("Fe(OH)2", "Mn(NO3)2"), SMILES strings ("O=C=O"), or legacy ionic notation
// ("OH-", "Al3+") when `charge` is omitted.
//
// `charge` (optional integer) overrides charge extraction from the input string.
// When provided the formula must be a bare formula with no +/- suffix.
// Parenthesised formulas do not support an explicit charge override.
//
// Returns 422 with a descriptive error on failure.
func HandleLewis(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Input  string `json:"input"`
		Charge *int   `json:"charge"` // nil = extract from suffix; non-nil = explicit override
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.Input == "" {
		writeError(w, http.StatusBadRequest, "input is required")
		return
	}

	// Validate explicit charge before doing any resolution work.
	if req.Charge != nil {
		c := *req.Charge
		if c < -maxChargeMagnitude || c > maxChargeMagnitude {
			writeError(w, http.StatusBadRequest,
				fmt.Sprintf("charge %d is outside the supported range [%d, +%d]; check your formula for embedded charge notation",
					c, -maxChargeMagnitude, maxChargeMagnitude))
			return
		}
		if strings.ContainsAny(req.Input, "+-") {
			writeError(w, http.StatusBadRequest,
				"when 'charge' is provided the formula must not contain '+' or '-'; remove the charge suffix from the formula")
			return
		}
		if strings.ContainsRune(req.Input, '(') {
			writeError(w, http.StatusBadRequest,
				"explicit 'charge' is not supported for parenthesised formulas; omit the charge field and use ionic suffix notation instead")
			return
		}
	}

	formula := req.Input
	if isSMILES(req.Input) {
		if req.Charge != nil {
			writeError(w, http.StatusBadRequest,
				"explicit 'charge' cannot be combined with SMILES input; omit the charge field")
			return
		}
		props, err := smiles.Resolve(req.Input)
		if err != nil {
			writeError(w, http.StatusUnprocessableEntity, "SMILES resolution failed: "+err.Error())
			return
		}
		formula = props.MolecularFormula
	}

	var ls *structure.LewisStructure
	var err error

	if req.Charge != nil {
		ls, err = structure.LookupLewisWithCharge(formula, *req.Charge)
	} else {
		ls, err = structure.LookupLewisWithError(formula)
	}
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}

	writeJSON(w, http.StatusOK, ls)
}
