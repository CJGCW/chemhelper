package api

import (
	"fmt"
	"net/http"
	"sort"
	"strings"

	"chemhelper/thermo"

	"github.com/shopspring/decimal"
)

// knownSolvents maps lowercase names to their Solvent constants.
var knownSolvents = map[string]thermo.Solvent{
	"water":               thermo.Water,
	"benzene":             thermo.Benzene,
	"ethanol":             thermo.Ethanol,
	"cyclohexane":         thermo.Cyclohexane,
	"carbontetrachloride": thermo.CarbonTetrachloride,
	"chloroform":          thermo.Chloroform,
}

func parseSolvent(name string) (thermo.Solvent, error) {
	key := strings.ToLower(strings.ReplaceAll(name, " ", ""))
	s, ok := knownSolvents[key]
	if !ok {
		return thermo.Solvent{}, &unknownSolventError{name}
	}
	return s, nil
}

type unknownSolventError struct{ name string }

func (e *unknownSolventError) Error() string {
	return "unknown solvent \"" + e.name + "\"; see GET /api/thermo/solvents for valid options"
}

// colligativeRequest is the shared body shape for BPE and FPD endpoints.
type colligativeRequest struct {
	Solvent        string `json:"solvent"`
	Molality       string `json:"molality"`
	VantHoffFactor string `json:"vant_hoff_factor,omitempty"`
}

func (req colligativeRequest) parse() (thermo.Solvent, decimal.Decimal, decimal.Decimal, error) {
	solvent, err := parseSolvent(req.Solvent)
	if err != nil {
		return thermo.Solvent{}, decimal.Zero, decimal.Zero, err
	}
	molality, err := decimal.NewFromString(req.Molality)
	if err != nil {
		return thermo.Solvent{}, decimal.Zero, decimal.Zero,
			fmt.Errorf("invalid molality %q: %w", req.Molality, err)
	}
	i, err := parseVantHoff(req.VantHoffFactor)
	if err != nil {
		return thermo.Solvent{}, decimal.Zero, decimal.Zero, err
	}
	return solvent, molality, i, nil
}

// GET /api/thermo/solvents
// Returns the list of known solvents sorted by name.
func HandleSolvents(w http.ResponseWriter, r *http.Request) {
	resp := make([]SolventResponse, 0, len(knownSolvents))
	for _, s := range knownSolvents {
		resp = append(resp, SolventResponse{
			Name:          s.Name,
			BoilingPoint:  s.BoilingPoint.String(),
			FreezingPoint: s.FreezingPoint.String(),
			Kb:            s.Kb.String(),
			Kf:            s.Kf.String(),
		})
	}
	sort.Slice(resp, func(i, j int) bool { return resp[i].Name < resp[j].Name })
	writeJSON(w, http.StatusOK, resp)
}

// POST /api/thermo/bpe
// Body: { "solvent": "water", "molality": "1.0", "vant_hoff_factor": "2" }
// Returns both ΔTb and the resulting new boiling point.
func HandleBPE(w http.ResponseWriter, r *http.Request) {
	var req colligativeRequest
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	solvent, molality, i, err := req.parse()
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	calc := thermo.BoilingPointElevation{
		Solvent:        solvent,
		Molality:       molality,
		VantHoffFactor: i,
	}
	delta, err := calc.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	newBP, err := calc.NewBoilingPoint()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, ColligativeResponse{
		Delta:    delta.Value.String(),
		NewPoint: newBP.Value.String(),
		Unit:     delta.Unit,
		Steps:    delta.Steps,
	})
}

// POST /api/thermo/fpd
// Body: { "solvent": "water", "molality": "1.0", "vant_hoff_factor": "2" }
// Returns both ΔTf and the resulting new freezing point.
func HandleFPD(w http.ResponseWriter, r *http.Request) {
	var req colligativeRequest
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	solvent, molality, i, err := req.parse()
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	calc := thermo.FreezingPointDepression{
		Solvent:        solvent,
		Molality:       molality,
		VantHoffFactor: i,
	}
	delta, err := calc.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	newFP, err := calc.NewFreezingPoint()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, ColligativeResponse{
		Delta:    delta.Value.String(),
		NewPoint: newFP.Value.String(),
		Unit:     delta.Unit,
		Steps:    delta.Steps,
	})
}

// POST /api/thermo/molality-from-bpe
// Body: { "solvent": "water", "delta_tb": "1.024", "vant_hoff_factor": "2" }
func HandleMolalityFromBPE(w http.ResponseWriter, r *http.Request) {
	var req BackCalcTbRequest
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	solvent, err := parseSolvent(req.Solvent)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	deltaTb, err := decimal.NewFromString(req.DeltaTb)
	if err != nil {
		writeError(w, http.StatusBadRequest, fmt.Sprintf("invalid delta_tb %q: %v", req.DeltaTb, err))
		return
	}
	i, err := parseVantHoff(req.VantHoffFactor)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	result, err := thermo.FindMolalityFromBPE{
		Solvent:        solvent,
		DeltaTb:        deltaTb,
		VantHoffFactor: i,
	}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}

// POST /api/thermo/molality-from-fpd
// Body: { "solvent": "water", "delta_tf": "3.72", "vant_hoff_factor": "2" }
func HandleMolalityFromFPD(w http.ResponseWriter, r *http.Request) {
	var req BackCalcTfRequest
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	solvent, err := parseSolvent(req.Solvent)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	deltaTf, err := decimal.NewFromString(req.DeltaTf)
	if err != nil {
		writeError(w, http.StatusBadRequest, fmt.Sprintf("invalid delta_tf %q: %v", req.DeltaTf, err))
		return
	}
	i, err := parseVantHoff(req.VantHoffFactor)
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	result, err := thermo.FindMolalityFromFPD{
		Solvent:        solvent,
		DeltaTf:        deltaTf,
		VantHoffFactor: i,
	}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}
