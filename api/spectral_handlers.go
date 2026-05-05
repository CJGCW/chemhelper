package api

import (
	"net/http"

	"chemhelper/spectral"
	"chemhelper/spectral/estimate"
)

// POST /api/spectral/predict-ir
// Body: { "smiles": "CC(=O)O" }
func HandlePredictIR(w http.ResponseWriter, r *http.Request) {
	var req struct {
		SMILES string `json:"smiles"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.SMILES == "" {
		writeError(w, http.StatusUnprocessableEntity, "smiles must not be empty")
		return
	}
	pred, err := spectral.PredictIR(req.SMILES)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, pred)
}

// POST /api/spectral/predict-hnmr
// Body: { "smiles": "CC(=O)O" }
func HandlePredictHNMR(w http.ResponseWriter, r *http.Request) {
	var req struct {
		SMILES string `json:"smiles"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.SMILES == "" {
		writeError(w, http.StatusUnprocessableEntity, "smiles must not be empty")
		return
	}
	pred, err := spectral.PredictHNMR(req.SMILES)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, pred)
}

// POST /api/spectral/predict-cnmr
// Body: { "smiles": "CC(=O)O" }
func HandlePredictCNMR(w http.ResponseWriter, r *http.Request) {
	var req struct {
		SMILES string `json:"smiles"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.SMILES == "" {
		writeError(w, http.StatusUnprocessableEntity, "smiles must not be empty")
		return
	}
	pred, err := spectral.PredictCNMR(req.SMILES)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, pred)
}

// POST /api/spectral/predict-ms
// Body: { "smiles": "CC(=O)O" }
func HandlePredictMS(w http.ResponseWriter, r *http.Request) {
	var req struct {
		SMILES string `json:"smiles"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.SMILES == "" {
		writeError(w, http.StatusUnprocessableEntity, "smiles must not be empty")
		return
	}
	pred, err := spectral.PredictMS(req.SMILES)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, pred)
}

// POST /api/spectral/analyze
// Body: { "type": "ir", "peaks": [ { "x": 1715, "y": 0.1, "label": "", "width": 20 } ] }
func HandleSpectralAnalyze(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Type  spectral.SpectrumType `json:"type"`
		Peaks []struct {
			X     float64 `json:"x"`
			Y     float64 `json:"y"`
			Label string  `json:"label"`
			Width float64 `json:"width"`
		} `json:"peaks"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.Type == "" {
		writeError(w, http.StatusUnprocessableEntity, "type must be one of: ir, 1h_nmr, 13c_nmr, mass_spec")
		return
	}

	peaks := make([]spectral.Peak, len(req.Peaks))
	for i, p := range req.Peaks {
		peaks[i] = spectral.Peak{X: p.X, Y: p.Y, Label: p.Label, Width: p.Width}
	}

	result, err := spectral.AnalyzeSpectrum(req.Type, peaks)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, result)
}

// POST /api/spectral/estimate
// Body: { "smiles": "CC(=O)C" }
func HandleSpectralEstimate(w http.ResponseWriter, r *http.Request) {
	var req struct {
		SMILES string `json:"smiles"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.SMILES == "" {
		writeError(w, http.StatusUnprocessableEntity, "smiles must not be empty")
		return
	}
	result, err := estimate.Estimate(req.SMILES)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, result)
}
