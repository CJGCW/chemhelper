package api

import (
	"net/http"

	"chemhelper/solution"

	"github.com/shopspring/decimal"
)

// POST /api/solution/molarity
// Body: { "moles": "1.5", "volume": { "value": "0.500", "prefix": "none" } }
func HandleMolarity(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Moles  string        `json:"moles"`
		Volume VolumeRequest `json:"volume"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	moles, err := decimal.NewFromString(req.Moles)
	if err != nil {
		writeError(w, http.StatusBadRequest, "invalid moles value: "+err.Error())
		return
	}
	vol, err := req.Volume.ToVolume()
	if err != nil {
		writeError(w, http.StatusBadRequest, err.Error())
		return
	}
	result, err := solution.FindMolarity{Moles: moles, Volume: vol}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}

// POST /api/solution/molarity-from-mass
// Body (SMILES):     { "mass": { "value": "58.44", "unit": "gram" }, "compound": {"smiles": "[Na+].[Cl-]"}, "volume": { "value": "1" } }
// Body (molar mass): { "mass": { "value": "58.44", "unit": "gram" }, "compound": {"molar_mass": "58.44"}, "volume": { "value": "1" } }
func HandleMolarityFromMass(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Mass     MassRequest     `json:"mass"`
		Compound CompoundRequest `json:"compound"`
		Volume   VolumeRequest   `json:"volume"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	mass, err := req.Mass.ToMass()
	if err != nil {
		writeError(w, http.StatusBadRequest, "mass: "+err.Error())
		return
	}
	molarMass, err := req.Compound.ResolveMolarMass()
	if err != nil {
		writeError(w, http.StatusBadRequest, "compound: "+err.Error())
		return
	}
	vol, err := req.Volume.ToVolume()
	if err != nil {
		writeError(w, http.StatusBadRequest, "volume: "+err.Error())
		return
	}
	result, err := solution.FindMolarityFromMass{Mass: mass, MolarMass: molarMass, Volume: vol}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}

// POST /api/solution/moles-from-molarity
// Body: { "molarity": "2.0", "volume": { "value": "0.250", "prefix": "none" } }
func HandleMolesFromMolarity(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Molarity string        `json:"molarity"`
		Volume   VolumeRequest `json:"volume"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	molarity, err := decimal.NewFromString(req.Molarity)
	if err != nil {
		writeError(w, http.StatusBadRequest, "invalid molarity: "+err.Error())
		return
	}
	vol, err := req.Volume.ToVolume()
	if err != nil {
		writeError(w, http.StatusBadRequest, "volume: "+err.Error())
		return
	}
	result, err := solution.FindMolesFromMolarity{Molarity: molarity, Volume: vol}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}

// POST /api/solution/molality
// Body: { "moles": "1.0", "solvent_mass": { "value": "1.250", "unit": "gram", "prefix": "kilo" } }
func HandleMolality(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Moles       string      `json:"moles"`
		SolventMass MassRequest `json:"solvent_mass"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	moles, err := decimal.NewFromString(req.Moles)
	if err != nil {
		writeError(w, http.StatusBadRequest, "invalid moles: "+err.Error())
		return
	}
	solventMass, err := req.SolventMass.ToMass()
	if err != nil {
		writeError(w, http.StatusBadRequest, "solvent_mass: "+err.Error())
		return
	}
	result, err := solution.FindMolality{Moles: moles, SolventMass: solventMass}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}

// POST /api/solution/dilution/concentration
// Body: { "initial_concentration": "6.0", "initial_volume": { "value": "0.050", "prefix": "none" }, "final_volume": { "value": "0.300", "prefix": "none" } }
func HandleDilutionConcentration(w http.ResponseWriter, r *http.Request) {
	var req struct {
		InitialConcentration string        `json:"initial_concentration"`
		InitialVolume        VolumeRequest `json:"initial_volume"`
		FinalVolume          VolumeRequest `json:"final_volume"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	c1, err := decimal.NewFromString(req.InitialConcentration)
	if err != nil {
		writeError(w, http.StatusBadRequest, "invalid initial_concentration: "+err.Error())
		return
	}
	v1, err := req.InitialVolume.ToVolume()
	if err != nil {
		writeError(w, http.StatusBadRequest, "initial_volume: "+err.Error())
		return
	}
	v2, err := req.FinalVolume.ToVolume()
	if err != nil {
		writeError(w, http.StatusBadRequest, "final_volume: "+err.Error())
		return
	}
	result, err := solution.DilutionFindFinalConcentration{
		InitialConcentration: c1,
		InitialVolume:        v1,
		FinalVolume:          v2,
	}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}

// POST /api/solution/dilution/volume
// Body: { "initial_concentration": "6.0", "initial_volume": { "value": "0.050", "prefix": "none" }, "final_concentration": "1.0" }
func HandleDilutionVolume(w http.ResponseWriter, r *http.Request) {
	var req struct {
		InitialConcentration string        `json:"initial_concentration"`
		InitialVolume        VolumeRequest `json:"initial_volume"`
		FinalConcentration   string        `json:"final_concentration"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	c1, err := decimal.NewFromString(req.InitialConcentration)
	if err != nil {
		writeError(w, http.StatusBadRequest, "invalid initial_concentration: "+err.Error())
		return
	}
	v1, err := req.InitialVolume.ToVolume()
	if err != nil {
		writeError(w, http.StatusBadRequest, "initial_volume: "+err.Error())
		return
	}
	c2, err := decimal.NewFromString(req.FinalConcentration)
	if err != nil {
		writeError(w, http.StatusBadRequest, "invalid final_concentration: "+err.Error())
		return
	}
	result, err := solution.DilutionFindFinalVolume{
		InitialConcentration: c1,
		InitialVolume:        v1,
		FinalConcentration:   c2,
	}.Calculate()
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	writeJSON(w, http.StatusOK, resultToResponse(result))
}
