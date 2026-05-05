package api

import (
	"net/http"

	"chemhelper/organic"
)

// POST /api/organic/name
// Body: { "carbons": [ { "position": 2, "substituents": [{ "name": "methyl", "position": 1 }], "bond_order": 1 }, ... ] }
func HandleOrganicName(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Carbons []struct {
			Position     int `json:"position"`
			BondOrder    int `json:"bond_order"`
			Substituents []struct {
				Name     string `json:"name"`
				Position int    `json:"position"`
			} `json:"substituents"`
		} `json:"carbons"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if len(req.Carbons) == 0 {
		writeError(w, http.StatusUnprocessableEntity, "carbons array must not be empty")
		return
	}

	carbons := make([]organic.Carbon, len(req.Carbons))
	for i, c := range req.Carbons {
		carbons[i] = organic.Carbon{
			Position:  c.Position,
			BondOrder: c.BondOrder,
		}
		for _, s := range c.Substituents {
			carbons[i].Substituents = append(carbons[i].Substituents,
				organic.Substituent{Name: s.Name, Position: s.Position})
		}
	}

	result, err := organic.NameAlkane(organic.Chain{Carbons: carbons})
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}

	writeJSON(w, http.StatusOK, map[string]any{
		"name":  result.Name,
		"steps": result.Steps,
	})
}

// POST /api/organic/validate-name
// Body: { "student": "2-methylbutane", "correct": "2-methylbutane" }
func HandleOrganicValidateName(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Student string `json:"student"`
		Correct string `json:"correct"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	if req.Correct == "" {
		writeError(w, http.StatusUnprocessableEntity, "correct name must not be empty")
		return
	}

	writeJSON(w, http.StatusOK, map[string]bool{
		"valid": organic.ValidateName(req.Student, req.Correct),
	})
}
