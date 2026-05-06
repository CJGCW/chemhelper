package api

import (
	"net/http"

	"chemhelper/structure/render"
)

type structureRenderRequest struct {
	SMILES        string `json:"smiles"`
	Width         int    `json:"width"`
	Height        int    `json:"height"`
	ShowLonePairs bool   `json:"showLonePairs"`
}

type structureRenderResponse struct {
	SVG string `json:"svg"`
}

// HandleStructureRender handles POST /api/structure/render.
// Request:  { "smiles": "CC(=O)O", "width": 240, "height": 180 }
// Response: { "svg": "<svg ...>...</svg>" }
func HandleStructureRender(w http.ResponseWriter, r *http.Request) {
	var req structureRenderRequest
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusUnprocessableEntity, "invalid request body")
		return
	}
	if req.SMILES == "" {
		writeError(w, http.StatusUnprocessableEntity, "smiles is required")
		return
	}
	if req.Width <= 0 {
		req.Width = 240
	}
	if req.Height <= 0 {
		req.Height = 180
	}

	svg, err := render.Render(req.SMILES, req.Width, req.Height, render.RenderOptions{
		ShowLonePairs: req.ShowLonePairs,
	})
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}

	writeJSON(w, http.StatusOK, structureRenderResponse{SVG: svg})
}
