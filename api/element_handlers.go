package api

import (
	"net/http"

	"chemhelper/element"

	"github.com/go-chi/chi/v5"
)

// pt is the single shared periodic table instance.
var pt = element.NewPeriodicTable()

// GET /api/elements
func HandleElements(w http.ResponseWriter, r *http.Request) {
	resp := make([]ElementResponse, 0, len(pt.Elements))
	for _, e := range pt.Elements {
		group, _ := e.GetGroup()
		resp = append(resp, ElementResponse{
			AtomicNumber:      e.AtomicNumber,
			Symbol:            e.Symbol,
			Name:              e.Name,
			AtomicWeight:      e.AtomicWeight.String(),
			Electronegativity: e.Electronegativity,
			VanDerWaalsRadius: e.VanDerWaalsRadius,
			Group:             e.Group,
			Period:            e.Period,
			GroupName:         group,
		})
	}
	writeJSON(w, http.StatusOK, resp)
}

// GET /api/elements/{symbol}
func HandleElement(w http.ResponseWriter, r *http.Request) {
	symbol := chi.URLParam(r, "symbol")
	e, found := pt.FindElementBySymbol(symbol)
	if !found {
		writeError(w, http.StatusNotFound, "element not found: "+symbol)
		return
	}
	group, _ := e.GetGroup()
	writeJSON(w, http.StatusOK, ElementResponse{
		AtomicNumber:      e.AtomicNumber,
		Symbol:            e.Symbol,
		Name:              e.Name,
		AtomicWeight:      e.AtomicWeight.String(),
		Electronegativity: e.Electronegativity,
		VanDerWaalsRadius: e.VanDerWaalsRadius,
		Group:             e.Group,
		Period:            e.Period,
		GroupName:         group,
	})
}

// POST /api/elements/compound
// Body: { "formula": "H2O" }
// Returns the list of elements and their atom counts in the formula.
func HandleCompound(w http.ResponseWriter, r *http.Request) {
	var req struct {
		Formula string `json:"formula"`
	}
	if err := decode(r, &req); err != nil {
		writeError(w, http.StatusBadRequest, "invalid request body: "+err.Error())
		return
	}
	elements, err := element.ParseCompoundElements(req.Formula, pt)
	if err != nil {
		writeError(w, http.StatusUnprocessableEntity, err.Error())
		return
	}
	resp := make([]CompoundElementResponse, 0, len(elements))
	for _, em := range elements {
		resp = append(resp, CompoundElementResponse{
			Symbol: em.Element.Symbol,
			Name:   em.Element.Name,
			Moles:  em.Moles.String(),
		})
	}
	writeJSON(w, http.StatusOK, resp)
}
