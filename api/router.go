package api

import (
	"net/http"

	"github.com/go-chi/chi/v5"
	"github.com/go-chi/chi/v5/middleware"
)

// NewRouter builds and returns the fully wired chi router.
func NewRouter() http.Handler {
	r := chi.NewRouter()

	r.Use(middleware.Logger)
	r.Use(middleware.Recoverer)
	r.Use(corsMiddleware)

	r.Route("/api", func(r chi.Router) {
		// Compound (SMILES resolution)
		r.Post("/compound/resolve", HandleCompoundResolve)

		// Elements
		r.Get("/elements", HandleElements)
		r.Get("/elements/{symbol}", HandleElement)
		r.Post("/elements/compound", HandleCompound)

		// Solution
		r.Route("/solution", func(r chi.Router) {
			r.Post("/molarity", HandleMolarity)
			r.Post("/molarity-from-mass", HandleMolarityFromMass)
			r.Post("/moles-from-molarity", HandleMolesFromMolarity)
			r.Post("/molality", HandleMolality)
			r.Post("/dilution/concentration", HandleDilutionConcentration)
			r.Post("/dilution/volume", HandleDilutionVolume)
		})

		// Thermo
		r.Route("/thermo", func(r chi.Router) {
			r.Get("/solvents", HandleSolvents)
			r.Post("/bpe", HandleBPE)
			r.Post("/fpd", HandleFPD)
			r.Post("/molality-from-bpe", HandleMolalityFromBPE)
			r.Post("/molality-from-fpd", HandleMolalityFromFPD)
		})
	})

	return r
}

// corsMiddleware allows requests from the React dev server.
func corsMiddleware(next http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Access-Control-Allow-Origin", "http://localhost:5173")
		w.Header().Set("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
		w.Header().Set("Access-Control-Allow-Headers", "Content-Type")
		if r.Method == http.MethodOptions {
			w.WriteHeader(http.StatusNoContent)
			return
		}
		next.ServeHTTP(w, r)
	})
}
