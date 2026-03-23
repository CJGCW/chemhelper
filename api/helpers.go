package api

import (
	"encoding/json"
	"net/http"

	"chemhelper/calc"
)

const maxBodyBytes = 1 << 20 // 1 MB

// writeJSON serialises v to JSON and writes it with the given status code.
func writeJSON(w http.ResponseWriter, status int, v any) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(status)
	json.NewEncoder(w).Encode(v)
}

// writeError writes a JSON error response.
func writeError(w http.ResponseWriter, status int, msg string) {
	writeJSON(w, status, ErrorResponse{Error: msg})
}

// decode reads and decodes a JSON request body into dst.
// Request bodies are capped at 1 MB.
func decode(r *http.Request, dst any) error {
	r.Body = http.MaxBytesReader(nil, r.Body, maxBodyBytes)
	return json.NewDecoder(r.Body).Decode(dst)
}

// resultToResponse converts a calc.Result to a CalcResponse.
func resultToResponse(r calc.Result) CalcResponse {
	return CalcResponse{
		Value:   r.Value.String(),
		Unit:    r.Unit,
		SigFigs: r.SigFigs,
		Steps:   r.Steps,
	}
}
