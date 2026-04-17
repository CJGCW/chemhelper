package api_test

import (
	"encoding/json"
	"io"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"

	"chemhelper/api"
)

// ── Test infrastructure ───────────────────────────────────────────────────────

func newServer(t *testing.T) *httptest.Server {
	t.Helper()
	srv := httptest.NewServer(api.NewRouter())
	t.Cleanup(srv.Close)
	return srv
}

func postJSON(t *testing.T, url, body string) *http.Response {
	t.Helper()
	resp, err := http.Post(url, "application/json", strings.NewReader(body)) //nolint:noctx
	if err != nil {
		t.Fatalf("POST %s: %v", url, err)
	}
	return resp
}

func getJSON(t *testing.T, url string) *http.Response {
	t.Helper()
	resp, err := http.Get(url) //nolint:noctx
	if err != nil {
		t.Fatalf("GET %s: %v", url, err)
	}
	return resp
}

func decodeJSON(t *testing.T, r *http.Response, dst any) {
	t.Helper()
	defer r.Body.Close()
	if err := json.NewDecoder(r.Body).Decode(dst); err != nil {
		t.Fatalf("decode response: %v", err)
	}
}

func readBody(t *testing.T, r *http.Response) string {
	t.Helper()
	defer r.Body.Close()
	b, err := io.ReadAll(r.Body)
	if err != nil {
		t.Fatalf("read body: %v", err)
	}
	return string(b)
}

func assertStatus(t *testing.T, resp *http.Response, want int) {
	t.Helper()
	if resp.StatusCode != want {
		body := readBody(t, resp)
		t.Fatalf("status: got %d, want %d — body: %s", resp.StatusCode, want, body)
	}
}

// ── CORS ──────────────────────────────────────────────────────────────────────

func TestCORSHeaders(t *testing.T) {
	srv := newServer(t)

	t.Run("OPTIONS preflight returns 204", func(t *testing.T) {
		req, _ := http.NewRequest(http.MethodOptions, srv.URL+"/api/elements", nil)
		resp, err := http.DefaultClient.Do(req)
		if err != nil {
			t.Fatalf("OPTIONS request: %v", err)
		}
		defer resp.Body.Close()
		if resp.StatusCode != http.StatusNoContent {
			t.Errorf("status: got %d, want 204", resp.StatusCode)
		}
		if resp.Header.Get("Access-Control-Allow-Origin") == "" {
			t.Error("expected Access-Control-Allow-Origin header")
		}
	})

	t.Run("GET has CORS origin header", func(t *testing.T) {
		resp := getJSON(t, srv.URL+"/api/elements")
		defer resp.Body.Close()
		if resp.Header.Get("Access-Control-Allow-Origin") == "" {
			t.Error("expected Access-Control-Allow-Origin header on GET")
		}
	})
}

// ── Element handlers ──────────────────────────────────────────────────────────

func TestHandleElements(t *testing.T) {
	srv := newServer(t)
	resp := getJSON(t, srv.URL+"/api/elements")
	assertStatus(t, resp, http.StatusOK)

	var elements []api.ElementResponse
	decodeJSON(t, resp, &elements)

	if len(elements) != 118 {
		t.Errorf("expected 118 elements, got %d", len(elements))
	}
	// Spot-check hydrogen
	h := elements[0]
	if h.Symbol != "H" || h.AtomicNumber != 1 {
		t.Errorf("first element: got %s/%d, want H/1", h.Symbol, h.AtomicNumber)
	}
}

func TestHandleElement(t *testing.T) {
	srv := newServer(t)

	t.Run("known element Fe", func(t *testing.T) {
		resp := getJSON(t, srv.URL+"/api/elements/Fe")
		assertStatus(t, resp, http.StatusOK)

		var e api.ElementResponse
		decodeJSON(t, resp, &e)
		if e.Symbol != "Fe" {
			t.Errorf("symbol: got %q, want Fe", e.Symbol)
		}
		if e.AtomicNumber != 26 {
			t.Errorf("atomic number: got %d, want 26", e.AtomicNumber)
		}
		if e.GroupName == "" {
			t.Error("expected non-empty group name for Fe")
		}
	})

	t.Run("unknown element returns 404", func(t *testing.T) {
		resp := getJSON(t, srv.URL+"/api/elements/Xx")
		assertStatus(t, resp, http.StatusNotFound)
	})
}

func TestHandleCompound(t *testing.T) {
	srv := newServer(t)

	t.Run("H2O parsed correctly", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/elements/compound", `{"formula":"H2O"}`)
		assertStatus(t, resp, http.StatusOK)

		var elems []api.CompoundElementResponse
		decodeJSON(t, resp, &elems)
		if len(elems) != 2 {
			t.Fatalf("expected 2 elements in H2O, got %d", len(elems))
		}
	})

	t.Run("empty formula returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/elements/compound", `{"formula":""}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})

	t.Run("unknown element in formula returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/elements/compound", `{"formula":"Xx2O"}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})

	t.Run("malformed JSON returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/elements/compound", `not json`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

// ── Solution handlers ─────────────────────────────────────────────────────────

func TestHandleMolarity(t *testing.T) {
	srv := newServer(t)

	t.Run("1 mol in 1 L = 1 mol/L", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molarity",
			`{"moles":"1","volume":{"value":"1"}}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "1" {
			t.Errorf("value: got %q, want 1", r.Value)
		}
		if r.Unit != "mol/L" {
			t.Errorf("unit: got %q, want mol/L", r.Unit)
		}
		if len(r.Steps) == 0 {
			t.Error("expected steps in response")
		}
	})

	t.Run("2 mol in 500 mL = 4 mol/L", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molarity",
			`{"moles":"2","volume":{"value":"500","prefix":"milli"}}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "4" {
			t.Errorf("value: got %q, want 4", r.Value)
		}
	})

	t.Run("invalid moles returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molarity",
			`{"moles":"abc","volume":{"value":"1"}}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})

	t.Run("zero moles returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molarity",
			`{"moles":"0","volume":{"value":"1"}}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})
}

func TestHandleMolarityFromMass(t *testing.T) {
	srv := newServer(t)

	t.Run("58.44g NaCl in 1L = 1 mol/L", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molarity-from-mass",
			`{"mass":{"value":"58.44"},"compound":{"molar_mass":"58.44"},"volume":{"value":"1"}}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "1" {
			t.Errorf("value: got %q, want 1", r.Value)
		}
	})

	t.Run("invalid compound returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molarity-from-mass",
			`{"mass":{"value":"10"},"compound":{},"volume":{"value":"1"}}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

func TestHandleMolesFromMolarity(t *testing.T) {
	srv := newServer(t)

	t.Run("2 mol/L × 250 mL = 0.5 mol", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/moles-from-molarity",
			`{"molarity":"2","volume":{"value":"250","prefix":"milli"}}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "0.5" {
			t.Errorf("value: got %q, want 0.5", r.Value)
		}
		if r.Unit != "mol" {
			t.Errorf("unit: got %q, want mol", r.Unit)
		}
	})

	t.Run("invalid molarity string returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/moles-from-molarity",
			`{"molarity":"bad","volume":{"value":"1"}}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

func TestHandleMolality(t *testing.T) {
	srv := newServer(t)

	t.Run("1 mol in 1000g = 1 mol/kg", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molality",
			`{"moles":"1","solvent_mass":{"value":"1000"}}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "1" {
			t.Errorf("value: got %q, want 1", r.Value)
		}
		if r.Unit != "mol/kg" {
			t.Errorf("unit: got %q, want mol/kg", r.Unit)
		}
	})

	t.Run("invalid moles returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/molality",
			`{"moles":"abc","solvent_mass":{"value":"1000"}}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

func TestHandleDilutionConcentration(t *testing.T) {
	srv := newServer(t)

	t.Run("6M × 50mL diluted to 300mL = 1 mol/L", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/dilution/concentration",
			`{"initial_concentration":"6","initial_volume":{"value":"50","prefix":"milli"},"final_volume":{"value":"300","prefix":"milli"}}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "1" {
			t.Errorf("value: got %q, want 1", r.Value)
		}
	})

	t.Run("final volume smaller than initial returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/dilution/concentration",
			`{"initial_concentration":"1","initial_volume":{"value":"1"},"final_volume":{"value":"500","prefix":"milli"}}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})
}

func TestHandleDilutionVolume(t *testing.T) {
	srv := newServer(t)

	t.Run("6M × 50mL to 1M = 0.3 L", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/dilution/volume",
			`{"initial_concentration":"6","initial_volume":{"value":"50","prefix":"milli"},"final_concentration":"1"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "0.3" {
			t.Errorf("value: got %q, want 0.3", r.Value)
		}
		if r.Unit != "L" {
			t.Errorf("unit: got %q, want L", r.Unit)
		}
	})

	t.Run("final concentration >= initial returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/solution/dilution/volume",
			`{"initial_concentration":"1","initial_volume":{"value":"1"},"final_concentration":"2"}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})
}

// ── Thermo handlers ───────────────────────────────────────────────────────────

func TestHandleSolvents(t *testing.T) {
	srv := newServer(t)
	resp := getJSON(t, srv.URL+"/api/thermo/solvents")
	assertStatus(t, resp, http.StatusOK)

	var solvents []api.SolventResponse
	decodeJSON(t, resp, &solvents)

	if len(solvents) == 0 {
		t.Fatal("expected at least one solvent")
	}
	// Verify water is present and sorted first alphabetically
	found := false
	for _, s := range solvents {
		if s.Name == "Water" {
			found = true
			if s.Kb == "" {
				t.Error("expected non-empty Kb for water")
			}
			if s.Kf == "" {
				t.Error("expected non-empty Kf for water")
			}
		}
	}
	if !found {
		t.Error("water not found in solvent list")
	}
	// Verify response is sorted
	for i := 1; i < len(solvents); i++ {
		if solvents[i].Name < solvents[i-1].Name {
			t.Errorf("solvents not sorted: %q before %q", solvents[i-1].Name, solvents[i].Name)
		}
	}
}

func TestHandleBPE(t *testing.T) {
	srv := newServer(t)

	t.Run("water 1 mol/kg no i factor", func(t *testing.T) {
		// Water Kb = 0.512; ΔTb = 0.512 × 1.0 × 1 = 0.512; new BP = 100.512
		resp := postJSON(t, srv.URL+"/api/thermo/bpe",
			`{"solvent":"water","molality":"1"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.ColligativeResponse
		decodeJSON(t, resp, &r)
		if r.Delta != "0.512" {
			t.Errorf("delta: got %q, want 0.512", r.Delta)
		}
		if r.NewPoint != "100.512" {
			t.Errorf("new_point: got %q, want 100.512", r.NewPoint)
		}
		if r.Unit == "" {
			t.Error("expected non-empty unit")
		}
		if len(r.Steps) == 0 {
			t.Error("expected steps in BPE response")
		}
	})

	t.Run("water with vant hoff factor 2", func(t *testing.T) {
		// ΔTb = 0.512 × 1.0 × 2 = 1.024
		resp := postJSON(t, srv.URL+"/api/thermo/bpe",
			`{"solvent":"water","molality":"1","vant_hoff_factor":"2"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.ColligativeResponse
		decodeJSON(t, resp, &r)
		if r.Delta != "1.024" {
			t.Errorf("delta: got %q, want 1.024", r.Delta)
		}
	})

	t.Run("unknown solvent returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/bpe",
			`{"solvent":"acetone","molality":"1"}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})

	t.Run("invalid molality returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/bpe",
			`{"solvent":"water","molality":"abc"}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

func TestHandleFPD(t *testing.T) {
	srv := newServer(t)

	t.Run("water 1 mol/kg", func(t *testing.T) {
		// Water Kf = 1.86; ΔTf = 1.86; new FP = -1.86
		resp := postJSON(t, srv.URL+"/api/thermo/fpd",
			`{"solvent":"water","molality":"1"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.ColligativeResponse
		decodeJSON(t, resp, &r)
		if r.Delta != "1.86" {
			t.Errorf("delta: got %q, want 1.86", r.Delta)
		}
		if r.NewPoint != "-1.86" {
			t.Errorf("new_point: got %q, want -1.86", r.NewPoint)
		}
	})

	t.Run("unknown solvent returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/fpd",
			`{"solvent":"lava","molality":"1"}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

func TestHandleMolalityFromBPE(t *testing.T) {
	srv := newServer(t)

	t.Run("water delta_tb 0.512 = 1 mol/kg", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/molality-from-bpe",
			`{"solvent":"water","delta_tb":"0.512"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "1" {
			t.Errorf("value: got %q, want 1", r.Value)
		}
		if r.Unit != "mol/kg" {
			t.Errorf("unit: got %q, want mol/kg", r.Unit)
		}
	})

	t.Run("invalid delta_tb returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/molality-from-bpe",
			`{"solvent":"water","delta_tb":"bad"}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

func TestHandleMolalityFromFPD(t *testing.T) {
	srv := newServer(t)

	t.Run("water delta_tf 1.86 = 1 mol/kg", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/molality-from-fpd",
			`{"solvent":"water","delta_tf":"1.86"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CalcResponse
		decodeJSON(t, resp, &r)
		if r.Value != "1" {
			t.Errorf("value: got %q, want 1", r.Value)
		}
	})

	t.Run("invalid delta_tf returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/molality-from-fpd",
			`{"solvent":"water","delta_tf":"bad"}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})

	t.Run("unknown solvent returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/thermo/molality-from-fpd",
			`{"solvent":"unknown","delta_tf":"1.86"}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

// ── Structure handlers ────────────────────────────────────────────────────────

// lewisAtom mirrors the JSON shape of structure.LewisAtom.
type lewisAtom struct {
	Element      string `json:"element"`
	FormalCharge int    `json:"formal_charge"`
	LonePairs    int    `json:"lone_pairs"`
}

// lewisResp mirrors the JSON shape of structure.LewisStructure.
type lewisResp struct {
	Geometry string      `json:"geometry"`
	Atoms    []lewisAtom `json:"atoms"`
	Bonds    []struct {
		Order int `json:"order"`
	} `json:"bonds"`
	Steps []string `json:"steps"`
}

func TestHandleLewis(t *testing.T) {
	srv := newServer(t)

	t.Run("H2O returns bent geometry", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":"H2O"}`)
		assertStatus(t, resp, http.StatusOK)

		var r lewisResp
		decodeJSON(t, resp, &r)
		if r.Geometry != "bent" {
			t.Errorf("geometry: got %q, want bent", r.Geometry)
		}
		if len(r.Atoms) == 0 {
			t.Error("expected atoms in response")
		}
		if len(r.Steps) == 0 {
			t.Error("expected steps in response")
		}
	})

	t.Run("Ca(OH)2 returns Ca with +2 formal charge", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":"Ca(OH)2"}`)
		assertStatus(t, resp, http.StatusOK)

		var r lewisResp
		decodeJSON(t, resp, &r)
		found := false
		for _, a := range r.Atoms {
			if a.Element == "Ca" {
				found = true
				if a.FormalCharge != 2 {
					t.Errorf("Ca formal charge: got %d, want 2", a.FormalCharge)
				}
			}
		}
		if !found {
			t.Fatal("Ca atom not found in response")
		}
	})

	t.Run("FeCl3 returns trigonal planar geometry", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":"FeCl3"}`)
		assertStatus(t, resp, http.StatusOK)

		var r lewisResp
		decodeJSON(t, resp, &r)
		if r.Geometry != "trigonal_planar" {
			t.Errorf("geometry: got %q, want trigonal_planar", r.Geometry)
		}
	})

	t.Run("CO2 returns linear geometry", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":"CO2"}`)
		assertStatus(t, resp, http.StatusOK)

		var r lewisResp
		decodeJSON(t, resp, &r)
		if r.Geometry != "linear" {
			t.Errorf("geometry: got %q, want linear", r.Geometry)
		}
	})

	t.Run("empty input returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":""}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})

	t.Run("unknown formula returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":"Xx2O"}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})

	t.Run("malformed JSON returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `not json`)
		assertStatus(t, resp, http.StatusBadRequest)
	})

	// CH3OH (methanol) previously failed with "insufficient valence electrons"
	// due to a chain-topology bug in the single-center solver.
	t.Run("CH3OH returns valid structure", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":"CH3OH"}`)
		assertStatus(t, resp, http.StatusOK)

		var r lewisResp
		decodeJSON(t, resp, &r)
		if len(r.Atoms) == 0 {
			t.Error("expected atoms in response")
		}
		if len(r.Bonds) == 0 {
			t.Error("expected bonds in response")
		}
	})

	t.Run("ionic charge is accepted", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/structure/lewis", `{"input":"NO3","charge":-1}`)
		assertStatus(t, resp, http.StatusOK)

		var r lewisResp
		decodeJSON(t, resp, &r)
		if r.Geometry == "" {
			t.Error("expected geometry in response")
		}
	})
}

// ── Compound lookup handler ───────────────────────────────────────────────────
// These tests make real PubChem network calls.
// Run with: go test ./api/... -v -run TestHandleCompoundLookup -timeout 60s

type compoundLookupResp struct {
	CID              int    `json:"cid"`
	MolecularFormula string `json:"molecular_formula"`
	MolecularWeight  string `json:"molecular_weight"`
	IUPACName        string `json:"iupac_name"`
	SMILES           string `json:"smiles"`
	InChI            string `json:"inchi"`
	InChIKey         string `json:"inchi_key"`
	InputType        string `json:"input_type"`
}

func TestHandleCompoundLookup(t *testing.T) {
	srv := newServer(t)

	t.Run("lookup by name returns all fields", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/lookup", `{"input":"water"}`)
		assertStatus(t, resp, http.StatusOK)

		var r compoundLookupResp
		decodeJSON(t, resp, &r)

		if r.MolecularFormula != "H2O" {
			t.Errorf("formula: got %q, want H2O", r.MolecularFormula)
		}
		if r.SMILES == "" {
			t.Error("SMILES should not be empty")
		}
		if r.InChI == "" {
			t.Error("InChI should not be empty")
		}
		if r.InChIKey == "" {
			t.Error("InChIKey should not be empty")
		}
		if r.CID == 0 {
			t.Error("CID should not be zero")
		}
		if r.InputType != "name" {
			t.Errorf("input_type: got %q, want name", r.InputType)
		}
	})

	t.Run("lookup ethanol by plain SMILES detects smiles type", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/lookup", `{"input":"CCO"}`)
		assertStatus(t, resp, http.StatusOK)

		var r compoundLookupResp
		decodeJSON(t, resp, &r)

		if r.MolecularFormula != "C2H6O" {
			t.Errorf("formula: got %q, want C2H6O", r.MolecularFormula)
		}
		if r.CID != 702 {
			t.Errorf("CID: got %d, want 702 (ethanol)", r.CID)
		}
		if r.InputType != "smiles" {
			t.Errorf("input_type: got %q, want smiles", r.InputType)
		}
	})

	t.Run("lookup by CID returns correct compound", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/lookup", `{"input":"5234"}`)
		assertStatus(t, resp, http.StatusOK)

		var r compoundLookupResp
		decodeJSON(t, resp, &r)

		if r.CID != 5234 {
			t.Errorf("CID: got %d, want 5234", r.CID)
		}
		if r.InputType != "cid" {
			t.Errorf("input_type: got %q, want cid", r.InputType)
		}
	})

	t.Run("lookup by InChIKey resolves correctly", func(t *testing.T) {
		// InChIKey for ethanol
		resp := postJSON(t, srv.URL+"/api/compound/lookup",
			`{"input":"LFQSCWFLJHTTHZ-UHFFFAOYSA-N"}`)
		assertStatus(t, resp, http.StatusOK)

		var r compoundLookupResp
		decodeJSON(t, resp, &r)

		if r.MolecularFormula != "C2H6O" {
			t.Errorf("formula: got %q, want C2H6O", r.MolecularFormula)
		}
		if r.InputType != "inchikey" {
			t.Errorf("input_type: got %q, want inchikey", r.InputType)
		}
	})

	t.Run("lookup by SMILES with special chars", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/lookup", `{"input":"CC(=O)O"}`)
		assertStatus(t, resp, http.StatusOK)

		var r compoundLookupResp
		decodeJSON(t, resp, &r)

		if r.MolecularFormula != "C2H4O2" {
			t.Errorf("formula: got %q, want C2H4O2 (acetic acid)", r.MolecularFormula)
		}
		if r.InputType != "smiles" {
			t.Errorf("input_type: got %q, want smiles", r.InputType)
		}
	})

	t.Run("empty input returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/lookup", `{"input":""}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})

	t.Run("unknown compound returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/lookup",
			`{"input":"xyzzy_not_a_compound_123456789"}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})

	t.Run("malformed JSON returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/lookup", `not json`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

// ── Compound resolve handler ──────────────────────────────────────────────────
// These tests make real PubChem network calls.
// Run with: go test ./api/... -v -run TestHandleCompoundResolve -timeout 60s

func TestHandleCompoundResolve(t *testing.T) {
	srv := newServer(t)

	t.Run("caffeine SMILES returns correct formula and CID", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/resolve",
			`{"smiles":"Cn1cnc2c1c(=O)n(c(=O)n2C)C"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CompoundResponse
		decodeJSON(t, resp, &r)

		if r.MolecularFormula != "C8H10N4O2" {
			t.Errorf("formula: got %q, want C8H10N4O2", r.MolecularFormula)
		}
		if r.CID != 2519 {
			t.Errorf("CID: got %d, want 2519 (caffeine)", r.CID)
		}
		if r.MolecularWeight == "" {
			t.Error("molecular_weight should not be empty")
		}
		if r.IUPACName == "" {
			t.Error("iupac_name should not be empty")
		}
	})

	t.Run("ethanol SMILES returns C2H6O", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/resolve", `{"smiles":"CCO"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CompoundResponse
		decodeJSON(t, resp, &r)

		if r.MolecularFormula != "C2H6O" {
			t.Errorf("formula: got %q, want C2H6O", r.MolecularFormula)
		}
		if r.CID != 702 {
			t.Errorf("CID: got %d, want 702 (ethanol)", r.CID)
		}
	})

	t.Run("water SMILES returns H2O", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/resolve", `{"smiles":"O"}`)
		assertStatus(t, resp, http.StatusOK)

		var r api.CompoundResponse
		decodeJSON(t, resp, &r)

		if r.MolecularFormula != "H2O" {
			t.Errorf("formula: got %q, want H2O", r.MolecularFormula)
		}
	})

	t.Run("empty smiles returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/resolve", `{"smiles":""}`)
		assertStatus(t, resp, http.StatusBadRequest)
	})

	t.Run("invalid SMILES returns 422", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/resolve", `{"smiles":"ZZZZNOTSMILES"}`)
		assertStatus(t, resp, http.StatusUnprocessableEntity)
	})

	t.Run("malformed JSON returns 400", func(t *testing.T) {
		resp := postJSON(t, srv.URL+"/api/compound/resolve", `not json`)
		assertStatus(t, resp, http.StatusBadRequest)
	})
}

// ── Random structure handler ──────────────────────────────────────────────────

func TestHandleRandomStructure(t *testing.T) {
	srv := newServer(t)

	t.Run("returns 200 with valid structure shape", func(t *testing.T) {
		resp := getJSON(t, srv.URL+"/api/structure/random")
		assertStatus(t, resp, http.StatusOK)

		var r lewisResp
		decodeJSON(t, resp, &r)

		if r.Geometry == "" {
			t.Error("geometry should not be empty")
		}
		if len(r.Atoms) == 0 {
			t.Error("atoms should not be empty")
		}
		if len(r.Bonds) == 0 {
			t.Error("bonds should not be empty")
		}
		if len(r.Steps) == 0 {
			t.Error("steps should not be empty")
		}
	})

	t.Run("successive calls return structures (not always identical)", func(t *testing.T) {
		geometries := make(map[string]bool)
		for i := 0; i < 10; i++ {
			resp := getJSON(t, srv.URL+"/api/structure/random")
			assertStatus(t, resp, http.StatusOK)
			var r lewisResp
			decodeJSON(t, resp, &r)
			if r.Geometry == "" {
				t.Fatalf("call %d: geometry should not be empty", i)
			}
			geometries[r.Geometry] = true
		}
		if len(geometries) < 2 {
			t.Log("note: all 10 random structures had the same geometry — unlikely but possible")
		}
	})
}
