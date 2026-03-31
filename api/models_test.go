package api

import (
	"testing"
)

// ── parsePrefix ───────────────────────────────────────────────────────────────

func TestParsePrefix(t *testing.T) {
	cases := []struct {
		input   string
		wantErr bool
	}{
		{"", false},
		{"none", false},
		{"kilo", false},
		{"hecto", false},
		{"deca", false},
		{"deci", false},
		{"centi", false},
		{"milli", false},
		{"micro", false},
		{"MILLI", false}, // case-insensitive
		{"nano", true},
		{"mega", true},
		{"unknown", true},
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			_, err := parsePrefix(c.input)
			if c.wantErr && err == nil {
				t.Errorf("expected error for prefix %q but got none", c.input)
			}
			if !c.wantErr && err != nil {
				t.Errorf("unexpected error for prefix %q: %v", c.input, err)
			}
		})
	}
}

// ── parseMassUnit ─────────────────────────────────────────────────────────────

func TestParseMassUnit(t *testing.T) {
	cases := []struct {
		input   string
		wantErr bool
	}{
		{"", false},
		{"gram", false},
		{"GRAM", false},
		{"pound", false},
		{"ounce", false},
		{"kilogram", true},
		{"mg", true},
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			_, err := parseMassUnit(c.input)
			if c.wantErr && err == nil {
				t.Errorf("expected error for unit %q but got none", c.input)
			}
			if !c.wantErr && err != nil {
				t.Errorf("unexpected error for unit %q: %v", c.input, err)
			}
		})
	}
}

// ── MassRequest.ToMass ────────────────────────────────────────────────────────

func TestMassRequestToMass(t *testing.T) {
	cases := []struct {
		name    string
		req     MassRequest
		wantErr bool
	}{
		{"default gram", MassRequest{Value: "10"}, false},
		{"explicit gram", MassRequest{Value: "58.44", Unit: "gram"}, false},
		{"kilogram prefix", MassRequest{Value: "1.5", Unit: "gram", Prefix: "kilo"}, false},
		{"milligram", MassRequest{Value: "500", Unit: "gram", Prefix: "milli"}, false},
		{"pound", MassRequest{Value: "2.2", Unit: "pound"}, false},
		{"invalid value", MassRequest{Value: "abc"}, true},
		{"invalid unit", MassRequest{Value: "10", Unit: "ton"}, true},
		{"invalid prefix", MassRequest{Value: "10", Prefix: "mega"}, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := c.req.ToMass()
			if c.wantErr && err == nil {
				t.Errorf("expected error but got none")
			}
			if !c.wantErr && err != nil {
				t.Errorf("unexpected error: %v", err)
			}
		})
	}
}

// ── VolumeRequest.ToVolume ────────────────────────────────────────────────────

func TestVolumeRequestToVolume(t *testing.T) {
	cases := []struct {
		name    string
		req     VolumeRequest
		wantErr bool
	}{
		{"default litre", VolumeRequest{Value: "1"}, false},
		{"millilitre", VolumeRequest{Value: "500", Prefix: "milli"}, false},
		{"microlitre", VolumeRequest{Value: "250", Prefix: "micro"}, false},
		{"kilolitre", VolumeRequest{Value: "2", Prefix: "kilo"}, false},
		{"invalid value", VolumeRequest{Value: "abc"}, true},
		{"invalid prefix", VolumeRequest{Value: "1", Prefix: "giga"}, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := c.req.ToVolume()
			if c.wantErr && err == nil {
				t.Errorf("expected error but got none")
			}
			if !c.wantErr && err != nil {
				t.Errorf("unexpected error: %v", err)
			}
		})
	}
}

// ── CompoundRequest.ResolveMolarMass ──────────────────────────────────────────

func TestCompoundRequestResolveMolarMass(t *testing.T) {
	t.Run("molar_mass path", func(t *testing.T) {
		c := CompoundRequest{MolarMass: "58.44"}
		mm, err := c.ResolveMolarMass()
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if mm.String() != "58.44" {
			t.Errorf("got %s, want 58.44", mm)
		}
	})
	t.Run("neither set", func(t *testing.T) {
		_, err := CompoundRequest{}.ResolveMolarMass()
		if err == nil {
			t.Error("expected error when neither smiles nor molar_mass is set")
		}
	})
	t.Run("both set", func(t *testing.T) {
		_, err := CompoundRequest{SMILES: "O", MolarMass: "18.015"}.ResolveMolarMass()
		if err == nil {
			t.Error("expected error when both smiles and molar_mass are set")
		}
	})
	t.Run("invalid molar_mass string", func(t *testing.T) {
		_, err := CompoundRequest{MolarMass: "not-a-number"}.ResolveMolarMass()
		if err == nil {
			t.Error("expected error for non-numeric molar_mass")
		}
	})
}

// ── parseVantHoff ─────────────────────────────────────────────────────────────

func TestParseVantHoff(t *testing.T) {
	cases := []struct {
		input   string
		want    string
		wantErr bool
	}{
		{"", "0", false},
		{"1", "1", false},
		{"2.5", "2.5", false},
		{"abc", "", true},
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			v, err := parseVantHoff(c.input)
			if c.wantErr {
				if err == nil {
					t.Errorf("expected error for %q but got none", c.input)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if v.String() != c.want {
				t.Errorf("got %s, want %s", v, c.want)
			}
		})
	}
}

// ── isSMILES ──────────────────────────────────────────────────────────────────

func TestIsSMILES(t *testing.T) {
	cases := []struct {
		input string
		want  bool
	}{
		{"H2O", false},
		{"Ca(OH)2", false},
		{"FeCl3", false},
		{"NaCl", false},
		{"O=C=O", true},        // double bond
		{"[Na+].[Cl-]", true},  // brackets and dot
		{"C#N", true},           // triple bond
		{"[Ca+2].[OH-]", true},  // brackets
		{"CC(=O)O", true},       // equals sign
		{"C@H", true},           // chirality
		{"C/C=C/C", true},       // stereo slash
	}
	for _, c := range cases {
		t.Run(c.input, func(t *testing.T) {
			if got := isSMILES(c.input); got != c.want {
				t.Errorf("isSMILES(%q) = %v, want %v", c.input, got, c.want)
			}
		})
	}
}
