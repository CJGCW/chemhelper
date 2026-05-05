package organic_test

import (
	"testing"

	"chemhelper/organic"
)

func makeChain(n int, subsByPos map[int][]string) organic.Chain {
	carbons := make([]organic.Carbon, n)
	for i := range carbons {
		pos := i + 1
		carbons[i] = organic.Carbon{Position: pos, BondOrder: 1}
		for _, name := range subsByPos[pos] {
			carbons[i].Substituents = append(carbons[i].Substituents,
				organic.Substituent{Name: name, Position: 1})
		}
	}
	return organic.Chain{Carbons: carbons}
}

func TestNameAlkane(t *testing.T) {
	tests := []struct {
		name  string
		chain organic.Chain
		want  string
	}{
		{
			name:  "2-methylbutane",
			chain: makeChain(4, map[int][]string{2: {"methyl"}}),
			want:  "2-methylbutane",
		},
		{
			name:  "2,3-dimethylpentane",
			chain: makeChain(5, map[int][]string{2: {"methyl"}, 3: {"methyl"}}),
			want:  "2,3-dimethylpentane",
		},
		{
			name:  "4-ethyl-2-methylhexane",
			chain: makeChain(6, map[int][]string{2: {"methyl"}, 4: {"ethyl"}}),
			want:  "4-ethyl-2-methylhexane",
		},
		{
			name:  "2,2,4-trimethylpentane",
			chain: makeChain(5, map[int][]string{2: {"methyl", "methyl"}, 4: {"methyl"}}),
			want:  "2,2,4-trimethylpentane",
		},
		{
			name:  "unsubstituted hexane",
			chain: makeChain(6, nil),
			want:  "hexane",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			result, err := organic.NameAlkane(tc.chain)
			if err != nil {
				t.Fatalf("NameAlkane error: %v", err)
			}
			if result.Name != tc.want {
				t.Errorf("got %q, want %q", result.Name, tc.want)
			}
			if len(result.Steps) == 0 {
				t.Error("expected non-empty steps")
			}
		})
	}
}

func TestNameAlkaneErrors(t *testing.T) {
	_, err := organic.NameAlkane(organic.Chain{})
	if err == nil {
		t.Error("expected error for empty chain, got nil")
	}

	_, err = organic.NameAlkane(organic.Chain{Carbons: make([]organic.Carbon, 15)})
	if err == nil {
		t.Error("expected error for chain length 15")
	}
}

func TestValidateName(t *testing.T) {
	tests := []struct {
		student string
		correct string
		want    bool
	}{
		{"2-methylbutane", "2-methylbutane", true},
		{"2-Methylbutane", "2-methylbutane", true},
		{"2 - methylbutane", "2-methylbutane", true},
		{"3-methylbutane", "2-methylbutane", false},
		{"2,3-dimethylpentane", "2,3-dimethylpentane", true},
		{"4-ethyl-2-methylhexane", "4-ethyl-2-methylhexane", true},
		{"2,2,4-trimethylpentane", "2,2,4-trimethylpentane", true},
		{"2,2,4-Trimethylpentane", "2,2,4-trimethylpentane", true},
		{"wrong", "2-methylbutane", false},
	}
	for _, tc := range tests {
		got := organic.ValidateName(tc.student, tc.correct)
		if got != tc.want {
			t.Errorf("ValidateName(%q, %q) = %v, want %v", tc.student, tc.correct, got, tc.want)
		}
	}
}
