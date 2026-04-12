package structure

import (
	"fmt"
	"math/rand"
)

// randomCombo is a valid symmetric AB_n molecule specification.
type randomCombo struct {
	Center   string
	Terminal string
	N        int
}

// validCombos is built once at startup from chemistry rules.
var validCombos []randomCombo

func init() {
	validCombos = buildValidCombos()
}

// buildValidCombos enumerates every symmetric AB_n molecule where:
//   - totalVE is even (no radicals)
//   - remaining electrons on center are non-negative and even
//   - electron groups (bonds + lone pairs on center) are 2–6
//   - period-2 centers cannot exceed 4 electron groups (no expanded octet)
func buildValidCombos() []randomCombo {
	centers   := []string{"Be", "B", "C", "N", "O", "Si", "P", "S", "Cl", "As", "Se", "Br", "I", "Xe", "Te"}
	terminals := []string{"H", "F", "Cl", "Br"}

	var combos []randomCombo
	for _, center := range centers {
		vc        := elemValence(center)
		period    := elemPeriod(center)
		canExpand := period >= 3
		maxN := 4
		if canExpand {
			maxN = 6
		}

		for _, terminal := range terminals {
			if center == terminal {
				continue
			}
			// The backend picks the least-electronegative heavy atom as centre.
			// If the terminal is less electronegative than the centre, the backend
			// would misidentify the terminal as centre (e.g. OBr2 → Br picked as
			// centre, two-centre solver, wrong geometry). Skip those combos.
			// H is always treated as a terminal by the backend regardless of EN.
			if terminal != "H" && elemEN(terminal) < elemEN(center) {
				continue
			}
			vt  := elemValence(terminal)
			tlp := 6 // lone electrons placed on terminal to complete its octet
			if terminal == "H" {
				tlp = 0 // hydrogen has no lone pairs (duet rule)
			}

			for n := 2; n <= maxN; n++ {
				totalVE := vc + n*vt
				if totalVE%2 != 0 {
					continue // radical — skip
				}
				onCenter := totalVE - 2*n - n*tlp
				if onCenter < 0 || onCenter%2 != 0 {
					continue
				}
				lpCenter := onCenter / 2
				groups   := n + lpCenter
				if groups < 2 || groups > 6 {
					continue
				}
				if !canExpand && groups > 4 {
					continue // period-2 elements cannot expand octet
				}
				combos = append(combos, randomCombo{center, terminal, n})
			}
		}
	}
	return combos
}

// GenerateRandom picks a random valid AB_n molecule and returns its Lewis structure.
func GenerateRandom() (*LewisStructure, error) {
	if len(validCombos) == 0 {
		return nil, fmt.Errorf("no valid combinations available")
	}
	c := validCombos[rand.Intn(len(validCombos))]
	formula := fmt.Sprintf("%s%s%d", c.Center, c.Terminal, c.N)
	return generateFromParsed(formula, 0)
}
