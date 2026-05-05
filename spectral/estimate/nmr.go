package estimate

import "sort"

// estimateHNMR returns estimated ¹H NMR peaks.
// Chemically equivalent H's are grouped into a single peak with summed integration.
func estimateHNMR(m *molecule) []Peak {
	type hGroup struct {
		ppm         float64
		splitting   string
		integration int
		label       string
	}

	sigToGroup := map[string]*hGroup{}

	for i := range m.atoms {
		a := &m.atoms[i]
		if a.hCount == 0 {
			continue
		}

		ppm, splitting, label := classifyHEnv(m, i)
		if label == "" {
			continue
		}

		sig := atomSig(m, i, 3, map[int]bool{})
		// Append label to sig so different environments at same ppm don't merge
		key := sig + "|" + label

		if g, ok := sigToGroup[key]; ok {
			g.integration += a.hCount
		} else {
			sigToGroup[key] = &hGroup{
				ppm:         ppm,
				splitting:   splitting,
				integration: a.hCount,
				label:       label,
			}
		}
	}

	peaks := make([]Peak, 0, len(sigToGroup))
	for _, g := range sigToGroup {
		peaks = append(peaks, Peak{
			X:           g.ppm,
			Y:           float64(g.integration) * 15,
			Label:       g.label,
			Width:       0.08,
			Splitting:   g.splitting,
			Integration: g.integration,
		})
	}

	// Sort by ppm descending (high-field right convention: decreasing ppm)
	sort.Slice(peaks, func(i, j int) bool {
		return peaks[i].X > peaks[j].X
	})

	return peaks
}

// classifyHEnv classifies the NMR environment of H's attached to atom i.
// Returns (ppm, splitting, label). Returns ("","","") to skip.
func classifyHEnv(m *molecule, i int) (float64, string, string) {
	a := &m.atoms[i]

	switch a.element {
	case "O":
		// Carboxylic acid OH
		for _, nb := range a.nbrs {
			if isCarboxylicAcidC(m, nb.idx) {
				return 11.5, "broad s", "COOH"
			}
		}
		// Phenol OH
		for _, nb := range a.nbrs {
			if m.atoms[nb.idx].aromatic {
				return 8.0, "broad s", "ArOH (phenol)"
			}
		}
		// Alcohol OH
		return 2.5, "broad s", "O–H (alcohol)"

	case "N":
		// Amide NH
		for _, nb := range a.nbrs {
			if isAmideC(m, nb.idx) {
				return 6.5, "broad s", "N–H (amide)"
			}
		}
		// Aromatic amine NH (aniline-type)
		for _, nb := range a.nbrs {
			if m.atoms[nb.idx].aromatic {
				return 3.5, "broad s", "N–H (ArNH)"
			}
		}
		return 1.5, "broad s", "N–H (amine)"

	case "C":
		return classifyCH(m, i)
	}

	return 0, "", ""
}

// classifyCH classifies the NMR environment of H on a carbon atom.
func classifyCH(m *molecule, i int) (float64, string, string) {
	a := &m.atoms[i]

	// Aromatic H
	if a.aromatic {
		return 7.3, "m", "Ar–H"
	}

	// Terminal alkyne H (sp C-H: C≡C-H)
	for _, nb := range a.nbrs {
		if nb.order == 3 && m.atoms[nb.idx].element == "C" {
			return 2.5, "s", "≡C–H (terminal alkyne)"
		}
	}

	// Aldehyde H (C=O on this carbon)
	for _, nb := range a.nbrs {
		if nb.order == 2 && m.atoms[nb.idx].element == "O" {
			return 9.7, "s", "CHO (aldehyde)"
		}
	}

	// Alkene H (sp2 C=C)
	for _, nb := range a.nbrs {
		if nb.order == 2 && m.atoms[nb.idx].element == "C" {
			return 5.5, "m", "=CH– (alkene)"
		}
	}

	// Sp3 — find most deshielding substituent on adjacent atoms
	return classifySp3CH(m, i)
}

// classifySp3CH returns the ppm, splitting, and label for an sp3 C-H environment.
func classifySp3CH(m *molecule, i int) (float64, string, string) {
	a := &m.atoms[i]

	// Scan neighbors for the most deshielding group
	maxShift := 0.0
	maxLabel := ""
	hasHalogen := false
	hasO := false
	hasN := false
	hasCarbonyl := false
	hasAromaticNbr := false

	for _, nb := range a.nbrs {
		nbAtom := &m.atoms[nb.idx]
		switch nbAtom.element {
		case "F":
			if 4.5 > maxShift {
				maxShift = 4.5
				maxLabel = "α-CH to F"
			}
			hasHalogen = true
		case "Cl":
			if 3.5 > maxShift {
				maxShift = 3.5
				maxLabel = "α-CH to Cl"
			}
			hasHalogen = true
		case "Br":
			if 3.3 > maxShift {
				maxShift = 3.3
				maxLabel = "α-CH to Br"
			}
			hasHalogen = true
		case "I":
			if 3.1 > maxShift {
				maxShift = 3.1
				maxLabel = "α-CH to I"
			}
			hasHalogen = true
		case "O":
			if !hasO {
				hasO = true
				// Distinguish ester α-C from ether/alcohol α-C
				shift := 3.5
				label := "α-CH to O"
				if 3.5 > maxShift {
					maxShift = shift
					maxLabel = label
				}
			}
		case "N":
			if !hasN {
				hasN = true
				if 2.5 > maxShift {
					maxShift = 2.5
					maxLabel = "α-CH to N"
				}
			}
		case "C":
			if nbAtom.aromatic && !hasAromaticNbr {
				hasAromaticNbr = true
				// Benzylic position
				if 2.5 > maxShift {
					maxShift = 2.5
					maxLabel = "benzylic CH"
				}
			}
			// Check if neighbor C is a carbonyl
			for _, nb2 := range nbAtom.nbrs {
				if nb2.order == 2 && m.atoms[nb2.idx].element == "O" && !hasCarbonyl {
					hasCarbonyl = true
					if 2.3 > maxShift {
						maxShift = 2.3
						maxLabel = "α-CH to C=O"
					}
				}
			}
		}
	}

	// Suppress unused variable warnings
	_ = hasHalogen

	sp := nPlusOne(m, i)

	if maxLabel != "" {
		return maxShift, sp, maxLabel
	}

	// Plain alkyl — classify by H count on this carbon
	switch {
	case a.hCount >= 3:
		return 0.9, sp, "CH₃ (alkyl)"
	case a.hCount == 2:
		return 1.3, sp, "CH₂ (alkyl)"
	default:
		return 1.5, sp, "CH (alkyl)"
	}
}

// nPlusOne returns the n+1 splitting pattern label by counting adjacent C-H's.
func nPlusOne(m *molecule, i int) string {
	adjH := 0
	for _, nb := range m.atoms[i].nbrs {
		nbAtom := &m.atoms[nb.idx]
		// Only count adjacent sp3 C H's for the n+1 rule
		if nbAtom.element == "C" && !nbAtom.aromatic && nb.order == 1 {
			adjH += nbAtom.hCount
		}
	}
	switch adjH {
	case 0:
		return "s"
	case 1:
		return "d"
	case 2:
		return "t"
	case 3:
		return "q"
	default:
		return "m"
	}
}

// ── ¹³C NMR ───────────────────────────────────────────────────────────────────

// estimateCNMR returns estimated ¹³C NMR peaks.
// Chemically equivalent carbons produce a single peak.
func estimateCNMR(m *molecule) []Peak {
	type cGroup struct {
		ppm   float64
		label string
	}

	sigToGroup := map[string]*cGroup{}

	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" {
			continue
		}

		ppm, label := classifyCNMR(m, i)
		sig := atomSig(m, i, 3, map[int]bool{})
		key := sig + "|" + label

		if _, ok := sigToGroup[key]; !ok {
			sigToGroup[key] = &cGroup{ppm: ppm, label: label}
		}
	}

	peaks := make([]Peak, 0, len(sigToGroup))
	for _, g := range sigToGroup {
		peaks = append(peaks, Peak{
			X:     g.ppm,
			Y:     70,
			Label: g.label,
			Width: 2,
		})
	}

	sort.Slice(peaks, func(i, j int) bool {
		return peaks[i].X > peaks[j].X
	})

	return peaks
}

// classifyCNMR returns the estimated ¹³C chemical shift and label for carbon i.
func classifyCNMR(m *molecule, i int) (float64, string) {
	a := &m.atoms[i]

	// Aromatic
	if a.aromatic {
		return 128, "aromatic C"
	}

	// Nitrile carbon (C≡N)
	for _, nb := range a.nbrs {
		if nb.order == 3 && m.atoms[nb.idx].element == "N" {
			return 118, "nitrile C (C≡N)"
		}
	}

	// Alkyne carbon (C≡C)
	for _, nb := range a.nbrs {
		if nb.order == 3 && m.atoms[nb.idx].element == "C" {
			return 75, "alkyne C"
		}
	}

	// Carbonyl carbons — check specific type
	for _, nb := range a.nbrs {
		if nb.order != 2 || m.atoms[nb.idx].element != "O" {
			continue
		}
		if isCarboxylicAcidC(m, i) {
			return 178, "carboxylic acid C=O"
		}
		if isEsterC(m, i) {
			return 170, "ester C=O"
		}
		if isAmideC(m, i) {
			return 168, "amide C=O"
		}
		if a.hCount > 0 {
			return 200, "aldehyde C=O"
		}
		return 205, "ketone C=O"
	}

	// Alkene carbon (C=C)
	for _, nb := range a.nbrs {
		if nb.order == 2 && m.atoms[nb.idx].element == "C" {
			return 125, "alkene C (C=C)"
		}
	}

	// Sp3 carbons — look for electron-withdrawing neighbors
	for _, nb := range a.nbrs {
		nbAtom := &m.atoms[nb.idx]
		switch nbAtom.element {
		case "O":
			return 65, "C bonded to O"
		case "N":
			return 50, "C bonded to N"
		case "F":
			return 90, "C bonded to F"
		case "Cl":
			return 40, "C bonded to Cl"
		case "Br":
			return 35, "C bonded to Br"
		case "I":
			return 20, "C bonded to I"
		}
	}

	// Plain alkyl: 1°/2°/3°/4° based on H count
	switch a.hCount {
	case 3:
		return 15, "alkyl C (CH₃)"
	case 2:
		return 25, "alkyl C (CH₂)"
	case 1:
		return 35, "alkyl C (CH)"
	default:
		return 40, "alkyl C (quaternary)"
	}
}
