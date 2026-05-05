package estimate

import "sort"

// estimateIR returns estimated IR peaks for a molecule.
// Source: Brown & Foote Organic Chemistry, Table 12.3.
func estimateIR(m *molecule) []Peak {
	var peaks []Peak

	// ── C-H stretches ─────────────────────────────────────────────────────────

	if hasTerminalAlkyneH(m) {
		peaks = append(peaks, Peak{X: 3300, Y: 90, Width: 30, Label: "≡C–H stretch (terminal alkyne)"})
	}
	if hasSp2CH(m) {
		peaks = append(peaks, Peak{X: 3050, Y: 65, Width: 40, Label: "sp² C–H stretch (vinyl/aryl)"})
	}
	if hasAldehyde(m) {
		// Aldehyde shows two weak C-H bands in this region
		peaks = append(peaks, Peak{X: 2820, Y: 40, Width: 30, Label: "aldehyde C–H (2720–2820, two bands)"})
	}
	if hasSp3CH(m) {
		peaks = append(peaks, Peak{X: 2920, Y: 80, Width: 60, Label: "sp³ C–H stretch"})
	}

	// ── O-H / N-H stretches ───────────────────────────────────────────────────

	if hasCarboxylicAcid(m) {
		peaks = append(peaks, Peak{X: 2700, Y: 60, Width: 600, Label: "O–H stretch (carboxylic acid, very broad)"})
	}
	if hasAlcohol(m) {
		peaks = append(peaks, Peak{X: 3400, Y: 80, Width: 200, Label: "O–H stretch (alcohol, broad)"})
	}
	if hasAmide(m) && hasPrimaryAmideNH(m) {
		peaks = append(peaks, Peak{X: 3300, Y: 70, Width: 80, Label: "N–H stretch (primary amide, two bands)"})
	} else if hasPrimaryAmineNH(m) {
		peaks = append(peaks, Peak{X: 3380, Y: 60, Width: 80, Label: "N–H stretch (primary amine, two bands)"})
		peaks = append(peaks, Peak{X: 3280, Y: 55, Width: 80, Label: "N–H stretch (primary amine)"})
	} else if hasSecondaryAmineNH(m) {
		peaks = append(peaks, Peak{X: 3330, Y: 55, Width: 70, Label: "N–H stretch (secondary amine)"})
	}

	// ── Triple bond region ────────────────────────────────────────────────────

	if hasAlkyne(m) {
		peaks = append(peaks, Peak{X: 2150, Y: 50, Width: 30, Label: "C≡C stretch"})
	}
	if hasNitrile(m) {
		peaks = append(peaks, Peak{X: 2240, Y: 85, Width: 25, Label: "C≡N stretch"})
	}

	// ── C=O stretches (most diagnostic) ──────────────────────────────────────

	if hasAcylChloride(m) {
		peaks = append(peaks, Peak{X: 1800, Y: 95, Width: 25, Label: "C=O stretch (acyl chloride)"})
	}
	if hasAnhydride(m) {
		peaks = append(peaks, Peak{X: 1820, Y: 90, Width: 30, Label: "C=O stretch (anhydride, high)"})
		peaks = append(peaks, Peak{X: 1760, Y: 85, Width: 30, Label: "C=O stretch (anhydride, low)"})
	}
	if hasEster(m) {
		peaks = append(peaks, Peak{X: 1740, Y: 95, Width: 25, Label: "C=O stretch (ester)"})
	}
	if hasAldehyde(m) {
		peaks = append(peaks, Peak{X: 1725, Y: 95, Width: 25, Label: "C=O stretch (aldehyde)"})
	}
	if hasCarboxylicAcid(m) {
		peaks = append(peaks, Peak{X: 1710, Y: 95, Width: 30, Label: "C=O stretch (carboxylic acid)"})
	}
	if hasKetone(m) {
		peaks = append(peaks, Peak{X: 1715, Y: 95, Width: 25, Label: "C=O stretch (ketone)"})
	}
	if hasAmide(m) {
		peaks = append(peaks, Peak{X: 1660, Y: 90, Width: 35, Label: "C=O stretch (amide)"})
	}

	// ── C=C stretches ─────────────────────────────────────────────────────────

	if hasAlkene(m) {
		peaks = append(peaks, Peak{X: 1650, Y: 45, Width: 30, Label: "C=C stretch (alkene)"})
	}
	if hasAromaticRing(m) {
		peaks = append(peaks, Peak{X: 1600, Y: 55, Width: 25, Label: "C=C stretch (aromatic)"})
		peaks = append(peaks, Peak{X: 1500, Y: 50, Width: 25, Label: "C=C stretch (aromatic)"})
	}

	// ── N=O ───────────────────────────────────────────────────────────────────

	if hasNitro(m) {
		peaks = append(peaks, Peak{X: 1530, Y: 85, Width: 25, Label: "N=O asymmetric stretch (nitro)"})
		peaks = append(peaks, Peak{X: 1350, Y: 80, Width: 25, Label: "N=O symmetric stretch (nitro)"})
	}

	// ── C-O single bond ───────────────────────────────────────────────────────

	if hasEster(m) {
		peaks = append(peaks, Peak{X: 1200, Y: 70, Width: 80, Label: "C–O stretch (ester)"})
	}
	if hasEther(m) {
		peaks = append(peaks, Peak{X: 1100, Y: 65, Width: 80, Label: "C–O stretch (ether)"})
	}
	if hasAlcohol(m) {
		peaks = append(peaks, Peak{X: 1050, Y: 60, Width: 80, Label: "C–O stretch (alcohol)"})
	}

	// ── Aromatic out-of-plane bend ────────────────────────────────────────────

	if hasAromaticRing(m) {
		peaks = append(peaks, Peak{X: 750, Y: 50, Width: 50, Label: "aromatic C–H out-of-plane bend"})
	}

	// Sort by wavenumber descending (standard IR presentation)
	sort.Slice(peaks, func(i, j int) bool {
		return peaks[i].X > peaks[j].X
	})

	return peaks
}

// ── Functional group detectors ────────────────────────────────────────────────

func hasSp3CH(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" || a.aromatic || a.hCount == 0 {
			continue
		}
		sp3 := true
		for _, nb := range a.nbrs {
			if nb.order > 1 {
				sp3 = false
				break
			}
		}
		if sp3 {
			return true
		}
	}
	return false
}

func hasSp2CH(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" || a.hCount == 0 {
			continue
		}
		if a.aromatic {
			return true
		}
		for _, nb := range a.nbrs {
			if nb.order == 2 && m.atoms[nb.idx].element == "C" {
				return true
			}
		}
	}
	return false
}

func hasTerminalAlkyneH(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" || a.hCount == 0 {
			continue
		}
		for _, nb := range a.nbrs {
			if nb.order == 3 && m.atoms[nb.idx].element == "C" {
				return true
			}
		}
	}
	return false
}

func hasAlkyne(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" {
			continue
		}
		for _, nb := range a.nbrs {
			if nb.order == 3 && m.atoms[nb.idx].element == "C" {
				return true
			}
		}
	}
	return false
}

func hasAlkene(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" || a.aromatic {
			continue
		}
		for _, nb := range a.nbrs {
			if nb.order == 2 && m.atoms[nb.idx].element == "C" && !m.atoms[nb.idx].aromatic {
				return true
			}
		}
	}
	return false
}

func hasAromaticRing(m *molecule) bool {
	for i := range m.atoms {
		if m.atoms[i].aromatic {
			return true
		}
	}
	return false
}

func hasNitrile(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" {
			continue
		}
		for _, nb := range a.nbrs {
			if nb.order == 3 && m.atoms[nb.idx].element == "N" {
				return true
			}
		}
	}
	return false
}

// hasCarbonylC returns true if atom i is a carbonyl carbon (C with =O neighbor).
func hasCarbonylC(m *molecule, i int) bool {
	a := &m.atoms[i]
	if a.element != "C" {
		return false
	}
	for _, nb := range a.nbrs {
		if nb.order == 2 && m.atoms[nb.idx].element == "O" {
			return true
		}
	}
	return false
}

func hasAldehyde(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" || a.hCount == 0 {
			continue
		}
		for _, nb := range a.nbrs {
			if nb.order == 2 && m.atoms[nb.idx].element == "O" {
				return true
			}
		}
	}
	return false
}

func hasKetone(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" || a.hCount > 0 {
			continue
		}
		cDoubleoO := false
		cNbrs := 0
		for _, nb := range a.nbrs {
			nbAtom := &m.atoms[nb.idx]
			if nb.order == 2 && nbAtom.element == "O" {
				cDoubleoO = true
			}
			if nbAtom.element == "C" {
				cNbrs++
			}
		}
		if cDoubleoO && cNbrs >= 2 {
			return true
		}
	}
	return false
}

// isCarboxylicAcidC returns true if atom i is the carbonyl C of a COOH group.
func isCarboxylicAcidC(m *molecule, i int) bool {
	a := &m.atoms[i]
	if a.element != "C" {
		return false
	}
	hasCarbonyl := false
	hasOH := false
	for _, nb := range a.nbrs {
		nbAtom := &m.atoms[nb.idx]
		if nb.order == 2 && nbAtom.element == "O" && nbAtom.hCount == 0 {
			hasCarbonyl = true
		}
		if nb.order == 1 && nbAtom.element == "O" && nbAtom.hCount > 0 {
			hasOH = true
		}
	}
	return hasCarbonyl && hasOH
}

func hasCarboxylicAcid(m *molecule) bool {
	for i := range m.atoms {
		if isCarboxylicAcidC(m, i) {
			return true
		}
	}
	return false
}

// isEsterC returns true if atom i is the carbonyl C of an ester (C(=O)OC).
func isEsterC(m *molecule, i int) bool {
	a := &m.atoms[i]
	if a.element != "C" {
		return false
	}
	hasCarbonyl := false
	hasOR := false
	for _, nb := range a.nbrs {
		nbAtom := &m.atoms[nb.idx]
		if nb.order == 2 && nbAtom.element == "O" {
			hasCarbonyl = true
		}
		if nb.order == 1 && nbAtom.element == "O" && nbAtom.hCount == 0 {
			// Check if this O is bonded to a C (not a second O)
			for _, nb2 := range nbAtom.nbrs {
				if nb2.idx != i && m.atoms[nb2.idx].element == "C" {
					hasOR = true
				}
			}
		}
	}
	return hasCarbonyl && hasOR
}

func hasEster(m *molecule) bool {
	for i := range m.atoms {
		if isEsterC(m, i) {
			return true
		}
	}
	return false
}

// isAmideC returns true if atom i is the carbonyl C of an amide (C(=O)N).
func isAmideC(m *molecule, i int) bool {
	a := &m.atoms[i]
	if a.element != "C" {
		return false
	}
	hasCarbonyl := false
	hasN := false
	for _, nb := range a.nbrs {
		nbAtom := &m.atoms[nb.idx]
		if nb.order == 2 && nbAtom.element == "O" {
			hasCarbonyl = true
		}
		if nbAtom.element == "N" {
			hasN = true
		}
	}
	return hasCarbonyl && hasN
}

func hasAmide(m *molecule) bool {
	for i := range m.atoms {
		if isAmideC(m, i) {
			return true
		}
	}
	return false
}

func hasPrimaryAmideNH(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "N" || a.hCount < 2 {
			continue
		}
		for _, nb := range a.nbrs {
			if isAmideC(m, nb.idx) {
				return true
			}
		}
	}
	return false
}

func hasAlcohol(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "O" || a.hCount == 0 {
			continue
		}
		// Must be bonded to C (not just floating)
		for _, nb := range a.nbrs {
			if m.atoms[nb.idx].element == "C" {
				// Make sure it's not a carboxylic acid OH (already handled separately)
				// We include COOH here too since we check separately above
				return true
			}
		}
	}
	return false
}

func hasPrimaryAmineNH(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "N" || a.hCount < 2 {
			continue
		}
		// Not an amide N
		isAmideN := false
		for _, nb := range a.nbrs {
			if isAmideC(m, nb.idx) {
				isAmideN = true
				break
			}
		}
		if !isAmideN {
			return true
		}
	}
	return false
}

func hasSecondaryAmineNH(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "N" || a.hCount != 1 {
			continue
		}
		isAmideN := false
		for _, nb := range a.nbrs {
			if isAmideC(m, nb.idx) {
				isAmideN = true
				break
			}
		}
		if !isAmideN {
			return true
		}
	}
	return false
}

func hasNitro(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "N" {
			continue
		}
		oNbrs := 0
		for _, nb := range a.nbrs {
			if m.atoms[nb.idx].element == "O" {
				oNbrs++
			}
		}
		if oNbrs >= 2 {
			return true
		}
	}
	return false
}

func hasAcylChloride(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" {
			continue
		}
		hasCarbonyl := false
		hasCl := false
		for _, nb := range a.nbrs {
			nbAtom := &m.atoms[nb.idx]
			if nb.order == 2 && nbAtom.element == "O" {
				hasCarbonyl = true
			}
			if nbAtom.element == "Cl" {
				hasCl = true
			}
		}
		if hasCarbonyl && hasCl {
			return true
		}
	}
	return false
}

func hasAnhydride(m *molecule) bool {
	// Pattern: C(=O)-O-C(=O)
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "O" || a.hCount > 0 {
			continue
		}
		// Count C neighbors that are carbonyl carbons
		carbonylNbrs := 0
		for _, nb := range a.nbrs {
			if hasCarbonylC(m, nb.idx) {
				carbonylNbrs++
			}
		}
		if carbonylNbrs >= 2 {
			return true
		}
	}
	return false
}

func hasEther(m *molecule) bool {
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "O" || a.hCount > 0 {
			continue
		}
		cNbrs := 0
		for _, nb := range a.nbrs {
			if m.atoms[nb.idx].element == "C" {
				cNbrs++
			}
		}
		if cNbrs >= 2 {
			return true
		}
	}
	return false
}
