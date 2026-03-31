package structure

import (
	"fmt"
	"strings"
	"unicode"

	"chemhelper/element"
)

// ── TM-Specific Data ──────────────────────────────────────────────────────────
// Only the data that cannot be derived from the element package lives here:
// common oxidation states and neutral d-electron counts.

type tmMeta struct {
	CommonOS []int // common oxidation states, ascending
	DNeutral int   // d-electron count for the neutral atom
}

var tmMetaData = map[string]tmMeta{
	"Sc": {[]int{3}, 1},
	"Ti": {[]int{2, 3, 4}, 2},
	"V":  {[]int{2, 3, 4, 5}, 3},
	"Cr": {[]int{2, 3, 6}, 5}, // anomalous 3d⁵4s¹
	"Mn": {[]int{2, 3, 4, 7}, 5},
	"Fe": {[]int{2, 3}, 6},
	"Co": {[]int{2, 3}, 7},
	"Ni": {[]int{2}, 8},
	"Cu": {[]int{1, 2}, 10}, // anomalous 3d¹⁰4s¹
	"Zn": {[]int{2}, 10},
	"Y":  {[]int{3}, 1},
	"Ag": {[]int{1}, 10},
	"Cd": {[]int{2}, 10},
	"Pt": {[]int{2, 4}, 9},
	"Au": {[]int{1, 3}, 10},
	"Hg": {[]int{1, 2}, 10},
}

// dElectronsForIon returns the d-electron count for a TM at a given oxidation state.
func dElectronsForIon(sym string, os int) int {
	meta, ok := tmMetaData[sym]
	if !ok {
		return 0
	}
	d := meta.DNeutral
	// Remove s-electrons first (up to 2), then d-electrons
	fromS := os
	if fromS > 2 {
		fromS = 2
	}
	d -= os - fromS
	if d < 0 {
		return 0
	}
	return d
}

// ── Group Formula Parser ──────────────────────────────────────────────────────

// formulaGroup is a component of a TM compound formula.
type formulaGroup struct {
	Symbol  string // element symbol, or inner formula for (...) groups
	Count   int
	IsGroup bool // true if this came from a (...) group
}

// parseGroups parses a formula that may contain parenthesised groups.
// It uses the element package to validate element symbols.
// Examples: "Fe(OH)2", "Mn(NO3)2", "Fe2(SO4)3", "FeCl3"
// The charge suffix must already be stripped before calling.
func parseGroups(s string) ([]formulaGroup, error) {
	var groups []formulaGroup
	i := 0
	for i < len(s) {
		switch {
		case s[i] == '(':
			depth := 1
			j := i + 1
			for j < len(s) && depth > 0 {
				if s[j] == '(' {
					depth++
				} else if s[j] == ')' {
					depth--
				}
				j++
			}
			if depth != 0 {
				return nil, fmt.Errorf("unmatched '(' in formula %q", s)
			}
			inner := s[i+1 : j-1]
			// Read count after ')'
			k := j
			for k < len(s) && s[k] >= '0' && s[k] <= '9' {
				k++
			}
			count := 1
			if k > j {
				count = 0
				for _, c := range s[j:k] {
					count = count*10 + int(c-'0')
				}
			}
			groups = append(groups, formulaGroup{Symbol: inner, Count: count, IsGroup: true})
			i = k

		case s[i] >= 'A' && s[i] <= 'Z':
			j := i + 1
			for j < len(s) && unicode.IsLower(rune(s[j])) {
				j++
			}
			sym := s[i:j]
			if _, ok := elemOf(sym); !ok {
				return nil, fmt.Errorf("unknown element %q in formula %q", sym, s)
			}
			k := j
			for k < len(s) && s[k] >= '0' && s[k] <= '9' {
				k++
			}
			count := 1
			if k > j {
				count = 0
				for _, c := range s[j:k] {
					count = count*10 + int(c-'0')
				}
			}
			groups = append(groups, formulaGroup{Symbol: sym, Count: count, IsGroup: false})
			i = k

		default:
			return nil, fmt.Errorf("unexpected character %q in formula", s[i])
		}
	}
	return groups, nil
}

// ── Ligand Charge Lookup ──────────────────────────────────────────────────────

// polyatomicLigandCharges maps normalised polyatomic ligand formulas to their
// ionic charges. Monatomic ligand charges are derived from element.MonatomicIonCharge().
var polyatomicLigandCharges = map[string]int{
	"oh":   -1,
	"cn":   -1,
	"no2":  -1,
	"no3":  -1,
	"scn":  -1,
	"hco3": -1,
	"clo":  -1,
	"clo3": -1,
	"clo4": -1,
	"co3":  -2,
	"so3":  -2,
	"so4":  -2,
	"c2o4": -2,
	"po4":  -3,
	// Neutral ligands
	"nh3": 0,
	"h2o": 0,
	"co":  0,
	"no":  0,
}

// ligandCharge returns the ionic charge for a ligand formula.
// For monatomic ligands it uses the element package; for polyatomic ones it
// uses the table above. Returns (charge, true) if known, (0, false) otherwise.
func ligandCharge(formula string) (int, bool) {
	key := strings.ToLower(formula)

	// Try polyatomic table first
	if ch, ok := polyatomicLigandCharges[key]; ok {
		return ch, true
	}

	// Single element: use element package
	if elem, ok := elemOf(formula); ok {
		return elem.MonatomicIonCharge(), true
	}

	// Try looking up the formula as a compound and infer from its standard ion
	// (falls back to 0 if unknown)
	return 0, false
}

// ── Oxidation State Inference ─────────────────────────────────────────────────

// inferOS determines the TM oxidation state from charge balance.
// overallCharge = tmCount × OS + Σ(ligand charges)
func inferOS(tmSym string, tmCount int, groups []formulaGroup, overallCharge int) (int, error) {
	totalLigandCharge := 0
	for _, g := range groups {
		if !g.IsGroup && g.Symbol == tmSym {
			continue
		}
		ch, ok := ligandCharge(g.Symbol)
		if !ok {
			return 0, fmt.Errorf("unknown ligand %q — cannot determine its charge for oxidation state inference", g.Symbol)
		}
		totalLigandCharge += ch * g.Count
	}

	numerator := overallCharge - totalLigandCharge
	if numerator%tmCount != 0 {
		return 0, fmt.Errorf(
			"cannot determine a whole-number oxidation state for %s: "+
				"(overall %+d − ligand %+d) / %d = %.1f",
			tmSym, overallCharge, totalLigandCharge, tmCount,
			float64(numerator)/float64(tmCount))
	}
	os := numerator / tmCount
	if os < -4 || os > 8 {
		return 0, fmt.Errorf("oxidation state %+d for %s is outside the supported range", os, tmSym)
	}
	return os, nil
}

// ── TM Entry Points ───────────────────────────────────────────────────────────

// generateTMGroup handles formulas with parentheses, e.g. "Fe(OH)2".
func generateTMGroup(input string) (*LewisStructure, error) {
	formula, overallCharge, err := extractCharge(input)
	if err != nil {
		return nil, err
	}
	groups, err := parseGroups(formula)
	if err != nil {
		return nil, err
	}
	tmSym, tmCount, err := findTMInGroups(groups)
	if err != nil {
		return nil, err
	}
	os, err := inferOS(tmSym, tmCount, groups, overallCharge)
	if err != nil {
		return nil, err
	}
	return buildTMStructure(tmSym, tmCount, groups, os, overallCharge)
}

// generateTMFlat handles flat TM formulas, e.g. "FeCl3", "Fe2O3".
func generateTMFlat(atoms []atomCount, overallCharge int) (*LewisStructure, error) {
	var groups []formulaGroup
	for _, ac := range atoms {
		groups = append(groups, formulaGroup{Symbol: ac.Symbol, Count: ac.Count, IsGroup: false})
	}
	tmSym, tmCount, err := findTMInGroups(groups)
	if err != nil {
		return nil, err
	}
	os, err := inferOS(tmSym, tmCount, groups, overallCharge)
	if err != nil {
		return nil, err
	}
	return buildTMStructure(tmSym, tmCount, groups, os, overallCharge)
}

// findTMInGroups locates the single metal in a group list.
// It accepts transition metals (groups 3–12) and main-group metals
// (alkali, alkaline earth, and post-transition metals).
func findTMInGroups(groups []formulaGroup) (sym string, count int, err error) {
	for _, g := range groups {
		if g.IsGroup {
			continue
		}
		if e, ok := elemOf(g.Symbol); ok && e.IsMetal() {
			if sym != "" && sym != g.Symbol {
				return "", 0, fmt.Errorf("multiple different metals in one formula are not supported")
			}
			sym = g.Symbol
			count += g.Count
		}
	}
	if sym == "" {
		return "", 0, fmt.Errorf("no metal found in formula")
	}
	return sym, count, nil
}

// ── Structure Builder ─────────────────────────────────────────────────────────

func buildTMStructure(tmSym string, tmCount int, groups []formulaGroup, os, overallCharge int) (*LewisStructure, error) {
	tmElem, _ := elemOf(tmSym)
	dElec := dElectronsForIon(tmSym, os)

	// Build display formula
	displayFormula := buildTMDisplayFormula(groups, overallCharge)

	var allAtoms []LewisAtom
	var allBonds []LewisBond

	// Place TM atoms
	for ti := 1; ti <= tmCount; ti++ {
		allAtoms = append(allAtoms, LewisAtom{
			ID:           fmt.Sprintf("%s%d", tmSym, ti),
			Element:      tmSym,
			LonePairs:    0,
			FormalCharge: os,
		})
	}

	// Distribute ligands across TM atoms
	totalLigandSlots := 0
	for _, g := range groups {
		if g.IsGroup || g.Symbol != tmSym {
			totalLigandSlots += g.Count
		}
	}
	ligandsPerTM := totalLigandSlots
	if tmCount > 1 && totalLigandSlots%tmCount == 0 {
		ligandsPerTM = totalLigandSlots / tmCount
	}

	// globalIdx tracks per-element atom numbering for unique IDs
	globalIdx := map[string]int{}
	// Reserve TM indices
	for ti := 1; ti <= tmCount; ti++ {
		globalIdx[tmSym] = ti
	}

	tmIdx := 1
	ligandAttached := 0

	for _, g := range groups {
		if !g.IsGroup && g.Symbol == tmSym {
			continue
		}
		for rep := 0; rep < g.Count; rep++ {
			if tmCount > 1 && ligandAttached > 0 && ligandAttached%ligandsPerTM == 0 && tmIdx < tmCount {
				tmIdx++
			}
			ligandAttached++
			currentTMID := fmt.Sprintf("%s%d", tmSym, tmIdx)

			ligAtoms, ligBonds, coordID, err := expandLigand(g, globalIdx)
			if err != nil {
				return nil, fmt.Errorf("cannot expand ligand %q: %w", g.Symbol, err)
			}
			allAtoms = append(allAtoms, ligAtoms...)
			allBonds = append(allBonds, ligBonds...)
			allBonds = append(allBonds, LewisBond{From: currentTMID, To: coordID, Order: 1})
		}
	}

	// Coordination number = bonds on TM1
	coordNum := 0
	tmID1 := fmt.Sprintf("%s1", tmSym)
	for _, b := range allBonds {
		if b.From == tmID1 || b.To == tmID1 {
			coordNum++
		}
	}
	geom := coordinationGeometry(coordNum)

	totalLigCharge, breakdown := ligandChargeBreakdown(groups, tmSym)
	notes := buildTMNotes(tmElem, tmSym, os, dElec, totalLigCharge, overallCharge, tmCount, breakdown)
	steps := buildTMSteps(tmElem, tmSym, tmCount, groups, os, overallCharge, dElec, coordNum, totalLigCharge, breakdown)

	return &LewisStructure{
		Name:                  displayFormula,
		Formula:               displayFormula,
		Charge:                overallCharge,
		TotalValenceElectrons: 0, // not meaningful for ionic model
		Geometry:              geom,
		Atoms:                 allAtoms,
		Bonds:                 allBonds,
		Steps:                 steps,
		Notes:                 notes,
	}, nil
}

// ── Ligand Expansion ──────────────────────────────────────────────────────────

// expandLigand produces the atoms, internal bonds, and coordinating atom ID for
// one ligand instance. For polyatomic ligands it calls generate() recursively.
// globalIdx is updated in-place to ensure unique atom IDs across the whole structure.
func expandLigand(g formulaGroup, globalIdx map[string]int) (atoms []LewisAtom, bonds []LewisBond, coordID string, err error) {
	if !g.IsGroup {
		// Monatomic ligand: derive LP and FC from element package
		elem, ok := elemOf(g.Symbol)
		if !ok {
			return nil, nil, "", fmt.Errorf("unknown element %q", g.Symbol)
		}
		globalIdx[g.Symbol]++
		atomID := fmt.Sprintf("%s%d", g.Symbol, globalIdx[g.Symbol])
		ionicCharge := elem.MonatomicIonCharge()
		// LP for a monatomic ion bonded once to TM:
		// After forming the M-L bond, the ligand atom has (8 - ionicCharge*(-2) - 2) / 2 LP
		// More simply: group 17 → 3 LP (halide); group 16 → 3 LP (one bond to TM); group 15 → 2 LP
		lonePairs := (elem.Valence() - 1 - ionicCharge) / 2
		// ionicCharge is negative for anions, so: Cl (7 - 1 - (-1)) / 2 = 7/2 → 3 (integer div)
		// O: (6 - 1 - (-2)) / 2 = 7/2 → 3 ✓
		// N: (5 - 1 - (-3)) / 2 = 7/2 → 3 ✓
		if lonePairs < 0 {
			lonePairs = 0
		}
		atoms = append(atoms, LewisAtom{
			ID:           atomID,
			Element:      g.Symbol,
			LonePairs:    lonePairs,
			FormalCharge: ionicCharge,
		})
		return atoms, nil, atomID, nil
	}

	// Polyatomic ligand: look up the charge and call generate() recursively
	ch, ok := ligandCharge(g.Symbol)
	if !ok {
		return nil, nil, "", fmt.Errorf("unknown polyatomic ligand %q", g.Symbol)
	}

	// Build the formula with charge suffix for the generator
	chargeStr := ""
	switch ch {
	case 0:
	case 1:
		chargeStr = "+"
	case -1:
		chargeStr = "-"
	default:
		if ch > 0 {
			chargeStr = fmt.Sprintf("%d+", ch)
		} else {
			chargeStr = fmt.Sprintf("%d-", -ch)
		}
	}
	ls, err := generate(g.Symbol + chargeStr)
	if err != nil {
		return nil, nil, "", fmt.Errorf("cannot generate Lewis structure for ligand %s: %w", g.Symbol, err)
	}

	// Re-number the atoms from the ligand structure to avoid ID collisions
	idRemap := map[string]string{}
	for _, a := range ls.Atoms {
		sym := symbolOf(a.ID)
		globalIdx[sym]++
		newID := fmt.Sprintf("%s%d", sym, globalIdx[sym])
		idRemap[a.ID] = newID
		atoms = append(atoms, LewisAtom{
			ID:           newID,
			Element:      a.Element,
			LonePairs:    a.LonePairs,
			FormalCharge: a.FormalCharge,
		})
	}
	for _, b := range ls.Bonds {
		bonds = append(bonds, LewisBond{
			From:  idRemap[b.From],
			To:    idRemap[b.To],
			Order: b.Order,
		})
	}

	// The coordinating atom is the least electronegative atom in the ligand
	// (i.e. the one that donates to the metal). For most ligands this is the
	// first heavy atom. We pick the atom with the lowest EN that isn't H.
	coordID = idRemap[ls.Atoms[0].ID]
	bestEN := 999.0
	for _, a := range ls.Atoms {
		if a.Element == "H" {
			continue
		}
		if en := elemEN(a.Element); en < bestEN {
			bestEN = en
			coordID = idRemap[a.ID]
		}
	}

	return atoms, bonds, coordID, nil
}

// ── Note and Step Builders ────────────────────────────────────────────────────

func buildTMNotes(tmElem *element.Element, tmSym string, os, dElec, totalLigCharge, overallCharge, tmCount int, breakdown string) string {
	var parts []string
	if tmElem.IsTransitionMetal() {
		parts = append(parts, fmt.Sprintf(
			"%s is in oxidation state %s (d%d electron configuration).",
			tmElem.Name, ionicChargeLabel(os), dElec))
	} else {
		parts = append(parts, fmt.Sprintf(
			"%s is in oxidation state %s.",
			tmElem.Name, ionicChargeLabel(os)))
	}
	parts = append(parts, fmt.Sprintf(
		"Lewis structures for ionic compounds use the ionic model: "+
			"the metal is shown as %s%s with formal charge %+d and no lone pairs displayed.",
		tmSym, ionicChargeLabel(os), os))
	if breakdown != "" {
		parts = append(parts, fmt.Sprintf(
			"Charge balance: %d × (%+d) + (%s) = %+d. ✓",
			tmCount, os, breakdown, overallCharge))
	}
	return strings.Join(parts, " ")
}

func buildTMSteps(tmElem *element.Element, tmSym string, tmCount int, groups []formulaGroup, os, overallCharge, dElec, coordNum, totalLigCharge int, breakdown string) []string {
	var steps []string
	if tmElem.IsTransitionMetal() {
		steps = append(steps, fmt.Sprintf("Identify the transition metal: %s (%s)", tmElem.Name, tmSym))
	} else {
		steps = append(steps, fmt.Sprintf("Identify the metal: %s (%s)", tmElem.Name, tmSym))
	}
	steps = append(steps, fmt.Sprintf("Determine %s oxidation state using charge balance:", tmSym))
	steps = append(steps, fmt.Sprintf("  Ligand charges: %s → total = %+d", breakdown, totalLigCharge))
	if tmCount == 1 {
		steps = append(steps, fmt.Sprintf(
			"  OS(%s) = overall − ligands = %+d − (%+d) = %+d",
			tmSym, overallCharge, totalLigCharge, os))
	} else {
		steps = append(steps, fmt.Sprintf(
			"  OS(%s) = (overall − ligands) / %d = (%+d − (%+d)) / %d = %+d",
			tmSym, tmCount, overallCharge, totalLigCharge, tmCount, os))
	}
	if tmElem.IsTransitionMetal() {
		steps = append(steps, fmt.Sprintf(
			"%s%s has a d%d electron configuration", tmSym, ionicChargeLabel(os), dElec))
	}
	steps = append(steps, fmt.Sprintf(
		"Apply the ionic model: draw %s as %s%s (formal charge %+d, no lone pairs)",
		tmSym, tmSym, ionicChargeLabel(os), os))

	for _, g := range groups {
		if !g.IsGroup && g.Symbol == tmSym {
			continue
		}
		ch, _ := ligandCharge(g.Symbol)
		label := ionicChargeLabel(ch)
		if g.Count == 1 {
			steps = append(steps, fmt.Sprintf(
				"Draw %s%s ligand bonded to %s via its Lewis structure", g.Symbol, label, tmSym))
		} else {
			steps = append(steps, fmt.Sprintf(
				"Draw %d × %s%s ligands, each bonded to %s", g.Count, g.Symbol, label, tmSym))
		}
	}

	steps = append(steps, fmt.Sprintf(
		"Verify charge balance: %d × (%+d OS) + (%+d ligand charges) = %+d overall ✓",
		tmCount, os, totalLigCharge, overallCharge))
	steps = append(steps, fmt.Sprintf(
		"Coordination number of %s: %d → %s geometry", tmSym, coordNum, coordinationGeometry(coordNum)))
	return steps
}

func ligandChargeBreakdown(groups []formulaGroup, tmSym string) (total int, breakdown string) {
	var parts []string
	for _, g := range groups {
		if !g.IsGroup && g.Symbol == tmSym {
			continue
		}
		ch, _ := ligandCharge(g.Symbol)
		total += ch * g.Count
		if g.Count == 1 {
			parts = append(parts, fmt.Sprintf("%s(%+d)", g.Symbol, ch))
		} else {
			parts = append(parts, fmt.Sprintf("%d×%s(%+d)", g.Count, g.Symbol, ch))
		}
	}
	return total, strings.Join(parts, " + ")
}

func buildTMDisplayFormula(groups []formulaGroup, charge int) string {
	var b strings.Builder
	for _, g := range groups {
		if g.IsGroup {
			b.WriteString("(")
			b.WriteString(g.Symbol)
			b.WriteString(")")
			if g.Count > 1 {
				b.WriteString(toSubscript(g.Count))
			}
		} else {
			b.WriteString(g.Symbol)
			if g.Count > 1 {
				b.WriteString(toSubscript(g.Count))
			}
		}
	}
	switch {
	case charge == 1:
		b.WriteRune('⁺')
	case charge == -1:
		b.WriteRune('⁻')
	case charge > 1:
		b.WriteString(toSuperscript(charge) + "⁺")
	case charge < -1:
		b.WriteString(toSuperscript(-charge) + "⁻")
	}
	return b.String()
}

func coordinationGeometry(n int) string {
	switch n {
	case 2:
		return "linear"
	case 3:
		return "trigonal_planar"
	case 4:
		return "tetrahedral"
	case 5:
		return "trigonal_bipyramidal"
	case 6:
		return "octahedral"
	default:
		return "unknown"
	}
}

func ionicChargeLabel(charge int) string {
	switch charge {
	case 0:
		return ""
	case 1:
		return "⁺"
	case -1:
		return "⁻"
	}
	if charge > 1 {
		return toSuperscript(charge) + "⁺"
	}
	return toSuperscript(-charge) + "⁻"
}
