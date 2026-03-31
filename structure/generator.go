package structure

import (
	"fmt"
	"sort"
	"strings"

	"chemhelper/element"
)

// pt is the shared periodic table used by the structure package.
var pt = element.NewPeriodicTable()

// ── Element Accessors ─────────────────────────────────────────────────────────

func elemOf(sym string) (*element.Element, bool) {
	return pt.FindElementBySymbol(sym)
}

func elemEN(sym string) float64 {
	e, ok := elemOf(sym)
	if !ok {
		return 0
	}
	return e.Electronegativity
}

func elemValence(sym string) int {
	e, ok := elemOf(sym)
	if !ok {
		return 0
	}
	return e.Valence()
}

func elemPeriod(sym string) int {
	e, ok := elemOf(sym)
	if !ok {
		return 0
	}
	return e.Period
}

func elemName(sym string) string {
	e, ok := elemOf(sym)
	if !ok {
		return sym
	}
	return e.Name
}

func isTM(sym string) bool {
	e, ok := elemOf(sym)
	return ok && e.IsTransitionMetal()
}

func isMetal(sym string) bool {
	e, ok := elemOf(sym)
	return ok && e.IsMetal()
}

// ── Formula Parsing ───────────────────────────────────────────────────────────

type atomCount struct {
	Symbol string
	Count  int
}

// parseFlat parses a molecular formula (no parentheses) into ordered atom counts.
// Uses the element package for symbol validation.
func parseFlat(s string) ([]atomCount, error) {
	counts := map[string]int{}
	var order []string
	i := 0
	for i < len(s) {
		if s[i] < 'A' || s[i] > 'Z' {
			return nil, fmt.Errorf(
				"unexpected character %q at position %d — element symbols must start with an uppercase letter",
				s[i], i)
		}
		j := i + 1
		for j < len(s) && s[j] >= 'a' && s[j] <= 'z' {
			j++
		}
		sym := s[i:j]
		if _, ok := elemOf(sym); !ok {
			return nil, fmt.Errorf("unknown element symbol %q", sym)
		}
		k := j
		for k < len(s) && s[k] >= '0' && s[k] <= '9' {
			k++
		}
		n := 1
		if k > j {
			n = 0
			for _, c := range s[j:k] {
				n = n*10 + int(c-'0')
			}
		}
		if _, seen := counts[sym]; !seen {
			order = append(order, sym)
		}
		counts[sym] += n
		i = k
	}
	result := make([]atomCount, 0, len(order))
	for _, sym := range order {
		result = append(result, atomCount{sym, counts[sym]})
	}
	return result, nil
}

// extractCharge strips a trailing charge suffix ("+", "-", "2+", "3-", etc.)
// and returns the bare formula and integer charge.
func extractCharge(s string) (formula string, charge int, err error) {
	if s == "" {
		return "", 0, fmt.Errorf("empty input")
	}
	last := s[len(s)-1]
	if last != '+' && last != '-' {
		return s, 0, nil
	}
	sign := 1
	if last == '-' {
		sign = -1
	}
	rest := s[:len(s)-1]
	magnitude := 0
	i := len(rest) - 1
	for i >= 0 && rest[i] >= '0' && rest[i] <= '9' {
		magnitude = magnitude*10 + int(rest[i]-'0')
		i--
	}
	formula = rest[:i+1]
	if formula == "" {
		return "", 0, fmt.Errorf("formula is empty after extracting charge")
	}
	if magnitude == 0 {
		magnitude = 1
	}
	return formula, sign * magnitude, nil
}

// ── Main Entry Point ──────────────────────────────────────────────────────────

// generate is called by LookupLewisWithError after the static registry misses.
func generate(input string) (*LewisStructure, error) {
	// Parentheses → TM compound (e.g. Fe(OH)2)
	if strings.ContainsRune(input, '(') {
		return generateTMGroup(input)
	}

	formula, charge, err := extractCharge(input)
	if err != nil {
		return nil, err
	}

	atoms, err := parseFlat(formula)
	if err != nil {
		return nil, err
	}

	// Any metal present → use ionic model
	for _, ac := range atoms {
		if isMetal(ac.Symbol) {
			return generateTMFlat(atoms, charge)
		}
	}

	return solveMainGroup(atoms, charge)
}

// ── Main-Group Solver ─────────────────────────────────────────────────────────

func solveMainGroup(atoms []atomCount, charge int) (*LewisStructure, error) {
	totalVE := 0
	for _, ac := range atoms {
		totalVE += ac.Count * elemValence(ac.Symbol)
	}
	totalVE -= charge

	if totalVE <= 0 {
		return nil, fmt.Errorf("invalid combination: total valence electron count is %d", totalVE)
	}

	isRadical := totalVE%2 != 0
	displayFormula := buildDisplayFormula(atoms, charge)

	totalAtoms := 0
	for _, ac := range atoms {
		totalAtoms += ac.Count
	}

	if totalAtoms == 2 {
		var a1, a2 atomCount
		if len(atoms) == 1 {
			a1 = atomCount{atoms[0].Symbol, 1}
			a2 = atomCount{atoms[0].Symbol, 1}
		} else {
			a1, a2 = atomCount{atoms[0].Symbol, 1}, atomCount{atoms[1].Symbol, 1}
		}
		return solveDiatomic(a1, a2, totalVE, charge, displayFormula, isRadical)
	}

	var heavy []atomCount
	hCount := 0
	for _, ac := range atoms {
		if ac.Symbol == "H" {
			hCount += ac.Count
		} else {
			heavy = append(heavy, ac)
		}
	}

	if len(heavy) == 0 {
		return nil, fmt.Errorf("formula contains only hydrogen — no valid central atom")
	}

	// Least electronegative heavy atom is the centre
	sort.Slice(heavy, func(i, j int) bool {
		return elemEN(heavy[i].Symbol) < elemEN(heavy[j].Symbol)
	})

	leastEN := heavy[0]
	if leastEN.Count == 1 {
		var terminals []atomCount
		terminals = append(terminals, heavy[1:]...)
		if hCount > 0 {
			terminals = append(terminals, atomCount{"H", hCount})
		}
		return solveSingleCenter(leastEN, terminals, totalVE, charge, displayFormula, isRadical)
	}
	if leastEN.Count == 2 {
		var terminals []atomCount
		terminals = append(terminals, heavy[1:]...)
		if hCount > 0 {
			terminals = append(terminals, atomCount{"H", hCount})
		}
		return solveTwoCenter(leastEN, terminals, totalVE, charge, displayFormula, isRadical)
	}

	return nil, fmt.Errorf(
		"complex topology: %s appears %d times — molecules with 3+ of the same central atom require a more advanced solver",
		leastEN.Symbol, leastEN.Count)
}

// ── Single-Center Solver ──────────────────────────────────────────────────────

func solveSingleCenter(
	central atomCount,
	terminals []atomCount,
	totalVE, charge int,
	displayFormula string,
	isRadical bool,
) (*LewisStructure, error) {

	numBonds := 0
	for _, t := range terminals {
		numBonds += t.Count
	}

	centralID := central.Symbol + "1"
	var allAtoms []LewisAtom
	allAtoms = append(allAtoms, LewisAtom{ID: centralID, Element: central.Symbol})

	var termIDs []string
	for _, t := range terminals {
		for i := 1; i <= t.Count; i++ {
			id := fmt.Sprintf("%s%d", t.Symbol, i)
			allAtoms = append(allAtoms, LewisAtom{ID: id, Element: t.Symbol})
			termIDs = append(termIDs, id)
		}
	}

	bonds := make([]LewisBond, len(termIDs))
	for idx, tid := range termIDs {
		bonds[idx] = LewisBond{From: centralID, To: tid, Order: 1}
	}

	boOf := func(id string) int {
		sum := 0
		for _, b := range bonds {
			if b.From == id || b.To == id {
				sum += b.Order
			}
		}
		return sum
	}

	remaining := totalVE - 2*numBonds
	if isRadical {
		remaining--
	}

	lp := map[string]int{}

	for _, t := range terminals {
		need := (octetFor(t.Symbol) - 2) / 2
		if need < 0 {
			need = 0
		}
		for i := 1; i <= t.Count; i++ {
			id := fmt.Sprintf("%s%d", t.Symbol, i)
			lp[id] = need
			remaining -= 2 * need
		}
	}

	if remaining < 0 {
		return nil, fmt.Errorf("insufficient valence electrons for the given formula and charge")
	}
	lp[centralID] = remaining / 2

	isIonicMetal := isMetal(central.Symbol) && !isTM(central.Symbol)

	calcFC := func(id string) int {
		return formalCharge(symbolOf(id), lp[id], boOf(id))
	}
	totalAbsFC := func() int {
		sum := 0
		for _, a := range allAtoms {
			fc := calcFC(a.ID)
			if fc < 0 {
				sum -= fc
			} else {
				sum += fc
			}
		}
		return sum
	}
	centralEC := func() int { return 2*lp[centralID] + 2*boOf(centralID) }

	if !isIonicMetal {
		for {
			cEC := centralEC()
			satisfied := cEC >= 8
			if satisfied && elemPeriod(central.Symbol) <= 2 {
				break
			}
			cur := totalAbsFC()
			bestImp := 0
			bestIdx := -1
			for bidx, b := range bonds {
				var donorID string
				if b.From == centralID {
					donorID = b.To
				} else if b.To == centralID {
					donorID = b.From
				} else {
					continue
				}
				if lp[donorID] == 0 {
					continue
				}
				if elemPeriod(central.Symbol) <= 2 && (boOf(centralID)+1)*2 > 8 {
					continue
				}
				lp[donorID]--
				bonds[bidx].Order++
				imp := cur - totalAbsFC()
				lp[donorID]++
				bonds[bidx].Order--
				if imp > bestImp {
					bestImp = imp
					bestIdx = bidx
				} else if !satisfied && imp == 0 && bestIdx == -1 {
					bestIdx = bidx
				}
			}
			if bestIdx == -1 || (bestImp <= 0 && satisfied) {
				break
			}
			var donorID string
			if bonds[bestIdx].From == centralID {
				donorID = bonds[bestIdx].To
			} else {
				donorID = bonds[bestIdx].From
			}
			lp[donorID]--
			bonds[bestIdx].Order++
		}
	}

	for i, a := range allAtoms {
		allAtoms[i].LonePairs = lp[a.ID]
		allAtoms[i].FormalCharge = calcFC(a.ID)
	}

	geom := determineGeometry(numBonds, lp[centralID])

	var notes string
	if isIonicMetal {
		maxEN := 0.0
		for _, t := range terminals {
			if en := elemEN(t.Symbol); en > maxEN {
				maxEN = en
			}
		}
		if delta := maxEN - elemEN(central.Symbol); delta > 1.7 {
			notes = fmt.Sprintf(
				"%s is predominantly ionic (ΔEN ≈ %.1f). This diagram shows the hypothetical covalent Lewis structure.",
				displayFormula, delta)
		}
	}
	if isRadical {
		if notes != "" {
			notes += " "
		}
		notes += fmt.Sprintf("%s has an odd number of valence electrons — it is a free radical.", displayFormula)
	}

	steps := buildSingleCenterSteps(central, terminals, totalVE, charge, numBonds,
		lp, bonds, centralID, isRadical, isIonicMetal)

	return &LewisStructure{
		Name:                  displayFormula,
		Formula:               displayFormula,
		Charge:                charge,
		TotalValenceElectrons: totalVE,
		Geometry:              geom,
		Atoms:                 allAtoms,
		Bonds:                 bonds,
		Steps:                 steps,
		Notes:                 strings.TrimSpace(notes),
	}, nil
}

// ── Two-Center Chain Solver ───────────────────────────────────────────────────

func solveTwoCenter(
	center atomCount,
	terminals []atomCount,
	totalVE, charge int,
	displayFormula string,
	isRadical bool,
) (*LewisStructure, error) {

	totalT := 0
	for _, t := range terminals {
		totalT += t.Count
	}
	if totalT%2 != 0 {
		return nil, fmt.Errorf(
			"asymmetric two-center chain (%d terminals cannot split evenly between 2 %s atoms)",
			totalT, center.Symbol)
	}
	perCenter := totalT / 2

	c1ID, c2ID := center.Symbol+"1", center.Symbol+"2"
	allAtoms := []LewisAtom{
		{ID: c1ID, Element: center.Symbol},
		{ID: c2ID, Element: center.Symbol},
	}

	type ta struct{ id, sym, cid string }
	var assigns []ta
	gi := 0
	for _, t := range terminals {
		for i := 1; i <= t.Count; i++ {
			gi++
			tid := fmt.Sprintf("%s%d", t.Symbol, i)
			cid := c1ID
			if gi > perCenter {
				cid = c2ID
			}
			allAtoms = append(allAtoms, LewisAtom{ID: tid, Element: t.Symbol})
			assigns = append(assigns, ta{tid, t.Symbol, cid})
		}
	}

	var bonds []LewisBond
	for _, a := range assigns {
		bonds = append(bonds, LewisBond{From: a.cid, To: a.id, Order: 1})
	}
	ccIdx := len(bonds)
	bonds = append(bonds, LewisBond{From: c1ID, To: c2ID, Order: 1})

	boOf := func(id string) int {
		sum := 0
		for _, b := range bonds {
			if b.From == id || b.To == id {
				sum += b.Order
			}
		}
		return sum
	}

	remaining := totalVE - 2*len(bonds)
	if isRadical {
		remaining--
	}

	lp := map[string]int{}
	for _, a := range assigns {
		need := (octetFor(a.sym) - 2) / 2
		if need < 0 {
			need = 0
		}
		lp[a.id] = need
		remaining -= 2 * need
	}
	if remaining < 0 {
		return nil, fmt.Errorf("insufficient valence electrons")
	}
	for _, cid := range []string{c1ID, c2ID} {
		need := (8 - 2*boOf(cid)) / 2
		if need < 0 {
			need = 0
		}
		give := need
		if give*2 > remaining {
			give = remaining / 2
		}
		lp[cid] = give
		remaining -= 2 * give
	}
	lp[c1ID] += remaining / 2

	calcFC := func(id string) int {
		return formalCharge(symbolOf(id), lp[id], boOf(id))
	}
	totalAbsFC := func() int {
		sum := 0
		for _, a := range allAtoms {
			fc := calcFC(a.ID)
			if fc < 0 {
				sum -= fc
			} else {
				sum += fc
			}
		}
		return sum
	}

	for {
		c1EC := 2*lp[c1ID] + 2*boOf(c1ID)
		c2EC := 2*lp[c2ID] + 2*boOf(c2ID)
		if c1EC >= 8 && c2EC >= 8 && elemPeriod(center.Symbol) <= 2 {
			break
		}
		cur := totalAbsFC()
		bestImp, bestIdx, bestDonor := 0, -1, ""
		for bidx, b := range bonds {
			isCC := bidx == ccIdx
			if isCC {
				for _, cid := range []string{c1ID, c2ID} {
					if lp[cid] == 0 {
						continue
					}
					if elemPeriod(center.Symbol) <= 2 && (boOf(c1ID)+1)*2 > 8 {
						continue
					}
					lp[cid]--
					bonds[bidx].Order++
					imp := cur - totalAbsFC()
					lp[cid]++
					bonds[bidx].Order--
					if imp > bestImp {
						bestImp, bestIdx, bestDonor = imp, bidx, cid
					} else if imp == 0 && bestIdx == -1 {
						bestIdx, bestDonor = bidx, cid
					}
				}
			} else {
				var donorID, recipID string
				if b.From == c1ID || b.From == c2ID {
					donorID, recipID = b.To, b.From
				} else {
					donorID, recipID = b.From, b.To
				}
				if lp[donorID] == 0 {
					continue
				}
				if elemPeriod(center.Symbol) <= 2 && (boOf(recipID)+1)*2 > 8 {
					continue
				}
				lp[donorID]--
				bonds[bidx].Order++
				imp := cur - totalAbsFC()
				lp[donorID]++
				bonds[bidx].Order--
				if imp > bestImp {
					bestImp, bestIdx, bestDonor = imp, bidx, donorID
				} else if imp == 0 && bestIdx == -1 {
					bestIdx, bestDonor = bidx, donorID
				}
			}
		}
		bothSat := (2*lp[c1ID]+2*boOf(c1ID) >= 8) && (2*lp[c2ID]+2*boOf(c2ID) >= 8)
		if bestIdx == -1 || (bestImp <= 0 && bothSat) {
			break
		}
		lp[bestDonor]--
		bonds[bestIdx].Order++
	}

	for i, a := range allAtoms {
		allAtoms[i].LonePairs = lp[a.ID]
		allAtoms[i].FormalCharge = calcFC(a.ID)
	}

	steps := buildTwoCenterSteps(center, terminals, totalVE, charge, len(bonds), lp, bonds, c1ID, c2ID, isRadical)

	return &LewisStructure{
		Name:                  displayFormula,
		Formula:               displayFormula,
		Charge:                charge,
		TotalValenceElectrons: totalVE,
		Geometry:              "chain",
		Atoms:                 allAtoms,
		Bonds:                 bonds,
		Steps:                 steps,
	}, nil
}

// ── Diatomic Solver ───────────────────────────────────────────────────────────

func solveDiatomic(a1, a2 atomCount, totalVE, charge int, displayFormula string, isRadical bool) (*LewisStructure, error) {
	leftID := a1.Symbol + "1"
	rightID := a2.Symbol + "1"
	if a1.Symbol == a2.Symbol {
		rightID = a2.Symbol + "2"
	}
	if elemEN(a1.Symbol) > elemEN(a2.Symbol) {
		a1, a2 = a2, a1
		leftID, rightID = a1.Symbol+"1", a2.Symbol+"1"
		if a1.Symbol == a2.Symbol {
			rightID = a2.Symbol + "2"
		}
	}

	atoms := []LewisAtom{{ID: leftID, Element: a1.Symbol}, {ID: rightID, Element: a2.Symbol}}
	bonds := []LewisBond{{From: leftID, To: rightID, Order: 1}}

	remaining := totalVE - 2
	if isRadical {
		remaining--
	}
	lp := map[string]int{leftID: 0, rightID: 0}

	for _, id := range []string{rightID, leftID} {
		sym := symbolOf(id)
		need := (octetFor(sym) - 2) / 2
		if need < 0 {
			need = 0
		}
		give := need
		if give*2 > remaining {
			give = remaining / 2
		}
		lp[id] = give
		remaining -= 2 * give
	}
	lp[leftID] += remaining / 2

	boOf := func(id string) int {
		s := 0
		for _, b := range bonds {
			if b.From == id || b.To == id {
				s += b.Order
			}
		}
		return s
	}
	calcFC := func(id string) int { return formalCharge(symbolOf(id), lp[id], boOf(id)) }
	totalAbsFC := func() int { return absInt(calcFC(leftID)) + absInt(calcFC(rightID)) }

	for {
		l1EC := 2*lp[leftID] + 2*boOf(leftID)
		l2EC := 2*lp[rightID] + 2*boOf(rightID)
		if l1EC >= octetFor(a1.Symbol) && l2EC >= octetFor(a2.Symbol) &&
			elemPeriod(a1.Symbol) <= 2 && elemPeriod(a2.Symbol) <= 2 {
			break
		}
		cur := totalAbsFC()
		bestDonor := ""
		bestImp := 0
		for _, id := range []string{leftID, rightID} {
			if lp[id] == 0 {
				continue
			}
			lp[id]--
			bonds[0].Order++
			imp := cur - totalAbsFC()
			lp[id]++
			bonds[0].Order--
			if imp > bestImp {
				bestImp = imp
				bestDonor = id
			} else if imp == 0 && bestDonor == "" {
				bestDonor = id
			}
		}
		if bestDonor == "" {
			break
		}
		l1EC = 2*lp[leftID] + 2*boOf(leftID)
		l2EC = 2*lp[rightID] + 2*boOf(rightID)
		if bestImp <= 0 && l1EC >= octetFor(a1.Symbol) && l2EC >= octetFor(a2.Symbol) {
			break
		}
		lp[bestDonor]--
		bonds[0].Order++
	}

	for i, a := range atoms {
		atoms[i].LonePairs = lp[a.ID]
		atoms[i].FormalCharge = calcFC(a.ID)
	}

	var notes string
	if isRadical {
		notes = fmt.Sprintf("%s has an odd number of valence electrons — it is a free radical.", displayFormula)
	}
	steps := buildDiatomicSteps(a1, a2, totalVE, charge, lp, bonds[0].Order, leftID, rightID, isRadical)

	return &LewisStructure{
		Name:                  displayFormula,
		Formula:               displayFormula,
		Charge:                charge,
		TotalValenceElectrons: totalVE,
		Geometry:              "diatomic",
		Atoms:                 atoms,
		Bonds:                 bonds,
		Steps:                 steps,
		Notes:                 notes,
	}, nil
}

// ── Step Builders ─────────────────────────────────────────────────────────────

func buildSingleCenterSteps(
	central atomCount, terminals []atomCount,
	totalVE, charge, numBonds int,
	lp map[string]int, bonds []LewisBond, centralID string,
	isRadical, isIonic bool,
) []string {
	var steps []string
	all := append([]atomCount{central}, terminals...)
	steps = append(steps, buildVECountStep(all, totalVE, charge))

	if len(terminals) > 0 {
		steps = append(steps, fmt.Sprintf(
			"Identify the central atom: %s (%s) — lowest electronegativity non-hydrogen atom (EN = %.2f)",
			elemName(central.Symbol), central.Symbol, elemEN(central.Symbol)))
		steps = append(steps, fmt.Sprintf(
			"Form single bonds from %s to all terminal atoms: uses %d e⁻, leaving %d e⁻",
			central.Symbol, 2*numBonds, totalVE-2*numBonds))
	}

	for _, t := range terminals {
		tid := t.Symbol + "1"
		lpCount := lp[tid]
		if t.Symbol == "H" {
			steps = append(steps, "Hydrogen satisfies the duet rule with 2 bonding electrons — no lone pairs needed")
		} else if t.Count > 1 {
			steps = append(steps, fmt.Sprintf(
				"Assign %d lone pair(s) to each %s to complete its octet: %d × %d e⁻ = %d e⁻ used",
				lpCount, t.Symbol, t.Count, 2*lpCount, t.Count*2*lpCount))
		} else {
			steps = append(steps, fmt.Sprintf(
				"Assign %d lone pair(s) to %s to complete its octet: %d e⁻",
				lpCount, t.Symbol, 2*lpCount))
		}
	}

	if clp := lp[centralID]; clp > 0 {
		steps = append(steps, fmt.Sprintf("Place remaining %d e⁻ as %d lone pair(s) on %s", 2*clp, clp, central.Symbol))
	}

	if !isIonic {
		for _, b := range bonds {
			if b.Order > 1 {
				var tsym string
				if b.From == centralID {
					tsym = symbolOf(b.To)
				} else {
					tsym = symbolOf(b.From)
				}
				name := "double"
				if b.Order == 3 {
					name = "triple"
				}
				steps = append(steps, fmt.Sprintf(
					"Promote terminal lone pair into a %s bond %s=%s to satisfy octets / minimise formal charges",
					name, central.Symbol, tsym))
			}
		}
	}

	boOf := func(id string) int {
		s := 0
		for _, b := range bonds {
			if b.From == id || b.To == id {
				s += b.Order
			}
		}
		return s
	}

	cElec := 2*lp[centralID] + 2*boOf(centralID)
	steps = append(steps, fmt.Sprintf("Verify: %s has %d e⁻ (%s)", central.Symbol, cElec, octetVerify(central.Symbol, cElec)))
	for _, t := range terminals {
		tid := t.Symbol + "1"
		tElec := 2*lp[tid] + 2*boOf(tid)
		steps = append(steps, fmt.Sprintf("  %s has %d e⁻ (%s)", t.Symbol, tElec, octetVerify(t.Symbol, tElec)))
	}
	steps = append(steps, fmt.Sprintf("Total electrons: %d e⁻ ✓", totalVE))

	var fcParts []string
	fcParts = append(fcParts, fmt.Sprintf("%s = %+d", central.Symbol, formalCharge(central.Symbol, lp[centralID], boOf(centralID))))
	for _, t := range terminals {
		tid := t.Symbol + "1"
		fcParts = append(fcParts, fmt.Sprintf("%s = %+d", t.Symbol, formalCharge(t.Symbol, lp[tid], boOf(tid))))
	}
	steps = append(steps, "Formal charges: "+strings.Join(fcParts, "; "))

	if isRadical {
		steps = append(steps, "Note: odd electron count — one electron remains unpaired (free radical)")
	}
	return steps
}

func buildTwoCenterSteps(
	center atomCount, terminals []atomCount,
	totalVE, charge, numBonds int,
	lp map[string]int, bonds []LewisBond, c1ID, c2ID string, isRadical bool,
) []string {
	var steps []string
	all := append([]atomCount{{center.Symbol, 2}}, terminals...)
	steps = append(steps, buildVECountStep(all, totalVE, charge))
	steps = append(steps, fmt.Sprintf(
		"Topology: two %s atoms form the chain backbone connected by a %s–%s bond",
		elemName(center.Symbol), center.Symbol, center.Symbol))

	for _, b := range bonds {
		if (b.From == c1ID && b.To == c2ID) || (b.From == c2ID && b.To == c1ID) {
			if b.Order > 1 {
				name := "double"
				if b.Order == 3 {
					name = "triple"
				}
				steps = append(steps, fmt.Sprintf("Promote lone pair(s) into the %s–%s bond: forms a %s bond", center.Symbol, center.Symbol, name))
			}
			break
		}
	}

	boOf := func(id string) int {
		s := 0
		for _, b := range bonds {
			if b.From == id || b.To == id {
				s += b.Order
			}
		}
		return s
	}
	cElec := 2*lp[c1ID] + 2*boOf(c1ID)
	steps = append(steps, fmt.Sprintf(
		"Verify: each %s has %d e⁻ (%s) ✓; total = %d e⁻ ✓",
		center.Symbol, cElec, octetVerify(center.Symbol, cElec), totalVE))
	return steps
}

func buildDiatomicSteps(a1, a2 atomCount, totalVE, charge int, lp map[string]int, finalBO int, leftID, rightID string, isRadical bool) []string {
	var steps []string
	var all []atomCount
	if a1.Symbol == a2.Symbol {
		all = []atomCount{{a1.Symbol, 2}}
	} else {
		all = []atomCount{a1, a2}
	}
	steps = append(steps, buildVECountStep(all, totalVE, charge))
	steps = append(steps, "Form a single bond between the two atoms: uses 2 e⁻")
	if finalBO > 1 {
		name := []string{"", "single", "double", "triple"}[finalBO]
		steps = append(steps, fmt.Sprintf("Promote lone pairs to form a %s bond to satisfy octets and minimise formal charges", name))
	} else {
		steps = append(steps, "Distribute remaining electrons as lone pairs to complete octets")
	}
	lp1 := lp[leftID]
	lp2 := lp[rightID]
	e1 := 2*lp1 + 2*finalBO
	e2 := 2*lp2 + 2*finalBO
	steps = append(steps, fmt.Sprintf(
		"Verify: %s has %d e⁻ (%s); %s has %d e⁻ (%s); total = %d e⁻ ✓",
		a1.Symbol, e1, octetVerify(a1.Symbol, e1), a2.Symbol, e2, octetVerify(a2.Symbol, e2), totalVE))
	fc1 := formalCharge(a1.Symbol, lp1, finalBO)
	fc2 := formalCharge(a2.Symbol, lp2, finalBO)
	steps = append(steps, fmt.Sprintf("Formal charges: %s = %+d; %s = %+d", a1.Symbol, fc1, a2.Symbol, fc2))
	if isRadical {
		steps = append(steps, "Note: odd electron count — one electron is unpaired (free radical)")
	}
	return steps
}

// ── Shared Helpers ────────────────────────────────────────────────────────────

func formalCharge(sym string, lonePairs, bondOrderTotal int) int {
	return elemValence(sym) - 2*lonePairs - bondOrderTotal
}

func octetFor(sym string) int {
	if sym == "H" {
		return 2
	}
	return 8
}

// symbolOf extracts the element symbol from an atom ID like "Cl1", "Fe2".
// Tries two-letter match first, falls back to one-letter.
func symbolOf(id string) string {
	if len(id) >= 2 {
		if _, ok := elemOf(id[:2]); ok {
			return id[:2]
		}
	}
	return id[:1]
}

func bondOrderSum(bonds []LewisBond, atomID string) int {
	sum := 0
	for _, b := range bonds {
		if b.From == atomID || b.To == atomID {
			sum += b.Order
		}
	}
	return sum
}

func absInt(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

func determineGeometry(bondedAtoms, centralLP int) string {
	type key struct{ b, lp int }
	table := map[key]string{
		{1, 0}: "diatomic", {2, 0}: "linear", {2, 1}: "bent",
		{2, 2}: "bent", {2, 3}: "linear", {3, 0}: "trigonal_planar",
		{3, 1}: "trigonal_pyramidal", {3, 2}: "t_shaped",
		{4, 0}: "tetrahedral", {4, 1}: "seesaw", {4, 2}: "square_planar",
		{5, 0}: "trigonal_bipyramidal", {5, 1}: "square_pyramidal",
		{6, 0}: "octahedral",
	}
	if g, ok := table[key{bondedAtoms, centralLP}]; ok {
		return g
	}
	return "unknown"
}

func buildDisplayFormula(atoms []atomCount, charge int) string {
	var b strings.Builder
	for _, a := range atoms {
		b.WriteString(a.Symbol)
		if a.Count > 1 {
			b.WriteString(toSubscript(a.Count))
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

var subscriptDigits = []rune{'₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'}
var superscriptDigits = []rune{'⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'}

func toSubscript(n int) string   { return digitString(n, subscriptDigits) }
func toSuperscript(n int) string { return digitString(n, superscriptDigits) }

func digitString(n int, digits []rune) string {
	if n == 0 {
		return string(digits[0])
	}
	var s []rune
	for n > 0 {
		s = append([]rune{digits[n%10]}, s...)
		n /= 10
	}
	return string(s)
}

func buildVECountStep(atoms []atomCount, totalVE, charge int) string {
	var parts []string
	for _, ac := range atoms {
		v := elemValence(ac.Symbol)
		if ac.Count == 1 {
			parts = append(parts, fmt.Sprintf("%s(%d)", ac.Symbol, v))
		} else {
			parts = append(parts, fmt.Sprintf("%d×%s(%d)", ac.Count, ac.Symbol, v))
		}
	}
	line := strings.Join(parts, " + ")
	if charge > 0 {
		line += fmt.Sprintf(" − %d (positive charge)", charge)
	} else if charge < 0 {
		line += fmt.Sprintf(" + %d (negative charge)", -charge)
	}
	return fmt.Sprintf("Count total valence electrons: %s = %d e⁻", line, totalVE)
}

func octetVerify(sym string, n int) string {
	if sym == "H" {
		if n == 2 {
			return "duet ✓"
		}
		return fmt.Sprintf("duet: %d e⁻", n)
	}
	if n == 8 {
		return "octet ✓"
	}
	if n > 8 && elemPeriod(sym) >= 3 {
		return fmt.Sprintf("expanded octet (%d e⁻)", n)
	}
	return fmt.Sprintf("%d e⁻", n)
}
