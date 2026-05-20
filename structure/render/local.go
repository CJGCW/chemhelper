// Local SMILES renderer for custom-label atoms ([R], [Nu], [LG], [X], [E], etc.)
// that PubChem cannot process.  Also used for simple charged species when
// showLonePairs is true, so the layout coordinates are available for dot placement.
package render

import (
	"fmt"
	"math"
	"strings"
)

// ── Custom-label detection ────────────────────────────────────────────────────

// customPrefixes are bracket atom symbols that PubChem doesn't know about.
var customPrefixes = []string{"R]", "Nu]", "LG]", "X]", "E+]", "Nu-]", "E]"}

// needsLocalRenderer reports whether the SMILES contains custom placeholder atoms.
func needsLocalRenderer(smiles string) bool {
	for _, p := range customPrefixes {
		if strings.Contains(smiles, "["+p) {
			return true
		}
	}
	return false
}

// ── Local molecule types ──────────────────────────────────────────────────────

type localAtom struct {
	symbol    string  // display label (e.g. "R", "Nu", "O", "Br")
	charge    int     // formal charge
	hCount    int     // explicit H count from bracket (for display, e.g. H3O+)
	x, y      float64 // layout position in internal units
	placed    bool
	bracketed bool
	aromatic  bool
}

type localBond struct {
	a1, a2  int
	order   int
	aromatic bool
}

type localMol struct {
	atoms []localAtom
	bonds []localBond
	adj   [][]int
}

// ── SMILES tokeniser ──────────────────────────────────────────────────────────

func parseSMILESLocal(smiles string) *localMol {
	var atoms []localAtom
	var bonds []localBond

	type stackEntry struct{ atomIdx, nextBond int }
	stack := []stackEntry{}
	ringMap := map[int]int{}
	lastAtom := -1
	nextBond := 1

	addAtom := func(sym string, chg, hcnt int, bracketed bool) int {
		idx := len(atoms)
		atoms = append(atoms, localAtom{symbol: sym, charge: chg, hCount: hcnt, bracketed: bracketed})
		if lastAtom >= 0 {
			bonds = append(bonds, localBond{a1: lastAtom, a2: idx, order: nextBond})
		}
		lastAtom = idx
		nextBond = 1
		return idx
	}

	i := 0
	for i < len(smiles) {
		c := smiles[i]
		switch {
		case c == '[':
			j := i + 1
			for j < len(smiles) && smiles[j] != ']' {
				j++
			}
			if j >= len(smiles) {
				return nil
			}
			sym, chg, hcnt := parseBracket(smiles[i+1 : j])
			addAtom(sym, chg, hcnt, true)
			i = j + 1

		case c == '(':
			stack = append(stack, stackEntry{lastAtom, nextBond})
			i++

		case c == ')':
			if len(stack) > 0 {
				top := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				lastAtom = top.atomIdx
				nextBond = top.nextBond
			}
			i++

		case c == '=':
			nextBond = 2
			i++
		case c == '#':
			nextBond = 3
			i++
		case c == '-':
			nextBond = 1
			i++
		case c == ':':
			nextBond = 1
			i++

		case c >= '1' && c <= '9':
			digit := int(c - '0')
			if prev, ok := ringMap[digit]; ok {
				bonds = append(bonds, localBond{a1: prev, a2: lastAtom, order: nextBond})
				delete(ringMap, digit)
			} else {
				ringMap[digit] = lastAtom
			}
			nextBond = 1
			i++

		case c == 'C' && i+1 < len(smiles) && smiles[i+1] == 'l':
			addAtom("Cl", 0, 0, false)
			i += 2

		case c == 'B' && i+1 < len(smiles) && smiles[i+1] == 'r':
			addAtom("Br", 0, 0, false)
			i += 2

		case isUpperLetter(c):
			sym := string(c)
			if i+1 < len(smiles) && isLowerLetter(smiles[i+1]) {
				sym += string(smiles[i+1])
				i++
			}
			addAtom(sym, 0, 0, false)
			i++

		case isLowerLetter(c): // aromatic atom
			idx := addAtom(strings.ToUpper(string(c)), 0, 0, false)
			atoms[idx].aromatic = true
			i++

		default:
			i++
		}
	}

	// Post-process: mark aromatic bonds
	for i := range bonds {
		if atoms[bonds[i].a1].aromatic && atoms[bonds[i].a2].aromatic {
			bonds[i].aromatic = true
		}
	}

	n := len(atoms)
	adj := make([][]int, n)
	for _, b := range bonds {
		if b.a1 < n && b.a2 < n {
			adj[b.a1] = append(adj[b.a1], b.a2)
			adj[b.a2] = append(adj[b.a2], b.a1)
		}
	}
	return &localMol{atoms: atoms, bonds: bonds, adj: adj}
}

// parseBracket parses the inside of [...] returning (symbol, charge, hCount).
// Handles: R, Nu, LG, X, E+, Nu-, H3O+, OH-, NH4+, CH3+, BH4-, etc.
func parseBracket(inner string) (sym string, charge, hCount int) {
	i := 0

	// Special case: H-prefix pattern like H3O+, H2O, H4N+
	// Detected when: first char is 'H', second char is a digit, followed by an uppercase element.
	if len(inner) > 2 && inner[0] == 'H' && inner[1] >= '1' && inner[1] <= '9' {
		i++ // skip H
		hcnt := 0
		for i < len(inner) && inner[i] >= '0' && inner[i] <= '9' {
			hcnt = hcnt*10 + int(inner[i]-'0')
			i++
		}
		if i < len(inner) && isUpperLetter(inner[i]) {
			// Looks like H3O+: hCount=3, element=O
			hCount = hcnt
			start := i
			i++
			for i < len(inner) && isLowerLetter(inner[i]) {
				i++
			}
			sym = inner[start:i]
			// Parse charge from the rest
			for i < len(inner) {
				if inner[i] == '+' {
					charge++
				} else if inner[i] == '-' {
					charge--
				}
				i++
			}
			return sym, charge, hCount
		}
		// Didn't match: reset and fall through
		i = 0
	}

	// General case: collect element symbol.
	// Standard SMILES: one uppercase + optional lowercase (Cl, Br, Na…)
	// Custom labels:   may be all-uppercase (LG, Nu where u happens lowercase, E, X…)
	// Strategy: collect one uppercase + following lowercase as the "element".
	// Then: if the next char is uppercase and is NOT 'H', collect it too (custom label).
	start := i
	if i < len(inner) && isUpperLetter(inner[i]) {
		i++
		for i < len(inner) && isLowerLetter(inner[i]) {
			i++
		}
	}
	sym = inner[start:i]

	// Extend for all-uppercase custom labels (LG, TS, etc.)
	for i < len(inner) && isUpperLetter(inner[i]) && inner[i] != 'H' {
		i++
		for i < len(inner) && isLowerLetter(inner[i]) {
			i++
		}
		sym = inner[start:i]
	}

	// Parse H-count suffix (e.g. OH- → hCount=1, NH4+ → hCount=4)
	if i < len(inner) && inner[i] == 'H' {
		i++
		d := 1
		if i < len(inner) && inner[i] >= '0' && inner[i] <= '9' {
			d = int(inner[i] - '0')
			i++
		}
		hCount = d
	}

	// Parse charge
	for i < len(inner) {
		if inner[i] == '+' {
			charge++
		} else if inner[i] == '-' {
			charge--
		}
		i++
	}
	return sym, charge, hCount
}

func isUpperLetter(c byte) bool { return c >= 'A' && c <= 'Z' }
func isLowerLetter(c byte) bool { return c >= 'a' && c <= 'z' }

// ── Implicit hydrogen computation ─────────────────────────────────────────────

var organicTargetValence = map[string]int{
	"C": 4, "N": 3, "O": 2, "S": 2, "P": 3,
	"F": 1, "Cl": 1, "Br": 1, "I": 1,
}

func isPlaceholder(sym string) bool {
	switch sym {
	case "R", "Nu", "LG", "X", "E":
		return true
	}
	return false
}

func computeImplicitH(mol *localMol) {
	for i := range mol.atoms {
		a := &mol.atoms[i]
		if a.bracketed {
			continue
		}
		targetVal, ok := organicTargetValence[a.symbol]
		if !ok {
			continue
		}
		var bondSum float64
		for _, b := range mol.bonds {
			if b.a1 != i && b.a2 != i {
				continue
			}
			if b.aromatic {
				bondSum += 1.5
			} else {
				bondSum += float64(b.order)
			}
		}
		h := targetVal - int(math.Round(bondSum))
		if h < 0 {
			h = 0
		}
		a.hCount = h
	}
}

// ── 2D Layout ─────────────────────────────────────────────────────────────────

const localBondLen = 1.0

// distributeAnglesRoot gives angles for bonds from the root atom (no parent).
func distributeAnglesRoot(n int) []float64 {
	switch n {
	case 1:
		return []float64{0}
	case 2:
		return []float64{0, math.Pi}
	case 3:
		return []float64{0, 2 * math.Pi / 3, 4 * math.Pi / 3}
	default:
		angles := make([]float64, n)
		for i := range angles {
			angles[i] = float64(i) * 2 * math.Pi / float64(n)
		}
		return angles
	}
}

// distributeAngles gives angles for `n` new bonds from an atom whose parent
// bond came from `forward` (the direction away from parent, i.e. the "free" side).
func distributeAngles(n int, forward float64) []float64 {
	switch n {
	case 1:
		// Chain continuation: 60° zigzag kink from forward
		// Use the "right turn" by default (each step kinks +60°)
		return []float64{forward - math.Pi/3}
	case 2:
		return []float64{
			forward + math.Pi/3,
			forward - math.Pi/3,
		}
	case 3:
		return []float64{
			forward,
			forward + 2*math.Pi/3,
			forward - 2*math.Pi/3,
		}
	default:
		angles := make([]float64, n)
		for i := range angles {
			angles[i] = forward + float64(i)*2*math.Pi/float64(n)
		}
		return angles
	}
}

// ── Collision avoidance helpers ───────────────────────────────────────────────

type pos2D struct{ x, y float64 }

func dist2D(a, b localAtom) float64 {
	dx, dy := a.x-b.x, a.y-b.y
	return math.Sqrt(dx*dx + dy*dy)
}

func bondedPair(mol *localMol, i, j int) bool {
	for _, b := range mol.bonds {
		if (b.a1 == i && b.a2 == j) || (b.a1 == j && b.a2 == i) {
			return true
		}
	}
	return false
}

func subtreeAtoms(mol *localMol, root, excludeParent int) []int {
	result := []int{root}
	visited := map[int]bool{root: true, excludeParent: true}
	queue := []int{root}
	for len(queue) > 0 {
		curr := queue[0]
		queue = queue[1:]
		for _, nb := range mol.adj[curr] {
			if !visited[nb] {
				visited[nb] = true
				result = append(result, nb)
				queue = append(queue, nb)
			}
		}
	}
	return result
}

func savePositions(mol *localMol, indices []int) []pos2D {
	saved := make([]pos2D, len(indices))
	for k, i := range indices {
		saved[k] = pos2D{mol.atoms[i].x, mol.atoms[i].y}
	}
	return saved
}

func restorePositions(mol *localMol, indices []int, saved []pos2D) {
	for k, i := range indices {
		mol.atoms[i].x = saved[k].x
		mol.atoms[i].y = saved[k].y
	}
}

func rotateBranchAround(mol *localMol, indices []int, pivot localAtom, angle float64) {
	c, s := math.Cos(angle), math.Sin(angle)
	for _, i := range indices {
		a := &mol.atoms[i]
		dx := a.x - pivot.x
		dy := a.y - pivot.y
		a.x = pivot.x + dx*c - dy*s
		a.y = pivot.y + dx*s + dy*c
	}
}

func countCollisions(mol *localMol, minDist float64) int {
	count := 0
	n := len(mol.atoms)
	for i := 0; i < n-1; i++ {
		for j := i + 1; j < n; j++ {
			if bondedPair(mol, i, j) {
				continue
			}
			if dist2D(mol.atoms[i], mol.atoms[j]) < minDist {
				count++
			}
		}
	}
	return count
}

func avoidCollisions(mol *localMol, parents []int) {
	const minDist = 0.7 * localBondLen
	tryAngles := []float64{
		2 * math.Pi / 3, -2 * math.Pi / 3,
		math.Pi / 3, -math.Pi / 3,
		math.Pi,
	}
	for range [5]struct{}{} {
		anyFixed := false
		for j := 1; j < len(mol.atoms); j++ {
			pj := parents[j]
			if pj < 0 {
				continue
			}
			for i := 0; i < len(mol.atoms); i++ {
				if i == j || bondedPair(mol, i, j) {
					continue
				}
				if dist2D(mol.atoms[i], mol.atoms[j]) >= minDist {
					continue
				}
				sub := subtreeAtoms(mol, j, pj)
				pivot := mol.atoms[pj]
				saved := savePositions(mol, sub)
				before := countCollisions(mol, minDist)
				best := before
				var bestA float64
				hasBest := false
				for _, a := range tryAngles {
					rotateBranchAround(mol, sub, pivot, a)
					score := countCollisions(mol, minDist)
					if score < best {
						best = score
						bestA = a
						hasBest = true
					}
					restorePositions(mol, sub, saved)
				}
				if hasBest {
					rotateBranchAround(mol, sub, pivot, bestA)
					anyFixed = true
				}
			}
		}
		if !anyFixed {
			break
		}
	}
}

// layoutMol assigns x,y coordinates to each atom using a DFS layout.
// Bond angles use standard organic chemistry conventions:
// chain → 60° zigzag, branches distribute evenly around the "forward" direction.
func layoutMol(mol *localMol) {
	n := len(mol.atoms)
	if n == 0 {
		return
	}
	mol.atoms[0].x = 0
	mol.atoms[0].y = 0
	mol.atoms[0].placed = true

	if n == 1 {
		return
	}

	parents := make([]int, n)
	for i := range parents {
		parents[i] = -1
	}

	var dfs func(atomIdx, parentIdx int, fromAngle float64)
	dfs = func(atomIdx, parentIdx int, fromAngle float64) {
		var toPlace []int
		for _, nb := range mol.adj[atomIdx] {
			if !mol.atoms[nb].placed {
				toPlace = append(toPlace, nb)
			}
		}
		if len(toPlace) == 0 {
			return
		}

		forward := fromAngle + math.Pi
		// Normalise to [0, 2π)
		for forward < 0 {
			forward += 2 * math.Pi
		}
		forward = math.Mod(forward, 2*math.Pi)

		angles := distributeAngles(len(toPlace), forward)

		a := &mol.atoms[atomIdx]
		for i, nb := range toPlace {
			parents[nb] = atomIdx
			angle := angles[i]
			b := &mol.atoms[nb]
			b.x = a.x + localBondLen*math.Cos(angle)
			b.y = a.y + localBondLen*math.Sin(angle)
			b.placed = true
			dfs(nb, atomIdx, angle)
		}
	}

	degree := len(mol.adj[0])
	if degree == 0 {
		return
	}

	// For root, start with first bond going right (0°)
	rootAngles := distributeAnglesRoot(degree)
	a0 := &mol.atoms[0]
	for i, nb := range mol.adj[0] {
		if mol.atoms[nb].placed {
			continue
		}
		parents[nb] = 0
		angle := rootAngles[i]
		b := &mol.atoms[nb]
		b.x = a0.x + localBondLen*math.Cos(angle)
		b.y = a0.y + localBondLen*math.Sin(angle)
		b.placed = true
		dfs(nb, 0, angle)
	}

	avoidCollisions(mol, parents)
}

// ── SVG generation ────────────────────────────────────────────────────────────

const localSVGPad = 10.0

// renderLocalSMILES renders a custom-label SMILES to SVG.
func renderLocalSMILES(smiles string, width, height int, showLonePairs bool) (string, error) {
	mol := parseSMILESLocal(smiles)
	if mol == nil || len(mol.atoms) == 0 {
		return renderLabelFallback(smiles, width, height), nil
	}
	computeImplicitH(mol)
	layoutMol(mol)
	return localToSVG(mol, width, height, showLonePairs), nil
}

// renderLabelFallback renders the SMILES as a plain text label.
func renderLabelFallback(smiles string, width, height int) string {
	return fmt.Sprintf(
		`<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg">
  <text x="%d" y="%d" text-anchor="middle" dominant-baseline="central" fill="currentColor" font-family="monospace" font-size="12">%s</text>
</svg>`,
		width, height, width/2, height/2, svgEscape(smiles))
}

func localToSVG(mol *localMol, width, height int, showLonePairs bool) string {
	n := len(mol.atoms)
	if n == 0 {
		return renderLabelFallback("", width, height)
	}

	// Compute bounding box
	minX, maxX := mol.atoms[0].x, mol.atoms[0].x
	minY, maxY := mol.atoms[0].y, mol.atoms[0].y
	for _, a := range mol.atoms {
		if a.x < minX {
			minX = a.x
		}
		if a.x > maxX {
			maxX = a.x
		}
		if a.y < minY {
			minY = a.y
		}
		if a.y > maxY {
			maxY = a.y
		}
	}

	drawW := float64(width) - 2*localSVGPad
	drawH := float64(height) - 2*localSVGPad
	dx := maxX - minX
	dy := maxY - minY

	var scale float64
	switch {
	case dx == 0 && dy == 0:
		scale = 40
	case dx == 0:
		scale = drawH / dy
	case dy == 0:
		scale = drawW / dx
	default:
		sx := drawW / dx
		sy := drawH / dy
		if sx < sy {
			scale = sx
		} else {
			scale = sy
		}
	}
	scale = math.Min(scale, 55) // cap scale for very small molecules

	toSX := func(x float64) float64 {
		return float64(width)/2 + (x-(minX+maxX)/2)*scale
	}
	toSY := func(y float64) float64 {
		return float64(height)/2 - (y-(minY+maxY)/2)*scale
	}

	fontSize := math.Max(10, math.Min(14, scale*0.45))
	labelRadius := fontSize * 0.65
	bondSepLocal := math.Max(1.5, math.Min(3.5, scale*0.07))

	subscripts := []string{"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"}

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf(
		`<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d">`,
		width, height, width, height,
	))
	sb.WriteString(`<g stroke="currentColor" stroke-linecap="round" stroke-linejoin="round">`)

	// Draw bonds
	for _, b := range mol.bonds {
		if b.a1 >= n || b.a2 >= n {
			continue
		}
		a1 := &mol.atoms[b.a1]
		a2 := &mol.atoms[b.a2]

		hn1 := len(mol.adj[b.a1])
		hn2 := len(mol.adj[b.a2])

		x1, y1 := toSX(a1.x), toSY(a1.y)
		x2, y2 := toSX(a2.x), toSY(a2.y)

		r1 := localLabelRadius(*a1, hn1, labelRadius)
		r2 := localLabelRadius(*a2, hn2, labelRadius)
		x1, y1, x2, y2 = shorten(x1, y1, x2, y2, r1, r2)

		switch b.order {
		case 2:
			drawDouble(&sb, x1, y1, x2, y2, bondSepLocal)
		case 3:
			drawTriple(&sb, x1, y1, x2, y2, bondSepLocal)
		default:
			sb.WriteString(fmt.Sprintf(
				`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.5"/>`,
				x1, y1, x2, y2,
			))
		}
	}

	// Draw atom labels (skeletal-aware: interior C and aromatic C rendered as vertices)
	for i, a := range mol.atoms {
		hn := len(mol.adj[i])
		label := localAtomLabel(a, hn, subscripts)
		if label == "" {
			continue // skeletal vertex: no text needed
		}
		sx, sy := toSX(a.x), toSY(a.y)

		sb.WriteString(fmt.Sprintf(
			`<text x="%.2f" y="%.2f" text-anchor="middle" dominant-baseline="central" `+
				`fill="currentColor" stroke="none" font-family="sans-serif" font-size="%.1f" font-weight="600">%s</text>`,
			sx, sy, fontSize, svgEscape(label),
		))
	}

	// Lone pair dots
	if showLonePairs {
		renderLocalLonePairs(&sb, mol, toSX, toSY, scale)
	}

	sb.WriteString(`</g></svg>`)
	return sb.String()
}

// localLabelRadius returns the radius to shorten bonds at this atom.
func localLabelRadius(a localAtom, heavyNeighborCount int, base float64) float64 {
	lbl := localAtomLabel(a, heavyNeighborCount, nil)
	if len(lbl) == 0 {
		return 0
	}
	extra := float64(len([]rune(lbl))-1) * base * 0.4
	return base + extra
}

// hLabel returns the hydrogen suffix string for a label.
func hLabel(hCount int, subscripts []string) string {
	if hCount <= 0 {
		return ""
	}
	if hCount == 1 {
		return "H"
	}
	if subscripts != nil && hCount < len(subscripts) {
		return "H" + subscripts[hCount]
	}
	return fmt.Sprintf("H%d", hCount)
}

// chargeLabel returns the superscript charge string.
func chargeLabel(charge int) string {
	switch charge {
	case 0:
		return ""
	case 1:
		return "⁺"
	case -1:
		return "⁻"
	case 2:
		return "²⁺"
	case -2:
		return "²⁻"
	default:
		if charge > 0 {
			return fmt.Sprintf("%d+", charge)
		}
		return fmt.Sprintf("%d-", -charge)
	}
}

// localAtomLabel constructs the skeletal-aware display label for an atom.
func localAtomLabel(a localAtom, heavyNeighborCount int, subscripts []string) string {
	sym := a.symbol
	if isPlaceholder(sym) {
		return sym + hLabel(a.hCount, subscripts) + chargeLabel(a.charge)
	}
	if a.aromatic && sym == "C" {
		return "" // aromatic carbon: skeletal vertex
	}
	if sym == "C" {
		if a.hCount == 0 && heavyNeighborCount >= 2 {
			return chargeLabel(a.charge) // interior quaternary C: vertex (just charge if charged)
		}
		return "C" + hLabel(a.hCount, subscripts) + chargeLabel(a.charge)
	}
	return sym + hLabel(a.hCount, subscripts) + chargeLabel(a.charge)
}

// ── Lone pair rendering ───────────────────────────────────────────────────────

// lonePairCount returns the number of lone pairs for an atom based on
// its symbol, formal charge, and bond count.
func lonePairCount(sym string, charge, bondCount int) int {
	valence := map[string]int{
		"O": 6, "S": 6, "Se": 6,
		"N": 5, "P": 5,
		"F": 7, "Cl": 7, "Br": 7, "I": 7,
		"B": 3,
		"C": 4,
	}
	v, ok := valence[sym]
	if !ok {
		return 0
	}
	// LP electrons = valence − charge − bond electrons (each bond uses 1 e from this atom)
	lpElectrons := v - charge - bondCount
	if lpElectrons < 0 {
		lpElectrons = 0
	}
	return lpElectrons / 2
}

// renderLocalLonePairs draws lone pair dots on heteroatoms in the local molecule.
func renderLocalLonePairs(sb *strings.Builder, mol *localMol,
	toSX, toSY func(float64) float64, scale float64) {

	dotR := math.Max(1.2, scale*0.04)
	dotSpacing := dotR * 2.5
	lpDist := math.Max(6, scale*0.28)

	for ai, a := range mol.atoms {
		lp := lonePairCount(a.symbol, a.charge, len(mol.adj[ai]))
		if lp <= 0 {
			continue
		}

		sx, sy := toSX(a.x), toSY(a.y)

		// Compute bond directions (unit vectors pointing away from this atom)
		var bx, by float64
		for _, nb := range mol.adj[ai] {
			dx := toSX(mol.atoms[nb].x) - sx
			dy := toSY(mol.atoms[nb].y) - sy
			d := math.Sqrt(dx*dx + dy*dy)
			if d > 0 {
				bx += dx / d
				by += dy / d
			}
		}

		// "Free" direction = opposite of the bond centroid
		var freeAngles []float64
		if bx == 0 && by == 0 {
			// No bonds (isolated atom): 4 directions
			freeAngles = []float64{math.Pi / 4, 3 * math.Pi / 4, 5 * math.Pi / 4, 7 * math.Pi / 4}
		} else {
			bondAngle := math.Atan2(by, bx)
			// Free direction is opposite
			freeAngle := bondAngle + math.Pi
			// For 2 LP: place above and below the free direction
			if lp >= 2 {
				freeAngles = []float64{freeAngle - math.Pi/4, freeAngle + math.Pi/4}
			} else {
				freeAngles = []float64{freeAngle}
			}
		}

		for lpIdx := 0; lpIdx < lp && lpIdx < len(freeAngles); lpIdx++ {
			angle := freeAngles[lpIdx]
			cx := sx + math.Cos(angle)*lpDist
			cy := sy + math.Sin(angle)*lpDist

			// Draw two dots side by side (perpendicular to direction)
			perpAngle := angle + math.Pi/2
			for _, sign := range []float64{-1, 1} {
				dx := cx + sign*math.Cos(perpAngle)*dotSpacing/2
				dy := cy + sign*math.Sin(perpAngle)*dotSpacing/2
				sb.WriteString(fmt.Sprintf(
					`<circle cx="%.2f" cy="%.2f" r="%.2f" fill="currentColor" stroke="none"/>`,
					dx, dy, dotR,
				))
			}
		}
	}
}

// ── Lone pairs on PubChem SVGs ────────────────────────────────────────────────

// addLonePairsToPubChemSVG appends lone pair dots to an SVG produced from
// PubChem SDF data.  mol contains the atom positions and adjacency info.
func addLonePairsToPubChemSVG(svg string, mol *molecule, width, height int) string {
	if len(mol.Atoms) == 0 {
		return svg
	}

	// Recompute the same coordinate transform used in toSVG
	minX, maxX := mol.Atoms[0].X, mol.Atoms[0].X
	minY, maxY := mol.Atoms[0].Y, mol.Atoms[0].Y
	for _, a := range mol.Atoms {
		if a.X < minX {
			minX = a.X
		}
		if a.X > maxX {
			maxX = a.X
		}
		if a.Y < minY {
			minY = a.Y
		}
		if a.Y > maxY {
			maxY = a.Y
		}
	}
	dx := maxX - minX
	dy := maxY - minY
	drawW := float64(width) - 2*svgPad
	drawH := float64(height) - 2*svgPad

	var scale float64
	switch {
	case dx == 0 && dy == 0:
		scale = 40
	case dx == 0:
		scale = drawH / dy
	case dy == 0:
		scale = drawW / dx
	default:
		sx := drawW / dx
		sy := drawH / dy
		if sx < sy {
			scale = sx
		} else {
			scale = sy
		}
	}

	cx := float64(width)/2 - ((minX+maxX)/2)*scale
	cy := float64(height)/2 + ((minY+maxY)/2)*scale

	toSX := func(x float64) float64 { return cx + x*scale }
	toSY := func(y float64) float64 { return cy - y*scale }

	// Build adjacency list from SDF bonds
	adj := make([][]int, len(mol.Atoms))
	for _, b := range mol.Bonds {
		if b.A1 < len(mol.Atoms) && b.A2 < len(mol.Atoms) {
			adj[b.A1] = append(adj[b.A1], b.A2)
			adj[b.A2] = append(adj[b.A2], b.A1)
		}
	}

	dotR := math.Max(1.2, scale*0.04)
	dotSpacing := dotR * 2.5
	lpDist := math.Max(5, scale*0.25)

	var dots strings.Builder
	for ai, a := range mol.Atoms {
		lp := lonePairCount(a.Symbol, a.Charge, len(adj[ai]))
		if lp <= 0 {
			continue
		}
		sx, sy := toSX(a.X), toSY(a.Y)

		var bxSum, bySum float64
		for _, nb := range adj[ai] {
			ddx := toSX(mol.Atoms[nb].X) - sx
			ddy := toSY(mol.Atoms[nb].Y) - sy
			d := math.Sqrt(ddx*ddx + ddy*ddy)
			if d > 0 {
				bxSum += ddx / d
				bySum += ddy / d
			}
		}

		var freeAngles []float64
		if bxSum == 0 && bySum == 0 {
			freeAngles = []float64{math.Pi / 4, 3 * math.Pi / 4, 5 * math.Pi / 4, 7 * math.Pi / 4}
		} else {
			bondAngle := math.Atan2(bySum, bxSum)
			freeAngle := bondAngle + math.Pi
			if lp >= 2 {
				freeAngles = []float64{freeAngle - math.Pi/4, freeAngle + math.Pi/4}
			} else {
				freeAngles = []float64{freeAngle}
			}
		}

		for lpIdx := 0; lpIdx < lp && lpIdx < len(freeAngles); lpIdx++ {
			angle := freeAngles[lpIdx]
			pcx := sx + math.Cos(angle)*lpDist
			pcy := sy + math.Sin(angle)*lpDist
			perpAngle := angle + math.Pi/2
			for _, sign := range []float64{-1, 1} {
				ddx := pcx + sign*math.Cos(perpAngle)*dotSpacing/2
				ddy := pcy + sign*math.Sin(perpAngle)*dotSpacing/2
				dots.WriteString(fmt.Sprintf(
					`<circle cx="%.2f" cy="%.2f" r="%.2f" fill="currentColor" stroke="none"/>`,
					ddx, ddy, dotR,
				))
			}
		}
	}

	if dots.Len() == 0 {
		return svg
	}

	// Insert dots before closing </g></svg>
	insert := dots.String()
	if idx := strings.LastIndex(svg, "</g>"); idx >= 0 {
		return svg[:idx] + insert + svg[idx:]
	}
	return svg + insert
}

// svgEscape escapes special XML characters in SVG text content.
func svgEscape(s string) string {
	s = strings.ReplaceAll(s, "&", "&amp;")
	s = strings.ReplaceAll(s, "<", "&lt;")
	s = strings.ReplaceAll(s, ">", "&gt;")
	return s
}
