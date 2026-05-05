// Package estimate provides rule-based spectral estimation from SMILES.
// Estimates are based on Brown & Foote Organic Chemistry correlation tables,
// not quantum-mechanical calculation.
package estimate

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"
)

// neighbor is an adjacency entry stored on each atom.
type neighbor struct {
	idx   int
	order float64 // 1.0=single, 1.5=aromatic, 2.0=double, 3.0=triple
}

// atom is a node in the molecular graph.
type atom struct {
	idx         int
	element     string
	aromatic    bool
	charge      int
	hCount      int // total attached H (explicit from brackets, or implicit)
	fromBracket bool
	nbrs        []neighbor
}

// molecule is a parsed SMILES molecular graph.
type molecule struct {
	atoms []atom
}

func (m *molecule) addBond(i, j int, order float64) {
	m.atoms[i].nbrs = append(m.atoms[i].nbrs, neighbor{idx: j, order: order})
	m.atoms[j].nbrs = append(m.atoms[j].nbrs, neighbor{idx: i, order: order})
}

// valenceFor returns the standard valence for an element.
func valenceFor(elem string) int {
	switch elem {
	case "C":
		return 4
	case "N":
		return 3
	case "O":
		return 2
	case "S":
		return 2
	case "P":
		return 3
	case "F", "Cl", "Br", "I":
		return 1
	case "B":
		return 3
	case "Si":
		return 4
	}
	return 0
}

// implicitH computes the implicit H count for atom i after all bonds are added.
func implicitH(m *molecule, i int) int {
	a := &m.atoms[i]
	v := valenceFor(a.element)
	if v == 0 {
		return 0
	}
	var bondSum float64
	for _, nb := range a.nbrs {
		bondSum += nb.order
	}
	h := v - int(math.Round(bondSum)) - a.charge
	if h < 0 {
		return 0
	}
	return h
}

// ringEntry stores state for an open ring closure.
type ringEntry struct {
	atomIdx int
	bond    float64 // -1 = unspecified
}

// smiParser parses a SMILES string character by character.
type smiParser struct {
	smi   string
	pos   int
	mol   *molecule
	prev  int            // most recent atom added (-1 = none yet)
	bond  float64        // pending explicit bond order (-1 = infer)
	stack []int          // branch stack: prev atom at each '(' level
	rings map[int]ringEntry
}

func (p *smiParser) peek() byte {
	if p.pos < len(p.smi) {
		return p.smi[p.pos]
	}
	return 0
}

func (p *smiParser) advance() byte {
	c := p.peek()
	if c != 0 {
		p.pos++
	}
	return c
}

func (p *smiParser) defaultBond(prevIdx, curIdx int) float64 {
	if p.mol.atoms[prevIdx].aromatic && p.mol.atoms[curIdx].aromatic {
		return 1.5
	}
	return 1.0
}

func (p *smiParser) addAtom(elem string, aromatic bool, charge, explicitH int, fromBracket bool) int {
	idx := len(p.mol.atoms)
	p.mol.atoms = append(p.mol.atoms, atom{
		idx:         idx,
		element:     elem,
		aromatic:    aromatic,
		charge:      charge,
		hCount:      explicitH,
		fromBracket: fromBracket,
	})
	return idx
}

func (p *smiParser) parseBracket() (elem string, aromatic bool, charge, explicitH int, err error) {
	p.advance() // skip '['

	// Isotope digits — ignore
	for p.pos < len(p.smi) && p.smi[p.pos] >= '0' && p.smi[p.pos] <= '9' {
		p.advance()
	}

	if p.pos >= len(p.smi) {
		return "", false, 0, 0, fmt.Errorf("unexpected end inside bracket atom")
	}

	c := p.smi[p.pos]
	if c >= 'a' && c <= 'z' {
		elem = strings.ToUpper(string(c))
		aromatic = true
		p.advance()
	} else if c >= 'A' && c <= 'Z' {
		elem = string(c)
		p.advance()
		if p.pos < len(p.smi) {
			c2 := p.smi[p.pos]
			if c2 >= 'a' && c2 <= 'z' {
				// Check for 2-char symbols: Cl, Br, Si, etc.
				two := string(c) + string(c2)
				switch two {
				case "Cl", "Br", "Si", "Na", "Mg", "Ca", "Al":
					elem = two
					p.advance()
				}
			}
		}
	} else {
		// * wildcard or other — treat as C
		elem = "C"
		p.advance()
	}

	// Stereo (@ and @@) — skip
	for p.pos < len(p.smi) && p.smi[p.pos] == '@' {
		p.advance()
	}

	// Explicit H
	if p.pos < len(p.smi) && p.smi[p.pos] == 'H' {
		p.advance()
		explicitH = 1
		if p.pos < len(p.smi) && p.smi[p.pos] >= '2' && p.smi[p.pos] <= '9' {
			n, _ := strconv.Atoi(string(p.smi[p.pos]))
			explicitH = n
			p.advance()
		}
	}

	// Charge
	if p.pos < len(p.smi) && (p.smi[p.pos] == '+' || p.smi[p.pos] == '-') {
		sign := 1
		if p.smi[p.pos] == '-' {
			sign = -1
		}
		p.advance()
		mag := 1
		// Multi-char charge like +2 or ++
		if p.pos < len(p.smi) && p.smi[p.pos] >= '2' && p.smi[p.pos] <= '9' {
			mag, _ = strconv.Atoi(string(p.smi[p.pos]))
			p.advance()
		} else if p.pos < len(p.smi) && (p.smi[p.pos] == '+' || p.smi[p.pos] == '-') {
			// ++ means +2, -- means -2
			mag = 2
			p.advance()
		}
		charge = sign * mag
	}

	// Skip to closing ']'
	for p.pos < len(p.smi) && p.smi[p.pos] != ']' {
		p.advance()
	}
	if p.pos < len(p.smi) {
		p.advance() // skip ']'
	}

	return elem, aromatic, charge, explicitH, nil
}

func (p *smiParser) parseAtom() (elem string, aromatic bool, err error) {
	if p.pos >= len(p.smi) {
		return "", false, fmt.Errorf("unexpected end of SMILES")
	}
	c := p.smi[p.pos]

	// Aromatic organic subset
	switch c {
	case 'c':
		p.advance()
		return "C", true, nil
	case 'n':
		p.advance()
		return "N", true, nil
	case 'o':
		p.advance()
		return "O", true, nil
	case 's':
		p.advance()
		return "S", true, nil
	case 'b':
		p.advance()
		return "B", true, nil
	case 'p':
		p.advance()
		return "P", true, nil
	}

	// 2-char aliphatic elements
	if p.pos+1 < len(p.smi) {
		two := string(p.smi[p.pos : p.pos+2])
		switch two {
		case "Cl", "Br", "Si", "Na", "Mg", "Ca", "Al":
			p.pos += 2
			return two, false, nil
		}
	}

	// Single-char aliphatic
	if c >= 'A' && c <= 'Z' {
		p.advance()
		return string(c), false, nil
	}

	return "", false, fmt.Errorf("unexpected character '%c' at position %d", c, p.pos)
}

func (p *smiParser) handleRingClosure(ringNum int) {
	if entry, exists := p.rings[ringNum]; exists {
		order := p.bond
		if order < 0 {
			if entry.bond >= 0 {
				order = entry.bond
			} else if p.prev >= 0 {
				order = p.defaultBond(entry.atomIdx, p.prev)
			} else {
				order = 1.0
			}
		}
		if p.prev >= 0 {
			p.mol.addBond(entry.atomIdx, p.prev, order)
		}
		delete(p.rings, ringNum)
	} else {
		p.rings[ringNum] = ringEntry{atomIdx: p.prev, bond: p.bond}
	}
	p.bond = -1
}

func (p *smiParser) parse() error {
	p.rings = make(map[int]ringEntry)

	for p.pos < len(p.smi) {
		c := p.smi[p.pos]

		switch {
		case c == '(':
			p.stack = append(p.stack, p.prev)
			p.advance()

		case c == ')':
			if len(p.stack) == 0 {
				return fmt.Errorf("unmatched ')' at position %d", p.pos)
			}
			p.prev = p.stack[len(p.stack)-1]
			p.stack = p.stack[:len(p.stack)-1]
			p.bond = -1
			p.advance()

		case c == '-':
			p.bond = 1.0
			p.advance()

		case c == '=':
			p.bond = 2.0
			p.advance()

		case c == '#':
			p.bond = 3.0
			p.advance()

		case c == ':':
			p.bond = 1.5
			p.advance()

		case c == '/' || c == '\\':
			p.bond = 1.0 // stereochemistry ignored, treat as single
			p.advance()

		case c == '.':
			p.prev = -1
			p.bond = -1
			p.advance()

		case c == '[':
			elem, aromatic, charge, explicitH, err := p.parseBracket()
			if err != nil {
				return err
			}
			idx := p.addAtom(elem, aromatic, charge, explicitH, true)
			if p.prev >= 0 {
				order := p.bond
				if order < 0 {
					order = p.defaultBond(p.prev, idx)
				}
				p.mol.addBond(p.prev, idx, order)
			}
			p.prev = idx
			p.bond = -1

		case c == '%':
			p.advance()
			if p.pos+1 >= len(p.smi) {
				return fmt.Errorf("invalid ring closure '%%' near position %d", p.pos)
			}
			n1, _ := strconv.Atoi(string(p.smi[p.pos]))
			p.advance()
			n2, _ := strconv.Atoi(string(p.smi[p.pos]))
			p.advance()
			p.handleRingClosure(n1*10 + n2)

		case c >= '0' && c <= '9':
			ringNum := int(c - '0')
			p.advance()
			p.handleRingClosure(ringNum)

		default:
			elem, aromatic, err := p.parseAtom()
			if err != nil {
				return err
			}
			idx := p.addAtom(elem, aromatic, 0, 0, false)
			if p.prev >= 0 {
				order := p.bond
				if order < 0 {
					order = p.defaultBond(p.prev, idx)
				}
				p.mol.addBond(p.prev, idx, order)
			}
			p.prev = idx
			p.bond = -1
		}
	}

	// Fill in implicit H for non-bracket atoms
	for i := range p.mol.atoms {
		if !p.mol.atoms[i].fromBracket {
			p.mol.atoms[i].hCount = implicitH(p.mol, i)
		}
	}

	return nil
}

// parseSMILES parses a SMILES string into a molecule.
func parseSMILES(smi string) (*molecule, error) {
	if strings.TrimSpace(smi) == "" {
		return nil, fmt.Errorf("empty SMILES string")
	}
	mol := &molecule{}
	p := &smiParser{smi: smi, prev: -1, bond: -1, mol: mol}
	if err := p.parse(); err != nil {
		return nil, fmt.Errorf("unparseable SMILES: %w", err)
	}
	if len(mol.atoms) == 0 {
		return nil, fmt.Errorf("SMILES produced no atoms")
	}
	return mol, nil
}

// ── Molecular properties ──────────────────────────────────────────────────────

var atomicWeights = map[string]float64{
	"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
	"F": 18.998, "Cl": 35.45, "Br": 79.904, "I": 126.904,
	"S": 32.06, "P": 30.974, "B": 10.81, "Si": 28.085,
}

var nominalMasses = map[string]int{
	"H": 1, "C": 12, "N": 14, "O": 16,
	"F": 19, "Cl": 35, "Br": 79, "I": 127,
	"S": 32, "P": 31, "B": 11, "Si": 28,
}

func (m *molecule) elemCounts() map[string]int {
	counts := make(map[string]int)
	for i := range m.atoms {
		a := &m.atoms[i]
		counts[a.element]++
		counts["H"] += a.hCount
	}
	return counts
}

func (m *molecule) molecularFormula() string {
	counts := m.elemCounts()
	var sb strings.Builder
	if c := counts["C"]; c > 0 {
		sb.WriteString("C")
		if c > 1 {
			sb.WriteString(strconv.Itoa(c))
		}
	}
	if h := counts["H"]; h > 0 {
		sb.WriteString("H")
		if h > 1 {
			sb.WriteString(strconv.Itoa(h))
		}
	}
	others := make([]string, 0, len(counts))
	for k := range counts {
		if k != "C" && k != "H" && counts[k] > 0 {
			others = append(others, k)
		}
	}
	sort.Strings(others)
	for _, k := range others {
		sb.WriteString(k)
		if counts[k] > 1 {
			sb.WriteString(strconv.Itoa(counts[k]))
		}
	}
	return sb.String()
}

func (m *molecule) molecularWeight() float64 {
	var total float64
	for i := range m.atoms {
		a := &m.atoms[i]
		if w, ok := atomicWeights[a.element]; ok {
			total += w
		}
		total += float64(a.hCount) * 1.008
	}
	return total
}

func (m *molecule) nominalMass() int {
	total := 0
	for i := range m.atoms {
		a := &m.atoms[i]
		if w, ok := nominalMasses[a.element]; ok {
			total += w
		}
		total += a.hCount * 1
	}
	return total
}

func (m *molecule) degreesUnsaturation() int {
	counts := m.elemCounts()
	c := counts["C"]
	h := counts["H"]
	n := counts["N"]
	x := counts["F"] + counts["Cl"] + counts["Br"] + counts["I"]
	dou := (2*c + 2 + n - h - x) / 2
	if dou < 0 {
		return 0
	}
	return dou
}

// ── Atom environment signature for NMR equivalence ───────────────────────────

// atomSig builds a topological signature for atom idx up to depth bonds.
// Chemically equivalent atoms share the same signature.
func atomSig(m *molecule, idx int, depth int, visited map[int]bool) string {
	a := &m.atoms[idx]
	elem := a.element
	if a.aromatic {
		elem = strings.ToLower(elem)
	}
	if depth == 0 {
		return elem
	}
	// Copy visited map for this branch
	newVisited := make(map[int]bool, len(visited)+1)
	for k, v := range visited {
		newVisited[k] = v
	}
	newVisited[idx] = true

	sigs := make([]string, 0, len(a.nbrs))
	for _, nb := range a.nbrs {
		if !newVisited[nb.idx] {
			sub := atomSig(m, nb.idx, depth-1, newVisited)
			sigs = append(sigs, fmt.Sprintf("%.0f:%s", nb.order, sub))
		}
	}
	sort.Strings(sigs)
	return fmt.Sprintf("%s(%s)", elem, strings.Join(sigs, ","))
}
