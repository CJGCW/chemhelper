// Package render fetches 2D compound structures from PubChem and converts
// them to dark-mode-compatible skeletal SVG using the molecule's SDF data.
package render

import (
	"fmt"
	"io"
	"math"
	"net/http"
	"net/url"
	"strconv"
	"strings"
	"sync"
	"time"
)

const pubchemBase = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles"

var httpClient = &http.Client{Timeout: 15 * time.Second}

// RenderOptions controls optional rendering features.
type RenderOptions struct {
	ShowLonePairs bool
}

type cacheKey struct {
	smiles        string
	w, h          int
	showLonePairs bool
}

var (
	mu    sync.RWMutex
	cache = make(map[cacheKey]cacheEntry)
)

type cacheEntry struct {
	svg string
	err error
}

// Render returns a skeletal SVG for the given SMILES string scaled to (width × height).
// Results are cached in memory; PubChem is only contacted once per unique (smiles, w, h, opts).
func Render(smiles string, width, height int, opts RenderOptions) (string, error) {
	k := cacheKey{smiles, width, height, opts.ShowLonePairs}

	mu.RLock()
	if e, ok := cache[k]; ok {
		mu.RUnlock()
		return e.svg, e.err
	}
	mu.RUnlock()

	svg, err := renderUncached(smiles, width, height, opts)

	mu.Lock()
	cache[k] = cacheEntry{svg, err}
	mu.Unlock()

	return svg, err
}

func renderUncached(smiles string, width, height int, opts RenderOptions) (string, error) {
	// Custom-label SMILES ([R], [Nu], [LG] etc.) cannot go through PubChem.
	if needsLocalRenderer(smiles) {
		return renderLocalSMILES(smiles, width, height, opts.ShowLonePairs)
	}

	sdfData, err := fetchSDF(smiles)
	if err != nil {
		return "", err
	}
	mol, err := parseSDF(sdfData)
	if err != nil {
		return "", fmt.Errorf("parse SDF: %w", err)
	}
	svg := toSVG(mol, width, height)
	if opts.ShowLonePairs {
		svg = addLonePairsToPubChemSVG(svg, mol, width, height)
	}
	return svg, nil
}

func fetchSDF(smiles string) (string, error) {
	encoded := url.PathEscape(smiles)
	u := fmt.Sprintf("%s/%s/SDF?record_type=2d", pubchemBase, encoded)

	resp, err := httpClient.Get(u)
	if err != nil {
		return "", fmt.Errorf("PubChem request failed: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode == http.StatusNotFound {
		return "", fmt.Errorf("compound not found in PubChem for SMILES %q", smiles)
	}
	if resp.StatusCode != http.StatusOK {
		return "", fmt.Errorf("PubChem returned status %d", resp.StatusCode)
	}

	data, err := io.ReadAll(resp.Body)
	if err != nil {
		return "", fmt.Errorf("read PubChem response: %w", err)
	}
	return string(data), nil
}

// ── SDF parsing ────────────────────────────────────────────────────────────────

type atom struct {
	X, Y    float64
	Symbol  string
	Charge  int
}

type bond struct {
	A1, A2 int     // 0-indexed
	Order  int     // 1, 2, 3, 4(aromatic)
	Stereo int     // 0=none, 1=Up(wedge), 6=Down(dash)
}

type molecule struct {
	Atoms []atom
	Bonds []bond
}

func parseSDF(sdf string) (*molecule, error) {
	// Normalise line endings
	sdf = strings.ReplaceAll(sdf, "\r\n", "\n")
	lines := strings.Split(sdf, "\n")

	if len(lines) < 4 {
		return nil, fmt.Errorf("SDF too short (%d lines)", len(lines))
	}

	// Line index 3 is the counts line (0-indexed)
	counts := lines[3]
	if len(counts) < 6 {
		return nil, fmt.Errorf("invalid counts line: %q", counts)
	}

	nAtoms, err := strconv.Atoi(strings.TrimSpace(counts[0:3]))
	if err != nil {
		return nil, fmt.Errorf("parse atom count: %w", err)
	}
	nBonds, err := strconv.Atoi(strings.TrimSpace(counts[3:6]))
	if err != nil {
		return nil, fmt.Errorf("parse bond count: %w", err)
	}
	if len(lines) < 4+nAtoms+nBonds {
		return nil, fmt.Errorf("SDF has %d lines, need at least %d", len(lines), 4+nAtoms+nBonds)
	}

	mol := &molecule{
		Atoms: make([]atom, 0, nAtoms),
		Bonds: make([]bond, 0, nBonds),
	}

	// Atom block starts at line index 4
	for i := 0; i < nAtoms; i++ {
		line := lines[4+i]
		if len(line) < 34 {
			return nil, fmt.Errorf("atom line %d too short: %q", i, line)
		}
		x, err := strconv.ParseFloat(strings.TrimSpace(line[0:10]), 64)
		if err != nil {
			return nil, fmt.Errorf("atom %d x coord: %w", i, err)
		}
		y, err := strconv.ParseFloat(strings.TrimSpace(line[10:20]), 64)
		if err != nil {
			return nil, fmt.Errorf("atom %d y coord: %w", i, err)
		}
		sym := strings.TrimSpace(line[31:34])

		var charge int
		if len(line) >= 39 {
			if cc, err2 := strconv.Atoi(strings.TrimSpace(line[36:39])); err2 == nil {
				switch cc {
				case 1:
					charge = 3
				case 2:
					charge = 2
				case 3:
					charge = 1
				case 5:
					charge = -1
				case 6:
					charge = -2
				case 7:
					charge = -3
				}
			}
		}

		mol.Atoms = append(mol.Atoms, atom{X: x, Y: y, Symbol: sym, Charge: charge})
	}

	// Bond block
	for i := 0; i < nBonds; i++ {
		line := lines[4+nAtoms+i]
		if len(line) < 9 {
			return nil, fmt.Errorf("bond line %d too short: %q", i, line)
		}
		a1, err := strconv.Atoi(strings.TrimSpace(line[0:3]))
		if err != nil {
			return nil, fmt.Errorf("bond %d a1: %w", i, err)
		}
		a2, err := strconv.Atoi(strings.TrimSpace(line[3:6]))
		if err != nil {
			return nil, fmt.Errorf("bond %d a2: %w", i, err)
		}
		order, err := strconv.Atoi(strings.TrimSpace(line[6:9]))
		if err != nil {
			return nil, fmt.Errorf("bond %d type: %w", i, err)
		}
		var stereo int
		if len(line) >= 12 {
			stereo, _ = strconv.Atoi(strings.TrimSpace(line[9:12]))
		}
		mol.Bonds = append(mol.Bonds, bond{
			A1: a1 - 1, A2: a2 - 1, // 1-indexed → 0-indexed
			Order:  order,
			Stereo: stereo,
		})
	}

	// Apply M  CHG lines (override atom charge for charged atoms)
	for _, line := range lines {
		if !strings.HasPrefix(line, "M  CHG") {
			continue
		}
		parts := strings.Fields(line)
		// parts: ["M", "CHG", n, a1, c1, a2, c2, ...]
		if len(parts) < 3 {
			continue
		}
		n, err := strconv.Atoi(parts[2])
		if err != nil {
			continue
		}
		for j := 0; j < n && 3+j*2+1 < len(parts); j++ {
			idx, err1 := strconv.Atoi(parts[3+j*2])
			chg, err2 := strconv.Atoi(parts[3+j*2+1])
			if err1 != nil || err2 != nil {
				continue
			}
			if idx-1 < len(mol.Atoms) {
				mol.Atoms[idx-1].Charge = chg
			}
		}
	}

	return mol, nil
}

// ── SVG generation ─────────────────────────────────────────────────────────────

const svgPad = 18.0 // pixels of padding around the structure

func toSVG(mol *molecule, width, height int) string {
	if len(mol.Atoms) == 0 {
		return fmt.Sprintf(`<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg"></svg>`, width, height)
	}

	// Bounding box of atom coordinates
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

	// Centre the molecule
	cx := float64(width)/2 - ((minX+maxX)/2)*scale
	cy := float64(height)/2 + ((minY+maxY)/2)*scale // flip Y: SDF y-up, SVG y-down

	toSX := func(x float64) float64 { return cx + x*scale }
	toSY := func(y float64) float64 { return cy - y*scale }

	// Double-bond parallel separation (capped for very small/large scales)
	bondSep := math.Max(1.5, math.Min(4.0, scale*0.08))

	// Font size for heteroatom labels
	fontSize := math.Max(10, math.Min(14, scale*0.5))

	// Radius to shorten bonds by when an endpoint is a visible atom label
	labelRadius := fontSize * 0.62

	var sb strings.Builder
	sb.WriteString(fmt.Sprintf(
		`<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d">`,
		width, height, width, height,
	))
	sb.WriteString(`<g stroke="currentColor" stroke-linecap="round" stroke-linejoin="round">`)

	// Build a map of implicit-H counts per heavy atom for label rendering.
	hCount := make([]int, len(mol.Atoms))
	for _, b := range mol.Bonds {
		if b.A1 >= len(mol.Atoms) || b.A2 >= len(mol.Atoms) {
			continue
		}
		if mol.Atoms[b.A2].Symbol == "H" {
			hCount[b.A1]++
		} else if mol.Atoms[b.A1].Symbol == "H" {
			hCount[b.A2]++
		}
	}

	// Draw bonds — skip any bond where either endpoint is H (skeletal style hides C–H bonds).
	for _, b := range mol.Bonds {
		if b.A1 >= len(mol.Atoms) || b.A2 >= len(mol.Atoms) {
			continue
		}
		a1 := mol.Atoms[b.A1]
		a2 := mol.Atoms[b.A2]
		if a1.Symbol == "H" || a2.Symbol == "H" {
			continue
		}

		x1, y1 := toSX(a1.X), toSY(a1.Y)
		x2, y2 := toSX(a2.X), toSY(a2.Y)

		// Shorten endpoints that have visible atom labels
		r1 := 0.0
		if isVisible(a1.Symbol) {
			r1 = labelRadius
		}
		r2 := 0.0
		if isVisible(a2.Symbol) {
			r2 = labelRadius
		}
		x1, y1, x2, y2 = shorten(x1, y1, x2, y2, r1, r2)

		switch {
		case b.Stereo == 1: // solid wedge
			drawWedge(&sb, x1, y1, x2, y2, false)
		case b.Stereo == 6: // dashed wedge
			drawWedge(&sb, x1, y1, x2, y2, true)
		case b.Order == 2:
			drawDouble(&sb, x1, y1, x2, y2, bondSep)
		case b.Order == 3:
			drawTriple(&sb, x1, y1, x2, y2, bondSep)
		default:
			sb.WriteString(fmt.Sprintf(
				`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.5"/>`,
				x1, y1, x2, y2,
			))
		}
	}

	// Draw atom labels for heteroatoms (not C, not H).
	// Append implicit H count as subscript so "O" with 1 H renders as "OH", etc.
	subscripts := []string{"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"}
	for ai, a := range mol.Atoms {
		if !isVisible(a.Symbol) {
			continue
		}
		sx, sy := toSX(a.X), toSY(a.Y)
		label := a.Symbol
		if n := hCount[ai]; n == 1 {
			label += "H"
		} else if n > 1 && n < len(subscripts) {
			label += "H" + subscripts[n]
		}
		if a.Charge != 0 {
			if a.Charge == 1 {
				label += "⁺"
			} else if a.Charge == -1 {
				label += "⁻"
			} else if a.Charge > 0 {
				label += fmt.Sprintf("%d+", a.Charge)
			} else {
				label += fmt.Sprintf("%d-", -a.Charge)
			}
		}
		sb.WriteString(fmt.Sprintf(
			`<text x="%.2f" y="%.2f" text-anchor="middle" dominant-baseline="central" `+
				`fill="currentColor" stroke="none" font-family="sans-serif" font-size="%.1f" font-weight="600">%s</text>`,
			sx, sy, fontSize, label,
		))
	}

	sb.WriteString(`</g></svg>`)
	return sb.String()
}

// isVisible returns true for atoms that get an explicit text label in skeletal style.
func isVisible(sym string) bool {
	return sym != "C" && sym != "H" && sym != ""
}

func shorten(x1, y1, x2, y2, r1, r2 float64) (float64, float64, float64, float64) {
	if r1 == 0 && r2 == 0 {
		return x1, y1, x2, y2
	}
	vx, vy := x2-x1, y2-y1
	d := math.Sqrt(vx*vx + vy*vy)
	if d < 1 {
		return x1, y1, x2, y2
	}
	ux, uy := vx/d, vy/d
	return x1 + ux*r1, y1 + uy*r1, x2 - ux*r2, y2 - uy*r2
}

func perp(x1, y1, x2, y2 float64) (float64, float64) {
	vx, vy := x2-x1, y2-y1
	d := math.Sqrt(vx*vx + vy*vy)
	if d < 1e-6 {
		return 0, 1
	}
	return -vy / d, vx / d
}

func drawDouble(sb *strings.Builder, x1, y1, x2, y2, sep float64) {
	px, py := perp(x1, y1, x2, y2)
	sb.WriteString(fmt.Sprintf(
		`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.5"/>`,
		x1+px*sep, y1+py*sep, x2+px*sep, y2+py*sep,
	))
	sb.WriteString(fmt.Sprintf(
		`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.5"/>`,
		x1-px*sep, y1-py*sep, x2-px*sep, y2-py*sep,
	))
}

func drawTriple(sb *strings.Builder, x1, y1, x2, y2, sep float64) {
	px, py := perp(x1, y1, x2, y2)
	// Centre line
	sb.WriteString(fmt.Sprintf(
		`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.5"/>`,
		x1, y1, x2, y2,
	))
	// Outer lines
	sb.WriteString(fmt.Sprintf(
		`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.5"/>`,
		x1+px*sep*1.8, y1+py*sep*1.8, x2+px*sep*1.8, y2+py*sep*1.8,
	))
	sb.WriteString(fmt.Sprintf(
		`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.5"/>`,
		x1-px*sep*1.8, y1-py*sep*1.8, x2-px*sep*1.8, y2-py*sep*1.8,
	))
}

func drawWedge(sb *strings.Builder, x1, y1, x2, y2 float64, isDash bool) {
	vx, vy := x2-x1, y2-y1
	d := math.Sqrt(vx*vx + vy*vy)
	if d < 1 {
		return
	}
	px, py := -vy/d, vx/d // perpendicular unit vector

	const tipHW = 4.0 // half-width at the wide end

	if isDash {
		const nLines = 6
		for i := 0; i <= nLines; i++ {
			t := float64(i) / float64(nLines)
			mx := x1 + vx*t
			my := y1 + vy*t
			hw := t * tipHW
			sb.WriteString(fmt.Sprintf(
				`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-width="1.3"/>`,
				mx-px*hw, my-py*hw, mx+px*hw, my+py*hw,
			))
		}
	} else {
		sb.WriteString(fmt.Sprintf(
			`<polygon points="%.2f,%.2f %.2f,%.2f %.2f,%.2f" fill="currentColor" stroke="none"/>`,
			x1, y1,
			x2+px*tipHW, y2+py*tipHW,
			x2-px*tipHW, y2-py*tipHW,
		))
	}
}
