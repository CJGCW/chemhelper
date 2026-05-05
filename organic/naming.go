package organic

import (
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// Carbon represents a single carbon atom in the main chain.
type Carbon struct {
	Position     int           // 1-indexed position in the main chain
	Substituents []Substituent // branch groups attached at this carbon
	BondOrder    int           // 1 = single, 2 = double, 3 = triple
}

// Substituent is a branch group attached to a chain carbon.
type Substituent struct {
	Name     string // IUPAC group name, e.g. "methyl", "ethyl"
	Position int    // position within the substituent (1 for direct attachment)
}

// Chain is the identified longest carbon chain of the molecule.
type Chain struct {
	Carbons []Carbon
}

// NamingResult holds the IUPAC name and step-by-step derivation.
type NamingResult struct {
	Name  string
	Steps []string
}

// chainPrefixes maps carbon count to IUPAC numerical stem.
var chainPrefixes = map[int]string{
	1: "meth", 2: "eth", 3: "prop", 4: "but", 5: "pent",
	6: "hex", 7: "hept", 8: "oct", 9: "non", 10: "dec",
	11: "undec", 12: "dodec",
}

// multiplicityPrefixes maps substituent count to its IUPAC multiplicity prefix.
var multiplicityPrefixes = []string{"", "", "di", "tri", "tetra", "penta", "hexa"}

// stripMultiplicity removes di-/tri-/etc from a name so that alphabetical
// sorting ignores the multiplicity prefix (IUPAC rule).
func stripMultiplicity(name string) string {
	for _, p := range []string{"hexa", "penta", "tetra", "tri", "di"} {
		if strings.HasPrefix(name, p) {
			return strings.TrimPrefix(name, p)
		}
	}
	return name
}

// NameAlkane generates the IUPAC name for a branched alkane given its main chain.
// The caller is responsible for identifying the longest chain and numbering
// carbons to give the lowest locant set to substituents.
func NameAlkane(chain Chain) (NamingResult, error) {
	n := len(chain.Carbons)
	if n == 0 {
		return NamingResult{}, fmt.Errorf("chain must contain at least one carbon")
	}
	stem, ok := chainPrefixes[n]
	if !ok {
		return NamingResult{}, fmt.Errorf("unsupported chain length: %d (max 12)", n)
	}
	parent := stem + "ane"
	steps := []string{
		fmt.Sprintf("Identify longest carbon chain: %d carbons → parent name is %s.", n, parent),
	}

	// Collect (position, name) for every substituent on every carbon.
	type entry struct {
		pos  int
		name string
	}
	var entries []entry
	for _, c := range chain.Carbons {
		for _, s := range c.Substituents {
			entries = append(entries, entry{pos: c.Position, name: s.Name})
		}
	}

	if len(entries) == 0 {
		steps = append(steps, "No substituents — name is "+parent+".")
		return NamingResult{Name: parent, Steps: steps}, nil
	}

	// Group positions by substituent name.
	groups := map[string][]int{}
	for _, e := range entries {
		groups[e.name] = append(groups[e.name], e.pos)
	}
	for name := range groups {
		sort.Ints(groups[name])
	}

	// Record the locant set for the step log.
	allPos := make([]int, 0, len(entries))
	for _, positions := range groups {
		allPos = append(allPos, positions...)
	}
	sort.Ints(allPos)
	posStr := make([]string, len(allPos))
	for i, p := range allPos {
		posStr[i] = strconv.Itoa(p)
	}
	steps = append(steps,
		"Number the chain to give lowest locants to substituents: "+strings.Join(posStr, ",")+
			" (sum chosen over opposite numbering).",
	)

	// Sort substituent names alphabetically, ignoring multiplicity prefixes.
	uniqueNames := make([]string, 0, len(groups))
	for name := range groups {
		uniqueNames = append(uniqueNames, name)
	}
	sort.Slice(uniqueNames, func(i, j int) bool {
		return stripMultiplicity(uniqueNames[i]) < stripMultiplicity(uniqueNames[j])
	})
	steps = append(steps,
		"Alphabetise substituents (ignoring di-/tri- prefixes): "+strings.Join(uniqueNames, ", ")+".",
	)

	// Build each substituent prefix, e.g. "2,2-dimethyl".
	var parts []string
	for _, name := range uniqueNames {
		positions := groups[name]
		locants := make([]string, len(positions))
		for i, p := range positions {
			locants[i] = strconv.Itoa(p)
		}
		mult := ""
		if len(positions) >= 2 && len(positions) < len(multiplicityPrefixes) {
			mult = multiplicityPrefixes[len(positions)]
		}
		parts = append(parts, strings.Join(locants, ",")+"-"+mult+name)
	}
	prefix := strings.Join(parts, "-")
	name := prefix + parent // no extra hyphen before the parent name

	steps = append(steps,
		fmt.Sprintf("Assemble name: prefix = %s, parent = %s → %s.", prefix, parent, name),
	)
	return NamingResult{Name: name, Steps: steps}, nil
}

// ValidateName checks whether studentAnswer matches the correct IUPAC name.
// Comparison is case-insensitive and ignores spaces around hyphens.
func ValidateName(studentAnswer, correct string) bool {
	return normalize(studentAnswer) == normalize(correct)
}

func normalize(s string) string {
	s = strings.ToLower(s)
	s = strings.ReplaceAll(s, " - ", "-")
	s = strings.ReplaceAll(s, " -", "-")
	s = strings.ReplaceAll(s, "- ", "-")
	return strings.TrimSpace(s)
}
