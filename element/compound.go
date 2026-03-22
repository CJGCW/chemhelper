package element

import (
	"fmt"
	"sort"
	"strconv"

	"github.com/shopspring/decimal"
)

type ElementMoles struct { // when creating compounds. I could just have moles be part of the element struct, but this is less confusing when balancing equations.
	Element Element
	Moles   decimal.Decimal
}

type Compound struct {
	Symbol    string
	Elements  []ElementMoles
	Mass      Mass
	Volume    Volume
	MolarMass decimal.Decimal
	Moles     decimal.Decimal
}

// Orders by symbol
func sortElementMoles(elements []ElementMoles) {
	sort.Slice(elements, func(i, j int) bool {
		return elements[i].Element.Symbol < elements[j].Element.Symbol
	})
}

// parseFragment recursively parses a compound fragment, returning a map of
// element symbol -> count. It handles nested parentheses such as Ca(OH)2 or
// Fe2(SO4)3.
func parseFragment(s string, pt *PeriodicTable) (map[string]int64, error) {
	counts := make(map[string]int64)
	i := 0
	for i < len(s) {
		switch {
		case s[i] == '(':
			// Find the matching closing paren
			depth := 1
			j := i + 1
			for j < len(s) && depth > 0 {
				switch s[j] {
				case '(':
					depth++
				case ')':
					depth--
				}
				j++
			}
			if depth != 0 {
				return nil, fmt.Errorf("unmatched '(' in compound")
			}
			// j now points one past the ')'
			inner := s[i+1 : j-1]
			// Read optional multiplier after ')'
			k := j
			for k < len(s) && s[k] >= '0' && s[k] <= '9' {
				k++
			}
			multiplier := int64(1)
			if k > j {
				var err error
				multiplier, err = strconv.ParseInt(s[j:k], 10, 64)
				if err != nil {
					return nil, err
				}
			}
			sub, err := parseFragment(inner, pt)
			if err != nil {
				return nil, err
			}
			for sym, cnt := range sub {
				counts[sym] += cnt * multiplier
			}
			i = k

		case s[i] >= 'A' && s[i] <= 'Z':
			// Read element symbol: one uppercase + optional lowercase
			j := i + 1
			for j < len(s) && s[j] >= 'a' && s[j] <= 'z' {
				j++
			}
			symbol := s[i:j]
			if _, found := pt.FindElementBySymbol(symbol); !found {
				return nil, fmt.Errorf("element %s not found in the periodic table", symbol)
			}
			// Read optional count
			k := j
			for k < len(s) && s[k] >= '0' && s[k] <= '9' {
				k++
			}
			count := int64(1)
			if k > j {
				var err error
				count, err = strconv.ParseInt(s[j:k], 10, 64)
				if err != nil {
					return nil, err
				}
			}
			counts[symbol] += count
			i = k

		default:
			return nil, fmt.Errorf("unexpected character %q in compound", s[i])
		}
	}
	return counts, nil
}

func ParseCompoundElements(compound string, pt *PeriodicTable) ([]ElementMoles, error) {
	if compound == "" {
		return nil, fmt.Errorf("no compound symbols passed")
	}
	counts, err := parseFragment(compound, pt)
	if err != nil {
		return nil, err
	}
	elements := make([]ElementMoles, 0, len(counts))
	for symbol, count := range counts {
		element, _ := pt.FindElementBySymbol(symbol) // existence already verified in parseFragment
		elements = append(elements, ElementMoles{
			Element: *element,
			Moles:   decimal.NewFromInt(count),
		})
	}
	return elements, nil
}
