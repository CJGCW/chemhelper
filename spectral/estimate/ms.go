package estimate

import "sort"

// estimateMS returns estimated mass spectrum peaks.
// Includes M+, isotope peaks (¹³C, Br/Cl), and diagnostic fragment losses.
func estimateMS(m *molecule) []Peak {
	mPlus := m.nominalMass()
	counts := m.elemCounts()

	var peaks []Peak

	// M+ (molecular ion)
	peaks = append(peaks, Peak{
		X: float64(mPlus), Y: 100,
		Label: "M⁺ (molecular ion)",
		Width: 0.3,
	})

	// M+1: ¹³C isotope (~1.1% per carbon)
	nC := counts["C"]
	if nC > 0 {
		mPlusOneIntensity := float64(nC) * 1.1
		peaks = append(peaks, Peak{
			X: float64(mPlus + 1), Y: mPlusOneIntensity,
			Label: "M+1 (¹³C isotope)",
			Width: 0.3,
		})
	}

	// Halogen isotope patterns
	nBr := counts["Br"]
	nCl := counts["Cl"]

	if nBr > 0 {
		// ⁷⁹Br/⁸¹Br: M+2 at ~100% relative to M+
		peaks = append(peaks, Peak{
			X: float64(mPlus + 2), Y: 97,
			Label: "M+2 (⁸¹Br isotope)",
			Width: 0.3,
		})
		if nBr >= 2 {
			peaks = append(peaks, Peak{
				X: float64(mPlus + 4), Y: 25,
				Label: "M+4 (di-Br)",
				Width: 0.3,
			})
		}
	}

	if nCl > 0 {
		// ³⁵Cl/³⁷Cl: M+2 at ~33%
		peaks = append(peaks, Peak{
			X: float64(mPlus + 2), Y: 33,
			Label: "M+2 (³⁷Cl isotope)",
			Width: 0.3,
		})
		if nCl >= 2 {
			peaks = append(peaks, Peak{
				X: float64(mPlus + 4), Y: 11,
				Label: "M+4 (di-Cl)",
				Width: 0.3,
			})
		}
	}

	// Diagnostic fragments
	frags := diagnosticFragments(m, mPlus)
	peaks = append(peaks, frags...)

	// Normalize: set base peak (highest Y) to 100
	maxY := 0.0
	for _, p := range peaks {
		if p.Y > maxY {
			maxY = p.Y
		}
	}
	if maxY > 0 {
		for i := range peaks {
			peaks[i].Y = peaks[i].Y / maxY * 100
		}
	}

	// Sort by m/z ascending
	sort.Slice(peaks, func(i, j int) bool {
		return peaks[i].X < peaks[j].X
	})

	return peaks
}

// diagnosticFragments returns characteristic fragment peaks based on functional groups.
func diagnosticFragments(m *molecule, mPlus int) []Peak {
	var frags []Peak

	// Loss of CH₃ (M-15): very common for branched aliphatics, methyls next to heteroatoms
	if mPlus-15 > 0 {
		frags = append(frags, Peak{
			X: float64(mPlus - 15), Y: 35,
			Label: "[M-CH₃]⁺", Width: 0.3,
		})
	}

	// Loss of OH (M-17): alcohols, carboxylic acids
	if hasAlcohol(m) || hasCarboxylicAcid(m) {
		if mPlus-17 > 0 {
			frags = append(frags, Peak{
				X: float64(mPlus - 17), Y: 50,
				Label: "[M-OH]⁺", Width: 0.3,
			})
		}
	}

	// Loss of H₂O (M-18): alcohols, carboxylic acids
	if hasAlcohol(m) || hasCarboxylicAcid(m) {
		if mPlus-18 > 0 {
			frags = append(frags, Peak{
				X: float64(mPlus - 18), Y: 60,
				Label: "[M-H₂O]⁺", Width: 0.3,
			})
		}
	}

	// Loss of CO (M-28): aldehydes, ketones
	if hasAldehyde(m) || hasKetone(m) {
		if mPlus-28 > 0 {
			frags = append(frags, Peak{
				X: float64(mPlus - 28), Y: 45,
				Label: "[M-CO]⁺", Width: 0.3,
			})
		}
	}

	// Loss of CHO (M-29): aldehydes (α-cleavage)
	if hasAldehyde(m) {
		if mPlus-29 > 0 {
			frags = append(frags, Peak{
				X: float64(mPlus - 29), Y: 70,
				Label: "[M-CHO]⁺ (α-cleavage)", Width: 0.3,
			})
		}
	}

	// Acylium ion from ketone α-cleavage: emit largest acyl fragment
	if hasKetone(m) {
		acylMass := largestAcylMass(m)
		if acylMass > 0 && acylMass < mPlus {
			frags = append(frags, Peak{
				X: float64(acylMass), Y: 80,
				Label: "acylium [RCO]⁺ (α-cleavage)", Width: 0.3,
			})
		}
	}

	// Loss of OR (M-OR): esters — loss of alkoxy group
	if hasEster(m) {
		alkOxMass := smallestAlkoxyMass(m)
		if alkOxMass > 0 && mPlus-alkOxMass > 0 {
			frags = append(frags, Peak{
				X: float64(mPlus - alkOxMass), Y: 65,
				Label: "[M-OR]⁺ (ester)", Width: 0.3,
			})
		}
	}

	// Loss of HCl (M-36/38): alkyl chlorides
	counts := m.elemCounts()
	if counts["Cl"] > 0 {
		if mPlus-36 > 0 {
			frags = append(frags, Peak{
				X: float64(mPlus - 36), Y: 40,
				Label: "[M-HCl]⁺", Width: 0.3,
			})
		}
	}

	// Loss of HBr (M-80): alkyl bromides
	if counts["Br"] > 0 {
		if mPlus-80 > 0 {
			frags = append(frags, Peak{
				X: float64(mPlus - 80), Y: 40,
				Label: "[M-HBr]⁺", Width: 0.3,
			})
		}
	}

	// m/z 77: phenyl cation (C₆H₅⁺) for aromatic compounds
	if hasAromaticRing(m) && mPlus > 77 {
		frags = append(frags, Peak{
			X: 77, Y: 55,
			Label: "C₆H₅⁺ (phenyl)", Width: 0.3,
		})
	}

	// m/z 91: tropylium (C₇H₇⁺) for benzyl/toluene systems
	counts2 := m.elemCounts()
	if hasAromaticRing(m) && counts2["C"] >= 7 && mPlus > 91 {
		frags = append(frags, Peak{
			X: 91, Y: 75,
			Label: "C₇H₇⁺ (tropylium)", Width: 0.3,
		})
	}

	// m/z 43: acetyl (CH₃CO⁺) for methyl ketones and esters
	if hasKetone(m) || hasEster(m) {
		if mPlus > 43 {
			frags = append(frags, Peak{
				X: 43, Y: 60,
				Label: "CH₃CO⁺ (acetyl)", Width: 0.3,
			})
		}
	}

	return frags
}

// largestAcylMass estimates the mass of the largest acyl fragment (RCO+) from ketone α-cleavage.
func largestAcylMass(m *molecule) int {
	best := 0
	for i := range m.atoms {
		a := &m.atoms[i]
		if a.element != "C" {
			continue
		}
		// Find carbonyl C's
		for _, nb := range a.nbrs {
			if nb.order == 2 && m.atoms[nb.idx].element == "O" {
				// Estimate acyl mass by counting atoms on one side of the C=O
				mass := 12 + 16 // C=O itself
				for _, nb2 := range a.nbrs {
					if nb2.idx != nb.idx {
						mass += fragmentMass(m, nb2.idx, i)
						break // take the larger side (first one found)
					}
				}
				if mass > best {
					best = mass
				}
			}
		}
	}
	return best
}

// fragmentMass estimates the nominal mass of the fragment rooted at atomIdx,
// not going back through parentIdx.
func fragmentMass(m *molecule, atomIdx, parentIdx int) int {
	nm, ok := nominalMasses[m.atoms[atomIdx].element]
	if !ok {
		return 0
	}
	total := nm + m.atoms[atomIdx].hCount
	for _, nb := range m.atoms[atomIdx].nbrs {
		if nb.idx == parentIdx {
			continue
		}
		total += fragmentMass(m, nb.idx, atomIdx)
	}
	return total
}

// smallestAlkoxyMass returns the nominal mass of the smallest OR group in an ester.
func smallestAlkoxyMass(m *molecule) int {
	smallest := 0
	for i := range m.atoms {
		if !isEsterC(m, i) {
			continue
		}
		// Find the -OR oxygen
		for _, nb := range m.atoms[i].nbrs {
			nbAtom := &m.atoms[nb.idx]
			if nb.order == 1 && nbAtom.element == "O" && nbAtom.hCount == 0 {
				// OR = O + R fragment
				orMass := 16 // O
				for _, nb2 := range nbAtom.nbrs {
					if nb2.idx != i {
						orMass += fragmentMass(m, nb2.idx, nb.idx)
					}
				}
				if smallest == 0 || orMass < smallest {
					smallest = orMass
				}
			}
		}
	}
	return smallest
}
