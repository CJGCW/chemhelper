package spectral

// IRCorrelation maps a functional group to its IR absorption characteristics.
// Source: Brown & Foote, Organic Chemistry, Table 12.3.
type IRCorrelation struct {
	Group      string  // functional group description
	MinWN      float64 // minimum wavenumber (cm⁻¹)
	MaxWN      float64 // maximum wavenumber (cm⁻¹)
	Intensity  string  // strong / medium / variable
	Shape      string  // broad / sharp / two bands
	Notes      string
}

// IRTable contains 21 diagnostic IR correlations from Brown Table 12.3.
var IRTable = []IRCorrelation{
	// O–H / N–H stretches
	{Group: "O–H (free)", MinWN: 3580, MaxWN: 3650, Intensity: "strong", Shape: "sharp", Notes: "Alcohol in dilute solution"},
	{Group: "O–H (H-bonded)", MinWN: 3200, MaxWN: 3550, Intensity: "strong", Shape: "broad", Notes: "Alcohol in neat liquid"},
	{Group: "O–H (acid)", MinWN: 2500, MaxWN: 3300, Intensity: "strong", Shape: "very broad", Notes: "Carboxylic acid; overlaps C–H region"},
	{Group: "N–H (1°)", MinWN: 3300, MaxWN: 3500, Intensity: "medium", Shape: "two bands", Notes: "Primary amine"},
	{Group: "N–H (2°)", MinWN: 3300, MaxWN: 3500, Intensity: "medium", Shape: "sharp", Notes: "Secondary amine; one band"},
	// C–H stretches
	{Group: "C–H (sp³)", MinWN: 2850, MaxWN: 2960, Intensity: "strong", Shape: "sharp", Notes: "Alkyl C–H"},
	{Group: "C–H (sp²)", MinWN: 3020, MaxWN: 3100, Intensity: "medium", Shape: "sharp", Notes: "Vinyl / aromatic C–H"},
	{Group: "C–H (sp, alkyne)", MinWN: 3260, MaxWN: 3330, Intensity: "strong", Shape: "sharp", Notes: "Terminal alkyne"},
	{Group: "C–H (aldehyde)", MinWN: 2700, MaxWN: 2850, Intensity: "medium", Shape: "two bands", Notes: "Weak; ~2720 and ~2820 cm⁻¹"},
	// Triple bonds
	{Group: "C≡C", MinWN: 2100, MaxWN: 2260, Intensity: "variable", Shape: "sharp", Notes: "Absent in symmetrical alkynes"},
	{Group: "C≡N", MinWN: 2200, MaxWN: 2260, Intensity: "strong", Shape: "sharp", Notes: "Nitrile"},
	// C=O stretches — most diagnostic region
	{Group: "C=O (ketone)", MinWN: 1705, MaxWN: 1720, Intensity: "strong", Shape: "sharp", Notes: ""},
	{Group: "C=O (aldehyde)", MinWN: 1720, MaxWN: 1740, Intensity: "strong", Shape: "sharp", Notes: "Also shows aldehyde C–H at 2700–2850"},
	{Group: "C=O (carboxylic acid)", MinWN: 1700, MaxWN: 1725, Intensity: "strong", Shape: "sharp", Notes: "Paired with broad O–H"},
	{Group: "C=O (ester)", MinWN: 1735, MaxWN: 1750, Intensity: "strong", Shape: "sharp", Notes: "Also strong C–O near 1200"},
	{Group: "C=O (amide)", MinWN: 1630, MaxWN: 1690, Intensity: "strong", Shape: "sharp", Notes: "Lower than typical C=O"},
	{Group: "C=O (acyl chloride)", MinWN: 1770, MaxWN: 1815, Intensity: "strong", Shape: "sharp", Notes: "Highest-frequency carbonyl"},
	{Group: "C=O (anhydride)", MinWN: 1800, MaxWN: 1850, Intensity: "strong", Shape: "two bands", Notes: "Two C=O bands ~1820 and ~1760"},
	// C=C stretches
	{Group: "C=C (alkene)", MinWN: 1620, MaxWN: 1680, Intensity: "variable", Shape: "sharp", Notes: "Weak or absent if symmetrical"},
	{Group: "C=C (aromatic)", MinWN: 1450, MaxWN: 1600, Intensity: "variable", Shape: "multiple", Notes: "Usually two bands"},
	// C–O and N=O
	{Group: "C–O", MinWN: 1000, MaxWN: 1260, Intensity: "strong", Shape: "sharp", Notes: "Ether / ester / alcohol"},
	{Group: "N=O (nitro)", MinWN: 1500, MaxWN: 1570, Intensity: "strong", Shape: "two bands", Notes: "Asymmetric stretch; symmetric at ~1350"},
}

// PredictIR returns an empty stub — full prediction requires SMILES parsing.
func PredictIR(_ string) (SpectrumPrediction, error) {
	return SpectrumPrediction{Type: IR}, nil
}
