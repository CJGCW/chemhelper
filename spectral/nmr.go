package spectral

// HNMRShift maps a proton environment to its typical ¹H chemical shift range.
// Source: Brown & Foote, Organic Chemistry, Table 13.3.
type HNMRShift struct {
	ProtonType string
	MinPPM     float64
	MaxPPM     float64
	Example    string
}

// NMRHTable is the reference table for ¹H NMR chemical shifts.
var NMRHTable = []HNMRShift{
	{ProtonType: "R–CH₃", MinPPM: 0.8, MaxPPM: 1.0, Example: "ethane"},
	{ProtonType: "R₂CH₂", MinPPM: 1.2, MaxPPM: 1.4, Example: "propane CH₂"},
	{ProtonType: "R₃CH", MinPPM: 1.4, MaxPPM: 1.7, Example: "cyclopentane"},
	{ProtonType: "Allylic –CH– (C=C–CH)", MinPPM: 1.6, MaxPPM: 1.9, Example: "propene CH₃"},
	{ProtonType: "C≡C–H", MinPPM: 1.8, MaxPPM: 2.5, Example: "propyne"},
	{ProtonType: "CO–CH (α to carbonyl)", MinPPM: 2.0, MaxPPM: 2.5, Example: "acetone CH₃"},
	{ProtonType: "X–CH (X = halogen)", MinPPM: 2.5, MaxPPM: 4.0, Example: "CH₃Cl"},
	{ProtonType: "RO–CH (ether)", MinPPM: 3.3, MaxPPM: 3.9, Example: "diethyl ether"},
	{ProtonType: "=CH– (vinyl)", MinPPM: 4.6, MaxPPM: 5.3, Example: "ethylene"},
	{ProtonType: "Ar–H", MinPPM: 6.5, MaxPPM: 8.5, Example: "benzene"},
	{ProtonType: "R–CHO (aldehyde)", MinPPM: 9.5, MaxPPM: 10.0, Example: "acetaldehyde"},
	{ProtonType: "RCOOH", MinPPM: 10.0, MaxPPM: 12.0, Example: "acetic acid"},
	{ProtonType: "R–OH (variable)", MinPPM: 1.0, MaxPPM: 5.0, Example: "ethanol; depends on conc."},
	{ProtonType: "Ar–OH (phenol)", MinPPM: 4.0, MaxPPM: 12.0, Example: "phenol"},
	{ProtonType: "R–NH (amine)", MinPPM: 0.5, MaxPPM: 5.0, Example: "variable; often broad"},
}

// CNMRShift maps a carbon environment to its typical ¹³C chemical shift range.
type CNMRShift struct {
	CarbonType string
	MinPPM     float64
	MaxPPM     float64
}

// NMRCTable is the reference table for ¹³C NMR chemical shifts.
var NMRCTable = []CNMRShift{
	{CarbonType: "R–CH₃ / R₂CH₂ / R₃CH (alkyl)", MinPPM: 0, MaxPPM: 50},
	{CarbonType: "C–O (ether, ester, alcohol)", MinPPM: 50, MaxPPM: 90},
	{CarbonType: "C=C (alkene)", MinPPM: 100, MaxPPM: 150},
	{CarbonType: "Aromatic C", MinPPM: 110, MaxPPM: 160},
	{CarbonType: "C=O (aldehyde / ketone)", MinPPM: 190, MaxPPM: 220},
	{CarbonType: "C=O (acid / ester / amide)", MinPPM: 160, MaxPPM: 185},
}

// PredictHNMR returns an empty stub.
func PredictHNMR(_ string) (NMRPrediction, error) {
	return NMRPrediction{Type: HNMR}, nil
}

// PredictCNMR returns an empty stub.
func PredictCNMR(_ string) (NMRPrediction, error) {
	return NMRPrediction{Type: CNMR}, nil
}
