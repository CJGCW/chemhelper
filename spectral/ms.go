package spectral

// FragmentLoss describes a common neutral loss in mass spectrometry.
type FragmentLoss struct {
	Mass        int
	Formula     string
	DiagnosticFor string
}

// FragmentTable contains common neutral losses sorted by mass.
// Useful for interpreting the difference between M⁺ and fragment ions.
var FragmentTable = []FragmentLoss{
	{Mass: 1, Formula: "H", DiagnosticFor: ""},
	{Mass: 15, Formula: "CH₃", DiagnosticFor: "methyl group"},
	{Mass: 17, Formula: "OH", DiagnosticFor: "alcohol (less common)"},
	{Mass: 18, Formula: "H₂O", DiagnosticFor: "alcohol, carboxylic acid"},
	{Mass: 27, Formula: "HCN", DiagnosticFor: "aromatic nitrile, pyridine"},
	{Mass: 28, Formula: "CO / C₂H₄", DiagnosticFor: "aldehyde, ketone / ethyl-bearing compound"},
	{Mass: 29, Formula: "CHO / C₂H₅", DiagnosticFor: "aldehyde / ethyl group"},
	{Mass: 31, Formula: "OCH₃", DiagnosticFor: "methyl ester, methyl ether"},
	{Mass: 32, Formula: "CH₃OH", DiagnosticFor: "methyl ester"},
	{Mass: 43, Formula: "COCH₃ / C₃H₇", DiagnosticFor: "methyl ketone / propyl group"},
	{Mass: 44, Formula: "CO₂ / CH₂=CHOH", DiagnosticFor: "carboxylic acid / enol"},
	{Mass: 45, Formula: "OC₂H₅", DiagnosticFor: "ethyl ester"},
	{Mass: 57, Formula: "C₄H₉", DiagnosticFor: "tert-butyl or n-butyl group"},
	{Mass: 77, Formula: "C₆H₅", DiagnosticFor: "monosubstituted benzene"},
	{Mass: 91, Formula: "C₇H₇⁺", DiagnosticFor: "tropylium / benzyl cation"},
}

// PredictMS returns an empty stub.
func PredictMS(_ string) (SpectrumPrediction, error) {
	return SpectrumPrediction{Type: MassSpec}, nil
}
