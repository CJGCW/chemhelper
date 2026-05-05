// Package organic provides typed constants for organic chemistry reaction classification.
package organic

// ReactionCategory classifies broad families of organic reactions.
type ReactionCategory string

const (
	SubstitutionCategory  ReactionCategory = "substitution"
	EliminationCategory   ReactionCategory = "elimination"
	AdditionCategory      ReactionCategory = "addition"
	RearrangementCategory ReactionCategory = "rearrangement"
	OxidationCategory     ReactionCategory = "oxidation"
	ReductionCategory     ReactionCategory = "reduction"
)

// ReactionType narrows a reaction to its specific mechanism.
type ReactionType string

const (
	SN1     ReactionType = "SN1"
	SN2     ReactionType = "SN2"
	E1      ReactionType = "E1"
	E2      ReactionType = "E2"
	E1cb    ReactionType = "E1cb"
	AdE     ReactionType = "AdE"     // electrophilic addition
	AdN     ReactionType = "AdN"     // nucleophilic addition
	EAS     ReactionType = "EAS"     // electrophilic aromatic substitution
	NAS     ReactionType = "NAS"     // nucleophilic aromatic substitution
	Radical ReactionType = "radical"
)

// Regiochemistry describes which carbon in an unsymmetrical reaction is favored.
type Regiochemistry string

const (
	Markovnikov     Regiochemistry = "Markovnikov"
	AntiMarkovnikov Regiochemistry = "anti-Markovnikov"
	Regiospecific   Regiochemistry = "regiospecific"
)

// Stereochemistry describes the stereo outcome of a reaction.
type Stereochemistry string

const (
	Inversion    Stereochemistry = "inversion"
	Retention    Stereochemistry = "retention"
	Racemization Stereochemistry = "racemization"
	Syn          Stereochemistry = "syn"
	Anti         Stereochemistry = "anti"
)
