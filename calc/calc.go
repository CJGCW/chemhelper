package calc

import "github.com/shopspring/decimal"

// Result holds the output of any calculation.
type Result struct {
	Value   decimal.Decimal
	Unit    string
	SigFigs int
	// Steps optionally records a worked solution, useful for a teaching tool.
	Steps []string
}

// Calculation is implemented by any struct that represents a chemistry
// calculation. Validate checks that inputs are well-formed before Calculate
// is called.
type Calculation interface {
	Validate() error
	Calculate() (Result, error)
}
