package thermo

import "github.com/shopspring/decimal"

// Solvent holds the physical constants needed for colligative property
// calculations. Kb and Kf are in °C·kg/mol. BoilingPoint and FreezingPoint
// are in °C at standard pressure.
type Solvent struct {
	Name         string
	BoilingPoint decimal.Decimal // °C
	FreezingPoint decimal.Decimal // °C
	Kb           decimal.Decimal // ebullioscopic constant, °C·kg/mol
	Kf           decimal.Decimal // cryoscopic constant, °C·kg/mol
}

// Common solvents with their colligative property constants.
var (
	Water = Solvent{
		Name:          "Water",
		BoilingPoint:  decimal.NewFromFloat(100.0),
		FreezingPoint: decimal.NewFromFloat(0.0),
		Kb:            decimal.NewFromFloat(0.512),
		Kf:            decimal.NewFromFloat(1.86),
	}
	Benzene = Solvent{
		Name:          "Benzene",
		BoilingPoint:  decimal.NewFromFloat(80.1),
		FreezingPoint: decimal.NewFromFloat(5.5),
		Kb:            decimal.NewFromFloat(2.53),
		Kf:            decimal.NewFromFloat(5.12),
	}
	Ethanol = Solvent{
		Name:          "Ethanol",
		BoilingPoint:  decimal.NewFromFloat(78.4),
		FreezingPoint: decimal.NewFromFloat(-114.6),
		Kb:            decimal.NewFromFloat(1.22),
		Kf:            decimal.NewFromFloat(1.99),
	}
	Cyclohexane = Solvent{
		Name:          "Cyclohexane",
		BoilingPoint:  decimal.NewFromFloat(80.7),
		FreezingPoint: decimal.NewFromFloat(6.5),
		Kb:            decimal.NewFromFloat(2.79),
		Kf:            decimal.NewFromFloat(20.2),
	}
	CarbonTetrachloride = Solvent{
		Name:          "Carbon Tetrachloride",
		BoilingPoint:  decimal.NewFromFloat(76.7),
		FreezingPoint: decimal.NewFromFloat(-22.9),
		Kb:            decimal.NewFromFloat(5.02),
		Kf:            decimal.NewFromFloat(29.8),
	}
	Chloroform = Solvent{
		Name:          "Chloroform",
		BoilingPoint:  decimal.NewFromFloat(61.2),
		FreezingPoint: decimal.NewFromFloat(-63.5),
		Kb:            decimal.NewFromFloat(3.63),
		Kf:            decimal.NewFromFloat(4.68),
	}
)
