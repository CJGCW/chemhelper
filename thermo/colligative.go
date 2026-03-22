// Package thermo implements thermodynamic and colligative property
// calculations. Each calculation type satisfies the calc.Calculation interface.
package thermo

import (
	"fmt"

	"chemhelper/calc"

	"github.com/shopspring/decimal"
)

var one = decimal.NewFromInt(1)

// BoilingPointElevation calculates the elevation of a solvent's boiling point
// caused by a dissolved solute.
//
//	ΔTb = Kb × m × i
//
// Molality is in mol/kg. VantHoffFactor defaults to 1 if not set (ideal
// non-electrolyte). For electrolytes, i reflects the number of particles the
// solute dissociates into (e.g. NaCl → i=2, CaCl2 → i=3).
type BoilingPointElevation struct {
	Solvent        Solvent
	Molality       decimal.Decimal // mol/kg
	VantHoffFactor decimal.Decimal // dimensionless; defaults to 1
}

func (b BoilingPointElevation) Validate() error {
	if b.Molality.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("molality must be a positive value, got %v", b.Molality)
	}
	if b.Solvent.Kb.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("solvent Kb must be a positive value, got %v", b.Solvent.Kb)
	}
	if b.VantHoffFactor.LessThan(decimal.Zero) {
		return fmt.Errorf("van't Hoff factor must be a positive value, got %v", b.VantHoffFactor)
	}
	return nil
}

func (b BoilingPointElevation) vantHoff() decimal.Decimal {
	if b.VantHoffFactor.Equal(decimal.Zero) {
		return one
	}
	return b.VantHoffFactor
}

func (b BoilingPointElevation) Calculate() (calc.Result, error) {
	if err := b.Validate(); err != nil {
		return calc.Result{}, err
	}
	i := b.vantHoff()
	deltaTb := b.Solvent.Kb.Mul(b.Molality).Mul(i)
	newBp := b.Solvent.BoilingPoint.Add(deltaTb)
	return calc.Result{
		Value: deltaTb,
		Unit:  "°C",
		Steps: []string{
			fmt.Sprintf("Identify Kb for %s: %v °C·kg/mol", b.Solvent.Name, b.Solvent.Kb),
			fmt.Sprintf("Van't Hoff factor (i): %v", i),
			fmt.Sprintf("Apply ΔTb = Kb × m × i: %v × %v × %v = %v °C", b.Solvent.Kb, b.Molality, i, deltaTb),
			fmt.Sprintf("New boiling point: %v + %v = %v °C", b.Solvent.BoilingPoint, deltaTb, newBp),
		},
	}, nil
}

// NewBoilingPoint returns the elevated boiling point of the solution (°C)
// rather than just the delta. It calls Calculate internally.
func (b BoilingPointElevation) NewBoilingPoint() (calc.Result, error) {
	result, err := b.Calculate()
	if err != nil {
		return calc.Result{}, err
	}
	newBp := b.Solvent.BoilingPoint.Add(result.Value)
	result.Value = newBp
	result.Unit = "°C"
	result.Steps = append(result.Steps,
		fmt.Sprintf("New boiling point: %v °C", newBp),
	)
	return result, nil
}

// FreezingPointDepression calculates the depression of a solvent's freezing
// point caused by a dissolved solute.
//
//	ΔTf = Kf × m × i
//
// Molality is in mol/kg. VantHoffFactor defaults to 1 if not set.
type FreezingPointDepression struct {
	Solvent        Solvent
	Molality       decimal.Decimal // mol/kg
	VantHoffFactor decimal.Decimal // dimensionless; defaults to 1
}

func (f FreezingPointDepression) Validate() error {
	if f.Molality.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("molality must be a positive value, got %v", f.Molality)
	}
	if f.Solvent.Kf.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("solvent Kf must be a positive value, got %v", f.Solvent.Kf)
	}
	if f.VantHoffFactor.LessThan(decimal.Zero) {
		return fmt.Errorf("van't Hoff factor must be a positive value, got %v", f.VantHoffFactor)
	}
	return nil
}

func (f FreezingPointDepression) vantHoff() decimal.Decimal {
	if f.VantHoffFactor.Equal(decimal.Zero) {
		return one
	}
	return f.VantHoffFactor
}

func (f FreezingPointDepression) Calculate() (calc.Result, error) {
	if err := f.Validate(); err != nil {
		return calc.Result{}, err
	}
	i := f.vantHoff()
	deltaTf := f.Solvent.Kf.Mul(f.Molality).Mul(i)
	newFp := f.Solvent.FreezingPoint.Sub(deltaTf)
	return calc.Result{
		Value: deltaTf,
		Unit:  "°C",
		Steps: []string{
			fmt.Sprintf("Identify Kf for %s: %v °C·kg/mol", f.Solvent.Name, f.Solvent.Kf),
			fmt.Sprintf("Van't Hoff factor (i): %v", i),
			fmt.Sprintf("Apply ΔTf = Kf × m × i: %v × %v × %v = %v °C", f.Solvent.Kf, f.Molality, i, deltaTf),
			fmt.Sprintf("New freezing point: %v - %v = %v °C", f.Solvent.FreezingPoint, deltaTf, newFp),
		},
	}, nil
}

// NewFreezingPoint returns the depressed freezing point of the solution (°C)
// rather than just the delta. It calls Calculate internally.
func (f FreezingPointDepression) NewFreezingPoint() (calc.Result, error) {
	result, err := f.Calculate()
	if err != nil {
		return calc.Result{}, err
	}
	newFp := f.Solvent.FreezingPoint.Sub(result.Value)
	result.Value = newFp
	result.Unit = "°C"
	result.Steps = append(result.Steps,
		fmt.Sprintf("New freezing point: %v °C", newFp),
	)
	return result, nil
}

// FindMolalityFromBPE back-calculates molality from a measured boiling point
// elevation. Useful when you know ΔTb and want to find the concentration.
//
//	m = ΔTb / (Kb × i)
type FindMolalityFromBPE struct {
	Solvent        Solvent
	DeltaTb        decimal.Decimal // °C, must be positive
	VantHoffFactor decimal.Decimal // defaults to 1
}

func (f FindMolalityFromBPE) Validate() error {
	if f.DeltaTb.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("ΔTb must be a positive value, got %v", f.DeltaTb)
	}
	if f.Solvent.Kb.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("solvent Kb must be a positive value, got %v", f.Solvent.Kb)
	}
	if f.VantHoffFactor.LessThan(decimal.Zero) {
		return fmt.Errorf("van't Hoff factor must be a positive value, got %v", f.VantHoffFactor)
	}
	return nil
}

func (f FindMolalityFromBPE) vantHoff() decimal.Decimal {
	if f.VantHoffFactor.Equal(decimal.Zero) {
		return one
	}
	return f.VantHoffFactor
}

func (f FindMolalityFromBPE) Calculate() (calc.Result, error) {
	if err := f.Validate(); err != nil {
		return calc.Result{}, err
	}
	i := f.vantHoff()
	molality := f.DeltaTb.Div(f.Solvent.Kb.Mul(i))
	return calc.Result{
		Value: molality,
		Unit:  "mol/kg",
		Steps: []string{
			fmt.Sprintf("Identify Kb for %s: %v °C·kg/mol", f.Solvent.Name, f.Solvent.Kb),
			fmt.Sprintf("Van't Hoff factor (i): %v", i),
			fmt.Sprintf("Rearrange ΔTb = Kb × m × i → m = ΔTb / (Kb × i)"),
			fmt.Sprintf("Calculate molality: %v / (%v × %v) = %v mol/kg", f.DeltaTb, f.Solvent.Kb, i, molality),
		},
	}, nil
}

// FindMolalityFromFPD back-calculates molality from a measured freezing point
// depression.
//
//	m = ΔTf / (Kf × i)
type FindMolalityFromFPD struct {
	Solvent        Solvent
	DeltaTf        decimal.Decimal // °C, must be positive
	VantHoffFactor decimal.Decimal // defaults to 1
}

func (f FindMolalityFromFPD) Validate() error {
	if f.DeltaTf.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("ΔTf must be a positive value, got %v", f.DeltaTf)
	}
	if f.Solvent.Kf.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("solvent Kf must be a positive value, got %v", f.Solvent.Kf)
	}
	if f.VantHoffFactor.LessThan(decimal.Zero) {
		return fmt.Errorf("van't Hoff factor must be a positive value, got %v", f.VantHoffFactor)
	}
	return nil
}

func (f FindMolalityFromFPD) vantHoff() decimal.Decimal {
	if f.VantHoffFactor.Equal(decimal.Zero) {
		return one
	}
	return f.VantHoffFactor
}

func (f FindMolalityFromFPD) Calculate() (calc.Result, error) {
	if err := f.Validate(); err != nil {
		return calc.Result{}, err
	}
	i := f.vantHoff()
	molality := f.DeltaTf.Div(f.Solvent.Kf.Mul(i))
	return calc.Result{
		Value: molality,
		Unit:  "mol/kg",
		Steps: []string{
			fmt.Sprintf("Identify Kf for %s: %v °C·kg/mol", f.Solvent.Name, f.Solvent.Kf),
			fmt.Sprintf("Van't Hoff factor (i): %v", i),
			fmt.Sprintf("Rearrange ΔTf = Kf × m × i → m = ΔTf / (Kf × i)"),
			fmt.Sprintf("Calculate molality: %v / (%v × %v) = %v mol/kg", f.DeltaTf, f.Solvent.Kf, i, molality),
		},
	}, nil
}
