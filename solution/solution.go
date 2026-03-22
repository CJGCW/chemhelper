// Package solution implements concentration and dilution calculations.
// Each calculation type satisfies the calc.Calculation interface.
package solution

import (
	"fmt"

	"chemhelper/calc"
	"chemhelper/units"

	"github.com/shopspring/decimal"
)

var (
	gramsPerKg = decimal.NewFromInt(1000)
)

// FindMolarity calculates molarity (mol/L) from moles of solute and volume
// of solution.
//
//	M = n / V
type FindMolarity struct {
	Moles  decimal.Decimal
	Volume units.Volume
}

func (f FindMolarity) Validate() error {
	if f.Moles.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("moles must be a positive value, got %v", f.Moles)
	}
	if _, err := f.Volume.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid volume: %w", err)
	}
	return nil
}

func (f FindMolarity) Calculate() (calc.Result, error) {
	if err := f.Validate(); err != nil {
		return calc.Result{}, err
	}
	volL, err := f.Volume.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	molarity := f.Moles.Div(volL)
	return calc.Result{
		Value: molarity,
		Unit:  "mol/L",
		Steps: []string{
			fmt.Sprintf("Convert volume to litres: %v L", volL),
			fmt.Sprintf("Divide moles by volume: %v mol ÷ %v L = %v mol/L", f.Moles, volL, molarity),
		},
	}, nil
}

// FindMolarityFromMass calculates molarity (mol/L) from the mass of a solute,
// its molar mass, and the volume of solution. Useful when moles haven't been
// computed yet.
//
//	n = mass / molarMass
//	M = n / V
type FindMolarityFromMass struct {
	Mass      units.Mass
	MolarMass decimal.Decimal // g/mol
	Volume    units.Volume
}

func (f FindMolarityFromMass) Validate() error {
	if f.MolarMass.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("molar mass must be a positive value, got %v", f.MolarMass)
	}
	if _, err := f.Mass.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid mass: %w", err)
	}
	if _, err := f.Volume.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid volume: %w", err)
	}
	return nil
}

func (f FindMolarityFromMass) Calculate() (calc.Result, error) {
	if err := f.Validate(); err != nil {
		return calc.Result{}, err
	}
	massG, err := f.Mass.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	volL, err := f.Volume.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	moles := massG.Div(f.MolarMass)
	molarity := moles.Div(volL)
	return calc.Result{
		Value: molarity,
		Unit:  "mol/L",
		Steps: []string{
			fmt.Sprintf("Convert mass to grams: %v g", massG),
			fmt.Sprintf("Convert volume to litres: %v L", volL),
			fmt.Sprintf("Calculate moles: %v g ÷ %v g/mol = %v mol", massG, f.MolarMass, moles),
			fmt.Sprintf("Calculate molarity: %v mol ÷ %v L = %v mol/L", moles, volL, molarity),
		},
	}, nil
}

// FindMolesFromMolarity calculates moles of solute from molarity and volume
// of solution.
//
//	n = M × V
type FindMolesFromMolarity struct {
	Molarity decimal.Decimal // mol/L
	Volume   units.Volume
}

func (f FindMolesFromMolarity) Validate() error {
	if f.Molarity.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("molarity must be a positive value, got %v", f.Molarity)
	}
	if _, err := f.Volume.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid volume: %w", err)
	}
	return nil
}

func (f FindMolesFromMolarity) Calculate() (calc.Result, error) {
	if err := f.Validate(); err != nil {
		return calc.Result{}, err
	}
	volL, err := f.Volume.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	moles := f.Molarity.Mul(volL)
	return calc.Result{
		Value: moles,
		Unit:  "mol",
		Steps: []string{
			fmt.Sprintf("Convert volume to litres: %v L", volL),
			fmt.Sprintf("Calculate moles: %v mol/L × %v L = %v mol", f.Molarity, volL, moles),
		},
	}, nil
}

// FindMolality calculates molality (mol/kg) from moles of solute and mass of
// solvent. Note: molality uses mass of solvent, not solution.
//
//	m = n / kg(solvent)
type FindMolality struct {
	Moles       decimal.Decimal
	SolventMass units.Mass
}

func (f FindMolality) Validate() error {
	if f.Moles.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("moles must be a positive value, got %v", f.Moles)
	}
	if _, err := f.SolventMass.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid solvent mass: %w", err)
	}
	return nil
}

func (f FindMolality) Calculate() (calc.Result, error) {
	if err := f.Validate(); err != nil {
		return calc.Result{}, err
	}
	massG, err := f.SolventMass.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	massKg := massG.Div(gramsPerKg)
	molality := f.Moles.Div(massKg)
	return calc.Result{
		Value: molality,
		Unit:  "mol/kg",
		Steps: []string{
			fmt.Sprintf("Convert solvent mass to grams: %v g", massG),
			fmt.Sprintf("Convert grams to kilograms: %v g ÷ 1000 = %v kg", massG, massKg),
			fmt.Sprintf("Calculate molality: %v mol ÷ %v kg = %v mol/kg", f.Moles, massKg, molality),
		},
	}, nil
}

// DilutionFindFinalConcentration solves for the final concentration when a
// solution is diluted to a known final volume.
//
//	C2 = (C1 × V1) / V2
type DilutionFindFinalConcentration struct {
	InitialConcentration decimal.Decimal // mol/L
	InitialVolume        units.Volume
	FinalVolume          units.Volume
}

func (d DilutionFindFinalConcentration) Validate() error {
	if d.InitialConcentration.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("initial concentration must be a positive value, got %v", d.InitialConcentration)
	}
	if _, err := d.InitialVolume.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid initial volume: %w", err)
	}
	if _, err := d.FinalVolume.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid final volume: %w", err)
	}
	return nil
}

func (d DilutionFindFinalConcentration) Calculate() (calc.Result, error) {
	if err := d.Validate(); err != nil {
		return calc.Result{}, err
	}
	v1, err := d.InitialVolume.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	v2, err := d.FinalVolume.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	if v2.LessThanOrEqual(v1) {
		return calc.Result{}, fmt.Errorf("final volume (%v L) must be greater than initial volume (%v L)", v2, v1)
	}
	c2 := d.InitialConcentration.Mul(v1).Div(v2)
	return calc.Result{
		Value: c2,
		Unit:  "mol/L",
		Steps: []string{
			fmt.Sprintf("Convert initial volume to litres: %v L", v1),
			fmt.Sprintf("Convert final volume to litres: %v L", v2),
			fmt.Sprintf("Apply C1V1 = C2V2: (%v mol/L × %v L) ÷ %v L = %v mol/L", d.InitialConcentration, v1, v2, c2),
		},
	}, nil
}

// DilutionFindFinalVolume solves for the final volume needed to reach a target
// concentration.
//
//	V2 = (C1 × V1) / C2
type DilutionFindFinalVolume struct {
	InitialConcentration decimal.Decimal // mol/L
	InitialVolume        units.Volume
	FinalConcentration   decimal.Decimal // mol/L
}

func (d DilutionFindFinalVolume) Validate() error {
	if d.InitialConcentration.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("initial concentration must be a positive value, got %v", d.InitialConcentration)
	}
	if d.FinalConcentration.LessThanOrEqual(decimal.Zero) {
		return fmt.Errorf("final concentration must be a positive value, got %v", d.FinalConcentration)
	}
	if d.FinalConcentration.GreaterThanOrEqual(d.InitialConcentration) {
		return fmt.Errorf("final concentration (%v mol/L) must be less than initial concentration (%v mol/L)", d.FinalConcentration, d.InitialConcentration)
	}
	if _, err := d.InitialVolume.ConvertToStandard(); err != nil {
		return fmt.Errorf("invalid initial volume: %w", err)
	}
	return nil
}

func (d DilutionFindFinalVolume) Calculate() (calc.Result, error) {
	if err := d.Validate(); err != nil {
		return calc.Result{}, err
	}
	v1, err := d.InitialVolume.ConvertToStandard()
	if err != nil {
		return calc.Result{}, err
	}
	v2 := d.InitialConcentration.Mul(v1).Div(d.FinalConcentration)
	return calc.Result{
		Value: v2,
		Unit:  "L",
		Steps: []string{
			fmt.Sprintf("Convert initial volume to litres: %v L", v1),
			fmt.Sprintf("Apply C1V1 = C2V2: (%v mol/L × %v L) ÷ %v mol/L = %v L", d.InitialConcentration, v1, d.FinalConcentration, v2),
		},
	}, nil
}
