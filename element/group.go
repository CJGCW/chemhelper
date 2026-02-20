package element

import (
	"fmt"
	"slices"
)

var metalloids = []int{5,14,32,33,51,52,84,85}
var chalcogens = []int{8,16,34}
var halogens = []int{9,17,35,53}
type Group int
const (
    Alkali Group = iota +1
    Alkaline
    Metal 
    Pnictogens = 15
    Chalcogens = 16
    Halogens = 17
    NobleGases = 18

)
func (e Element) GetGroup() (string, error) {
    //Atomic number specifics
     //Then generalize by group
     g := Group(e.Group)
    switch {
    case e.AtomicNumber == 1 :
        return "None", nil
    case (e.AtomicNumber==6):
        return "Carbon", nil
    case (e.AtomicNumber==7 || e.AtomicNumber==15):
        return "Pnictogens", nil
    case slices.Contains(chalcogens,e.AtomicNumber) :
        return "Chalcogens", nil
    case slices.Contains(halogens,e.AtomicNumber) :
        return "Halogens", nil
    case slices.Contains(metalloids,e.AtomicNumber) :
        return "Metalloid", nil
    case g == Alkali:
        return "Alkali Metals", nil
    case g == Alkaline:
        return "Alkaline Earth Metals", nil
    case g == NobleGases:
        return "Noble Gases", nil
    case g < 18:
        return "Metals", nil //It's an all encompassing Group 3-17 all contain metals or metalloids at some point
    default:
        return "Unknown", fmt.Errorf("unknown element passed")
    }   
}