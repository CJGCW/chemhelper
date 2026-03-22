package element

import (
	"errors"
	"fmt"
	"regexp"
	"strings"

	"github.com/shopspring/decimal"
)

func SetToSigFigs(val decimal.Decimal, sigfig int32) (decimal.Decimal, error) {
	if sigfig < 1 {
		return decimal.Zero, errors.New("significant figures must be greater than 0")
	}
	digits := val.NumDigits()
	scale := val.Exponent()
	d := int32(sigfig) - (int32(digits) + scale)
	return val.RoundBank(d), nil
}

//will need to be a string, Go doesn't handle floats as we want, i.e. 10.0 is treated as 10
func GetSignificantFigures(numStr string) (int, error) {
	re := regexp.MustCompile(`\d+\.*\d*`)
	matches := re.FindStringSubmatch(numStr)
	numStr = re.FindString(numStr)
	if matches == nil {
		return 0, fmt.Errorf("mass was not in a recognized format")
	}
	numStr = strings.TrimSpace(numStr)
	if strings.Contains(numStr, ".") {
		numStr = strings.Replace(numStr,".","",1)
		numStr = strings.TrimLeft(numStr, "0")
		return len(numStr), nil
	} else {
		numStr = strings.TrimLeft(strings.TrimRight(numStr, "0"), "0")
		return len(numStr), nil
	}
}

func GetLowestSignificantFigures(nums []string) (int, error) {
	if len(nums) == 0 {
		return 0, fmt.Errorf("no masses were passed")
	}
	var lowestNum int = -1;
	for _, num := range nums {
		newNum, err := GetSignificantFigures(num)
		if err != nil {
			return 0, fmt.Errorf("one or more masses was not in the correct format")
		}
		if newNum < lowestNum || lowestNum < 0{
			lowestNum = newNum
		}
	}
	return lowestNum, nil
}