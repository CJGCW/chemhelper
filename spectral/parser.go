package spectral

import (
	"fmt"
)

// ParseJCAMPDX parses a JCAMP-DX spectral file and returns extracted peaks.
// This is a stub — full JCAMP-DX parsing is not implemented.
func ParseJCAMPDX(data []byte) ([]Peak, error) {
	if len(data) == 0 {
		return nil, fmt.Errorf("empty JCAMP-DX data")
	}
	return nil, fmt.Errorf("JCAMP-DX parsing not yet implemented")
}

// ParseCSV parses a simple two-column (x,y) CSV and returns peaks.
// This is a stub — full CSV parsing is not implemented.
func ParseCSV(data []byte) ([]Peak, error) {
	if len(data) == 0 {
		return nil, fmt.Errorf("empty CSV data")
	}
	return nil, fmt.Errorf("CSV parsing not yet implemented")
}
