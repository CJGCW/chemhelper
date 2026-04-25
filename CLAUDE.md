# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
go build ./...          # build all packages
go test ./...           # run all tests
go test ./api/...       # run tests for a specific package
go test -run TestName   # run a single test by name
go test -v ./...        # verbose test output
```

No Makefile or task runner — standard Go tooling only.

## Architecture

This is a Go REST API (port 8080) providing chemistry calculations for a React frontend (expected at localhost:5173).

### Layer structure

```
main.go
  └── api/          HTTP handlers, routing (chi/v5), request/response models
        └── imports domain packages:
              calc/       Calculation interface and Result type
              units/      Mass, Volume, Prefix — type-safe unit conversion
              element/    Periodic table (118 elements), valence, properties
              smiles/     PubChem REST client with in-memory cache
              solution/   Molarity, molality, dilution calculations
              thermo/     Colligative properties, solvent database
              structure/  Lewis structure generator and static registry
```

### Key abstractions

**`calc.Calculation` interface** — every chemistry calculation implements:
- `Validate() error`
- `Calculate() (Result, error)` — returns Value, Unit, SigFigs, Steps

**`units.Property` interface** — polymorphic unit handling for Mass and Volume:
- `ConvertToStandard() (decimal.Decimal, error)`
- `GetMoles(decimal.Decimal) (decimal.Decimal, error)`

All calculations use `shopspring/decimal` for precision — avoid `float64` for chemistry values.

### API surface (api/router.go)

| Group | Prefix | Key endpoints |
|---|---|---|
| Compound | `/api/compound` | `/resolve`, `/lookup` (PubChem via SMILES) |
| Elements | `/api/elements` | `GET /`, `GET /{symbol}`, `POST /compound` |
| Solution | `/api/solution` | molarity, molality, dilution variants |
| Thermo | `/api/thermo` | BPE, FPD, back-calculate molality, solvents |
| Structure | `/api/structure` | `/lewis`, `/random` |

### Testing pattern

Tests live in `api/handlers_test.go` using `httptest` with black-box `package api_test` style. Domain packages each have their own `*_test.go`. Run a targeted package test with `go test ./structure/` etc.

### External dependency

`smiles/` calls the PubChem REST API and caches results in a `sync.RWMutex`-protected map. New compound lookups are slow (network); cached ones are instant.

### Structure registry

`structure/` maintains a hand-verified static registry of Lewis structures alongside an algorithmic generator. When adding a new molecule, prefer adding it to the static registry if the algorithm produces an incorrect result.
