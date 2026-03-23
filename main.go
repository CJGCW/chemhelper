package main

import (
	"fmt"
	"log"
	"net/http"

	"chemhelper/api"
)

func main() {
	router := api.NewRouter()
	addr := ":8080"
	fmt.Printf("chemhelper API listening on %s\n", addr)
	if err := http.ListenAndServe(addr, router); err != nil {
		log.Fatalf("server error: %v", err)
	}
}
