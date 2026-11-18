CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -fopenmp
BIN      := motor comparar generar arnes

all: $(BIN)

motor: src/main.cpp src/filtros.cpp src/imagen.cpp
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

comparar: tests/comparar.cpp
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

generar: bench/generar.cpp src/imagen.cpp
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

arnes: bench/arnes.cpp src/filtros.cpp src/imagen.cpp
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

prueba_filtros: tests/filtros.cpp src/filtros.cpp src/imagen.cpp
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

pruebas: all prueba_filtros
	./prueba_filtros
	bash tests/pruebas.sh

clean:
	rm -f $(BIN) prueba_filtros

.PHONY: all pruebas clean
