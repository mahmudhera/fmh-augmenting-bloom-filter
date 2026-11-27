# Rebuild the theory report and the aif tool.
#   make aif    -> build bin/aif (C++17; zlib for gzip input)
#   make pdf    -> rebuild the theory PDF from committed tables/figures (groff + ghostscript)
#   make        -> regenerate theory tables/figures (numpy + matplotlib) and the PDF
PYTHON ?= python3
CXX ?= g++
CXXFLAGS ?= -O3 -march=native -mprefer-vector-width=512 -std=c++17 -pthread -fno-math-errno -fassociative-math -fno-signed-zeros -fno-trapping-math -freciprocal-math

all: tables pdf

aif: bin/aif
bin/aif: src/aif.cpp src/ci_table.inc src/ci_table_db.inc
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -Wall -Wextra -o $@ $< -lz

report/gen/mc.csv: scripts/mc_check_theorem1.py
	$(PYTHON) scripts/mc_check_theorem1.py > $@

tables: report/gen/mc.csv scripts/compute_tables.py
	$(PYTHON) scripts/compute_tables.py report/gen report/figs report/gen/mc.csv

pdf:
	report/build.sh

.PHONY: all tables pdf aif

# per-stage cycle counters (stderr) for profiling
prof: bin/aif_prof
bin/aif_prof: src/aif.cpp src/ci_table.inc src/ci_table_db.inc
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -DAIF_PROF -Wall -Wextra -o $@ $< -lz
