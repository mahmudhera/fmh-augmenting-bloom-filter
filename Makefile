# Rebuild the theory report.
#   make        -> regenerate tables/figures (needs numpy + matplotlib) and the PDF
#   make pdf    -> rebuild the PDF only from the committed tables/figures (needs groff + ghostscript)
PYTHON ?= python3

all: tables pdf

report/gen/mc.csv: scripts/mc_check_theorem1.py
	$(PYTHON) scripts/mc_check_theorem1.py > $@

tables: report/gen/mc.csv scripts/compute_tables.py
	$(PYTHON) scripts/compute_tables.py report/gen report/figs report/gen/mc.csv

pdf:
	report/build.sh

.PHONY: all tables pdf
