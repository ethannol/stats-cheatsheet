# Run everything:        make all
# Run one topic:         make 04   (any two-digit prefix)
# Re-download data:      make refresh
PY := uv run python
SCRIPTS := $(sort $(wildcard analyses/*.py))

.PHONY: all data test refresh clean notebooks $(SCRIPTS)

all: data $(SCRIPTS)

data:
	$(PY) -m miningstats.data

$(SCRIPTS):
	@echo "\n=== $@ ==="
	$(PY) $@

# `make 03` runs analyses/03_*.py
%:
	$(PY) $(wildcard analyses/$@_*.py)

test:
	uv run pytest -q

# Convert the percent-format scripts to Jupyter notebooks under notebooks/
notebooks:
	mkdir -p notebooks
	for f in analyses/*.py; do uv run jupytext --to ipynb --output notebooks/$$(basename $${f%.py}).ipynb $$f; done

refresh:
	rm -rf data/*
	$(MAKE) data

clean:
	rm -rf figures/*
