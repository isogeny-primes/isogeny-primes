bin = venv/bin
env = source helper_scripts/bash_sage_init.sh && env PATH="${bin}:$$PATH"
sage_python = python3
pysrcdirs = sage_code/ tests/ isogeny_primes.py latex_helper.py plot_stats.py

# this script should automatically get the correct python from sage
# if not, you can try
# make venv sage_python=/path/to/sagemath/local/bin/python3
venv:
	${env} && ${sage_python} -m venv --system-site-packages venv
	. venv/bin/activate && ${env} pip install -U pip pip-tools

requirements.txt: venv requirements.in
	. venv/bin/activate && ${env} python3 -m piptools compile requirements.in

requirements-dev.txt: venv requirements-dev.in
	. venv/bin/activate && ${env} python3 -m piptools compile requirements-dev.in

.PHONY:	pip-compile
pip-compile: requirements.txt requirements-dev.txt

.PHONY: pip-install
pip-install: requirements.txt
	. venv/bin/activate && ${env} pip install -Ur requirements.txt

.PHONY: pip-install-dev
pip-install-dev: pip-install requirements-dev.txt
	. venv/bin/activate && ${env} pip install -Ur requirements-dev.txt

.PHONY: unittests
unittests: venv ## Run unittests using pytest
    # Runs all testcases and delivers a coverage report to your terminal
	. venv/bin/activate && ${env} coverage run -m pytest -vv --log-cli-level=DEBUG tests/unit_tests

.PHONY: integrationtests
integrationtests: venv ## Run integrationtests using pytest
	. venv/bin/activate && ${env} coverage run -m pytest -vv --log-cli-level=DEBUG tests/integration_tests

.PHONY: test
test: venv ## Run all tests using pytest
	. venv/bin/activate && ${env} coverage run -m pytest -vv -log-cli-level=DEBUG

.PHONY: test-report
test-report: venv
	. venv/bin/activate && ${env} coverage report

.PHONY: black
black: venv ## Check for source issues
	# verify that all pedantic source issues are resolved.
	@. venv/bin/activate && ${env} python3 -m black --check ${pysrcdirs}

.PHONY: check-types
check-types: venv ## Check for type issues with mypy
	@. venv/bin/activate && ${env} python3 -m mypy --check ${pysrcdirs}

.PHONY: isort
isort: venv
	@. venv/bin/activate && ${env} python3 -m isort ${pysrcdirs}

.PHONY: fix
fix: venv ## Automatically fix style issues
	# @. .venv/bin/activate && ${env} python3 -m isort ${pysrcdirs}
	@. venv/bin/activate && ${env} python3 -m black ${pysrcdirs}

	# autoflake removes unused imports and unused variables from Python code. It makes use of pyflakes to do this.
	@. venv/bin/activate && ${env} python3 -m autoflake -ri --remove-all-unused-imports ${pysrcdirs}
	${MAKE} black

.PHONY: vulture
vulture: venv
	@. venv/bin/activate && ${env} python3 -m vulture ${pysrcdirs} --min-confidence 100

.PHONY: lint
lint: venv  ## Do basic linting
	@. venv/bin/activate && ${env} pylint --extension-pkg-allow-list=sage ${pysrcdirs}

# should add lint at some point but still has to many failures at the moment
.PHONY: valid
valid: venv vulture fix test test-report

.PHONY: valid_fast
valid_fast: venv vulture fix unittests test-report


clean: ## Cleanup
clean: clean_venv

clean_venv:  # Remove venv
	@echo "Cleaning venv"
	@rm -rf venv

