SHELL := helper_scripts/bash_sage_shell

bin = venv/bin
env = env PATH="${bin}:$$PATH"
sage_python = python3
pysrcdirs = sage_code/ tests/ isogeny_primes.py latex_helper.py plot_stats.py
sage_version = >=9.4
version_command = "from sagemath.check_version import check_version;\
                   check_version(\"${sage_version}\")"

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
pip-install: venv/make_pip_install_complete
venv/make_pip_install_complete: requirements.txt
	. venv/bin/activate && ${env} pip install -Ur requirements.txt
	. venv/bin/activate && ${env} python -c ${version_command}
	touch venv/make_pip_install_complete

.PHONY: pip-install-dev
pip-install-dev: venv/make_pip_install_dev_complete
venv/make_pip_install_dev_complete: venv/make_pip_install_complete requirements-dev.txt
	. venv/bin/activate && ${env} pip install -Ur requirements-dev.txt
	touch venv/make_pip_install_dev_complete


.PHONY: unittests
unittests: pip-install-dev ## Run unittests using pytest
	. venv/bin/activate && ${env} coverage run -m pytest -vv --log-cli-level=DEBUG tests/fast_tests

.PHONY: integrationtests
integrationtests: pip-install-dev ## Run integrationtests using pytest
	. venv/bin/activate && ${env} coverage run -m pytest -vv --log-cli-level=DEBUG tests/slow_tests

.PHONY: test
test: pip-install-dev ## Run all tests using pytest
	. venv/bin/activate && ${env} coverage run -m pytest -vv -log-cli-level=DEBUG

.PHONY: test-report
test-report: pip-install-dev
	. venv/bin/activate && ${env} coverage report

.PHONY: black
black: pip-install-dev ## Check for source issues
	# verify that all pedantic source issues are resolved.
	@. venv/bin/activate && ${env} python3 -m black --check ${pysrcdirs}

.PHONY: check-types
check-types: pip-install-dev ## Check for type issues with mypy
	@. venv/bin/activate && ${env} python3 -m mypy --check ${pysrcdirs}

.PHONY: isort
isort: pip-install-dev
	@. venv/bin/activate && ${env} python3 -m isort ${pysrcdirs}

.PHONY: fix
fix: pip-install-dev ## Automatically fix style issues
	# @. .venv/bin/activate && ${env} python3 -m isort ${pysrcdirs}
	# format using black
	@. venv/bin/activate && ${env} python3 -m black ${pysrcdirs}

	# autoflake removes unused imports and unused variables from Python code. It makes use of pyflakes to do this.
	@. venv/bin/activate && ${env} python3 -m autoflake -ri --remove-all-unused-imports ${pysrcdirs}
	${MAKE} black

.PHONY: vulture
vulture: pip-install-dev
	# searching for unreachable code
	@. venv/bin/activate && ${env} python3 -m vulture ${pysrcdirs} --min-confidence 100

.PHONY: lint
lint: pip-install-dev  ## Do basic linting
	@. venv/bin/activate && ${env} pylint --extension-pkg-allow-list=sage ${pysrcdirs}

# should add lint at some point but still has to many failures at the moment
.PHONY: valid
valid: pip-install-dev vulture fix test test-report

.PHONY: valid_fast
valid_fast: pip-install-dev vulture fix unittests test-report


clean: ## Cleanup
clean: clean_venv clean_pytest clean_coverage clean_pycache

clean_venv:
	@echo "Removing venv"
	@rm -rf venv

clean_pytest:
	@echo "Removing pytest_cache"
	@rm -rf .pytest_cache

clean_coverage:
	@echo "Removing .coverage"
	@rm -f .coverage

clean_pycache:
	@echo "Removing all pyc and pyo files and all __pycache__ direcrories"
	@find . -type f -name '*.py[co]' -delete -o -type d -name __pycache__ -delete
