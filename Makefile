bin = venv/bin
env = env PATH="${bin}:$$PATH"
sage_python = python3
pysrcdirs = sage_code/
blackdirs =


# make sure that the sage system sage python is in the path under python3
# otherwise run
# make venv sage_python=python3
venv:
	${sage_python} -m venv --system-site-packages venv
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


test: venv ## Run unittests
    # Runs all testcases and delivers a coverage report to your terminal
	. venv/bin/activate && ${env} coverage run -m pytest -vv

test-report: venv
	. venv/bin/activate && ${env} coverage report

.PHONY: black
black: venv ## Check for source issues
	# verify that all pedantic source issues are resolved.
	@. venv/bin/activate && ${env} python3 -m black --check ${pysrcdirs} ${blackdirs}

check-types: venv ## Check for type issues with mypy
	@. venv/bin/activate && ${env} python3 -m mypy --check ${pysrcdirs}

isort: venv
	@. venv/bin/activate && ${env} python3 -m isort ${pysrcdirs}

fix: venv ## Automatically fix style issues
	# @. .venv/bin/activate && ${env} python3 -m isort ${pysrcdirs}

	@. venv/bin/activate && ${env} python3 -m black ${pysrcdirs} ${blackdirs}

	# autoflake removes unused imports and unused variables from Python code. It makes use of pyflakes to do this.
	@. venv/bin/activate && ${env} python3 -m autoflake -ri --remove-all-unused-imports ${pysrcdirs} ${blackdirs}
	${MAKE} check

vulture: venv
	@. venv/bin/activate && ${env} python3 -m vulture ${pysrcdirs} --min-confidence 100

lint: venv  ## Do basic linting
	@. venv/bin/activate && ${env} pylint ${pysrcdirs}


# isort and black linting have different ideas on correctness. isort is cleaner, and most of that is kept by black.
.PHONY: valid
valid: venv vulture isort fix lint check-types audit test test-report



clean: ## Cleanup
clean: clean_venv

clean_venv:  # Remove venv
	@echo "Cleaning venv"
	@rm -rf venv

