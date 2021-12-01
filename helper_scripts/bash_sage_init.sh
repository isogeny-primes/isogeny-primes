#!/bin/bash
# source this script to get an environment with a python 3
# where sourcing importing from sage works
HELPER_DIR=$(dirname -- "${BASH_SOURCE[0]}")
SAGE_BIN_DIR=$(sage -sh -c 'sage -root')/bin
source ${SAGE_BIN_DIR}/sage-env-config
source ${SAGE_BIN_DIR}/sage-env
export PATH=$(dirname -- ${HELPER_DIR})/venv/bin:${PATH}
