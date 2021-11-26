#!/bin/bash
# source this script to get an environment with a python 3
# where sourcing importing from sage works
SAGE_BIN_DIR=$(sage -sh -c 'sage -root')/bin
source ${SAGE_BIN_DIR}/sage-env-config
source ${SAGE_BIN_DIR}/sage-env
