#!/bin/bash
# source this script to get an environment with a python 3
# where sourcing importing from sage works
HELPER_DIR=$(dirname -- "${BASH_SOURCE[0]}")
SAGE_BIN_DIR=$(sage -sh -c 'echo $SAGE_LOCAL')/bin
if [ ! -f ${SAGE_BIN_DIR}/sage-env-config ]; then
    SAGE_BIN_DIR=$(sage -sh -c 'echo $SAGE_VENV')/bin
fi
source ${SAGE_BIN_DIR}/sage-env-config
if [ ! -f ${SAGE_BIN_DIR}/sage-env ]; then
    SAGE_BIN_DIR=$(sage -sh -c 'echo $SAGE_VENV')/bin
fi
source ${SAGE_BIN_DIR}/sage-env
export PATH=$(dirname -- ${HELPER_DIR})/venv/bin:$SAGE_BIN_DIR:${PATH}
