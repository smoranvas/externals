#!/bin/bash

echo ""
echo ">> Setting Environment Variables for CERNLIB"
echo ""

export SOFT=/user/a/alaoui/software # external, ahmed dir

export CERNVER=2005 # hardcoded
export CERN=${SOFT}/cern # ahmed dir
export CERN_LEVEL=$CERNVER
export CERN_ROOT=${CERN}/${CERN_LEVEL}
export CERN_LIB=${CERN_ROOT}/lib
export CERN_BIN=${CERN_ROOT}/bin
export CERN_INC=${CERN_ROOT}/include
export CERNLIBDIR=${CERN_LIB}
export CERN_INCLUDEDIR=${CERN_INC}

if [ ! -d $CERN_LIB ]; then
  echo "CERN Error: $CERN_LIB Not Found"
  return -1
fi

if [ ! -d $CERN_BIN ]; then
  echo "CERN Error: $CERN_BIN Not Found"
  return -1
fi

if [ ! -d $CERN_INC ]; then
  echo "CERN Error: $CERN_INC Not Found"
  return -1
fi

echo "CERN_LEVEL           is set to ${CERN_LEVEL}"
echo "CERN_ROOT            is set to ${CERN_ROOT}"
echo "CERN_LIB             is set to ${CERN_LIB}"
echo "CERN_BIN             is set to ${CERN_BIN}"
echo "CERN_INC             is set to ${CERN_INC}"

if [ -z "${PATH}" ]; then
  PATH=${CERN_BIN}
else
  PATH=${CERN_BIN}:${PATH}
fi

if [ -z "$LIBRARY_PATH" ]; then
  export LIBRARY_PATH="${CERN_LIB}"
else
  export LIBRARY_PATH="${CERN_LIB}:$LIBRARY_PATH"
fi

if [ -z "${LD_LIBRARY_PATH}" ]; then
  LD_LIBRARY_PATH=${CERN_LIB}
else
  LD_LIBRARY_PATH=${CERN_LIB}:${LD_LIBRARY_PATH}
fi
