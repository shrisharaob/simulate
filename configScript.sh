#!/bin/sh
# Shell script to be run for setup on a new pc to create the file structure for storing simulation results, supply src code dir full path
set -x
#SRC_CODE_DIR="$HOME/Documents/code"
SRC_CODE_DIR=$1
DIRECTORY="/Documents/cnrs/simResults/"
FILEBASE="$HOME$DIRECTORY"
echo "creating file tree : $FILEBASE"
mkdir -p $FILEBASE # -p creats the tree structure
echo "done"
echo "#define FILEBASE \"$FILEBASE\"" > "$SRC_CODE_DIR/config.h"
