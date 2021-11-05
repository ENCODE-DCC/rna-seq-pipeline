#!/bin/bash

set -e

# first argument is the input file that will be modified
# prefix contains the pipeline name
INPUT=$1
PREFIX=$2

cat $INPUT | jq ".+{\"${PREFIX}.docker\": \"${TAG}\"}" > tmp.json
cp tmp.json $INPUT && rm tmp.json
