#!/bin/bash

set -e

echo "Running cromwell $CROMWELL, womtool $WOMTOOL and image $TAG"

if [[ $# -eq 1 ]]; then
    caper run $1 --docker $TAG --cromwell $CROMWELL --womtool $WOMTOOL --options ./test/cromwell_options.json -m metadata.json
fi

if [[ $# -eq 2 ]]; then
    caper run $1 --docker $TAG --cromwell $CROMWELL --womtool $WOMTOOL -i $2 --options ./test/cromwell_options.json -m metadata.json
fi
