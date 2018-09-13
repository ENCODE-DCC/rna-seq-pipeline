#!/bin/bash

# from Jin Lee's excellent atac-seq-pipeline in https://www.github.com/encode-dcc/atac-seq-pipeline
# with minor changes

set -e # exit on error

if [ $# -lt 3 ]; then
  echo "Usage: ./test.sh [WDL] [INPUT_JSON] [DOCKER_IMAGE]"
  echo "Make sure to have cromwell-34.jar in your \$PATH as an executable (chmod +x)."
  exit 1
fi

WDL=$1
INPUT=$2
DOCKER_IMAGE=$3

if [ -f "cromwell-34.jar" ]; then
    echo "cromwell-34.jar already available, skipping download."
else
    wget -N -c https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
fi

CROMWELL_JAR=cromwell-34.jar
BACKEND_CONF=backends/backend.conf
RESULT_PREFIX=$(basename ${INPUT} .json)
METADATA=${RESULT_PREFIX}.metadata.json # metadata
RESULT=${RESULT_PREFIX}.result.json # output

if [ $4 = "docker" ]; then
    # Write workflow option JSON file for docker
    BACKEND=Local
    TMP_WF_OPT=$RESULT_PREFIX.test_rna_wf_opt.json
    cat > $TMP_WF_OPT << EOM
    {
        "default_runtime_attributes" : {
            "docker" : "$DOCKER_IMAGE"
        }
    }
EOM
fi

if [ $4 = "singularity" ]; then
    #build singularity container
    SINGULARITY_PULLFOLDER=~/.singularity singularity pull docker://${DOCKER_IMAGE}
    # Write workflow option JSON file for singularity
    BACKEND=singularity
    TMP_WF_OPT=$RESULT_PREFIX.test_rna_wf_opt.json
    SINGULARITY_IMAGE=$(echo ${DOCKER_IMAGE} | sed 's/quay\.io\/encode-dcc\///g' | sed 's/:/-/' | sed 's/$/\.simg/')
    cat > $TMP_WF_OPT << EOM
    {
        "default_runtime_attributes" : {
            "singularity_container" : "~/.singularity/$SINGULARITY_IMAGE"
        }
    }
EOM
fi

java -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} -jar ${CROMWELL_JAR} run ${WDL} -i ${INPUT} -o ${TMP_WF_OPT} -m ${METADATA}

rm -f ${TMP_WF_OPT}