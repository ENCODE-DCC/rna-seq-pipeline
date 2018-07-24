pipeline{

    agent none
    
    environment{
        QUAY_USER = credentials('quay-robot')
        QUAY_PASS = credentials('quay-robot-token')
    }

    stages{
        stage('Tag Image'){
            agent {label 'master-builder'}
            steps{
                echo "Building image tag.."
                script{
                    TAG = sh([script: "echo quay.io/encode-dcc/rna-seq-pipeline:${env.BRANCH_NAME}_${env.BUILD_NUMBER}", returnStdout: true]).trim()
                }
                echo "The image tag that I just built is $TAG"
            }
        }
        stage('Build Image'){
            agent{label 'slave-w-docker-cromwell-60GB-ebs'}
            steps{
                echo "Building image with tag $TAG from git commit ${env.GIT_COMMIT}"
                echo "The git commit hash will be available within the image in environment variable GIT_HASH"
                echo "The branch name will be available within the image in environment variable BUILD_BRANCH"
                echo "The image tag will be available within the image in environment variable MY_TAG"
                slackSend (color: '#7CFC00', message: "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}.")
                slackSend "The images will be tagged as $TAG"
                sh "docker pull quay.io/encode-dcc/rna-seq-pipeline:template"
                sh "docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io"
                sh "docker build --cache-from quay.io/encode-dcc/rna-seq-pipeline:template --build-arg GIT_COMMIT_HASH=${env.GIT_COMMIT} --build-arg BRANCH=${env.BRANCH_NAME} --build-arg BUILD_TAG=$TAG -t $TAG -t quay.io/encode-dcc/rna-seq-pipeline:template ."
                sh "docker push $TAG"
                sh "docker push quay.io/encode-dcc/rna-seq-pipeline:template"
                sh "docker logout"
                slackSend "Finished pushing $TAG into image repo."
            }
        }
        stage('Task level tests'){
            agent{label 'slave-w-docker-cromwell-60GB-ebs'}
            steps{
                echo "Running task level tests."
                sh "pwd"
                sh "ls -l"
                echo "Fetching chromosome 19 restricted index file for STAR from Google Cloud"
                sh "curl https://storage.googleapis.com/star-rsem-runs/reference-genomes/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz -o test_data/GRCh38_v24_ERCC_phiX_starIndex_chr19only.tgz"
                sh """test_wdl/test.sh test_wdl/verify_install.wdl test_wdl/verify_install_input.json $TAG
                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < verify_install.result.json
                   """
            }
        }
    }

    post{
        success{
                echo "Post build actions that run on success"
                slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished with"
                slackSend (color: '#7cfc00', message: "SUCCESS")
                slackSend "For details, visit ${env.BUILD_URL}"
            }
        failure{
                echo "Post build actions that run on failure"
                slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished with"
                slackSend (color: '#FF0000', message: "FAILURE")
                slackSend "For details, visit ${env.BUILD_URL}"
            }
        }
    }