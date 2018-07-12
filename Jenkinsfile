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
                echo "testing jenkins 2"
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
                echo "Building image with tag $TAG"
                slackSend (color: '#7CFC00', message: "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}.")
                slackSend "The images will be tagged as $TAG"
                sh "docker pull quay.io/encode-dcc/rna-seq-pipeline:template"
                sh "docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io"
                sh "docker build --cache-from quay.io/encode-dcc/rna-seq-pipeline:template -t rna-seq-pipeline ."
                sh "docker tag rna-seq-pipeline $TAG"
                sh "docker push $TAG"
                sh "docker logout"
                slackSend "Finished pushing $TAG into image repo."
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