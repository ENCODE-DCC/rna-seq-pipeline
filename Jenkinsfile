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
                echo "Building image with tag $TAG"
                echo "Hello Ana and Idan!"
            }
        }
    }
}