pipeline {
    agent any 
    stages {
        stage('Build') { 
            agent { 
                dockerfile {
                    filename 'Dockerfile'
                    label 'ubuntu-docker-python-venv'
                }
            }
            steps {
                sh 'make test_mix' 
                stash(name: 'assets-files', includes: 'assets/*.csv') 
            }
        }
    }
}
