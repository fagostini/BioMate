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
                sh 'make test_mix BCLCONVERT="assets/bcl-convert"' 
                stash(name: 'assets-files', includes: 'assets/*') 
            }
        }
    }
}
