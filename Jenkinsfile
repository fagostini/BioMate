pipeline {
    agent any 
    stages {
        stage('Build') { 
            steps {
                sh 'make test_mix' 
                stash(name: 'assets-files', includes: 'assets/*.csv') 
            }
        }
    }
}
