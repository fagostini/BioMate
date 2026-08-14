pipeline {
    agent any 
    stages {
        stage('Build') { 
            steps {
                sh 'make test_mix' 
                stash(name: 'compiled-results', includes: 'assets/*') 
            }
        }
    }
}
