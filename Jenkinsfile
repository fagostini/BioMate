pipeline {
    agent { dockerfile {
        filename 'Dockerfile'
        dir 'build'
        label 'my-defined-label'
        additionalBuildArgs  ''
        args '-v /tmp:/tmp'
    } } 
    stages {
        stage('Build') { 
            steps {
                sh 'make test_mix' 
                stash(name: 'assets-files', includes: 'assets/*.csv') 
            }
        }
    }
}
