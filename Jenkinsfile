pipeline {
    agent any 
    stages {
        stage('Build') { 
            agent { 
                dockerfile {
                    filename 'Dockerfile'
                }
            }
            steps {
                sh 'mkdir -p temp' 
                sh 'cp assets/SampleSheet_MixedIndexes.csv temp/SampleSheet.csv'
                sh 'biomate --verbose blabber --format fastq --sample-sheet temp/SampleSheet.csv --seq-number 10 --output temp --flowcell-id 20260310_LM43899_0385_A12GGASZR5'
                sh 'biomate --verbose fastrewind --input-path temp --output-path temp --threads 8'
                sh 'ulimit -Sn 65535 && assets/bcl-convert --output-directory temp/Demultiplexing --bcl-input-directory temp --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true --bcl-validate-sample-sheet-only true'
                sh 'ulimit -Sn 65535 && assets/bcl-convert --output-directory temp/Demultiplexing --bcl-input-directory temp --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true'
                stash(name: 'assets-files', includes: 'assets/*') 
            }
        }
    }
}
