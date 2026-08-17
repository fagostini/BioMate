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
                sh 'rm -fr temp/*'
                sh 'cp assets/SampleSheet_MixedIndexes.csv temp/SampleSheet.csv'
                sh 'biomate --verbose blabber --format fastq --sample-sheet temp/SampleSheet.csv --seq-number 10 --output temp --flowcell-id 20260310_LM43899_0385_A12GGASZR5 > temp/blabber.out 2> temp/blabber.err'
                sh 'biomate --verbose fastrewind --input-path temp --output-path temp --threads 8 > temp/fastrewind.out 2> temp/fastrewind.err'
                sh 'ulimit -Sn 65535 && assets/bcl-convert --output-directory temp/Demultiplexing --bcl-input-directory temp --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true --bcl-validate-sample-sheet-only true > temp/bcl-validate.out 2> temp/bcl-validate.err'
                sh 'ulimit -Sn 65535 && assets/bcl-convert --output-directory temp/Demultiplexing --bcl-input-directory temp --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true > temp/bcl-convert.out 2> temp/bcl-convert.err'
                sh 'find temp/Demultiplexing -name "*.fastq.gz" | grep -v -e "Undetermined" -e "_I1_" -e "_I2_" | xargs -r zgrep -c ^@ || true > temp/fastq-counts.txt 2> temp/fastq-counts.err'
                sh 'find temp/Demultiplexing -name "*.fastq.gz" | grep "Undetermined" | grep -v -e "_I1_" -e "_I2_" | xargs -r zgrep -c ^@ || true > temp/undetermined-counts.txt 2> temp/undetermined-counts.err'
                sh 'assets/compare_results.sh && echo "All files pairwise comparisons were successful!" || echo "WARNING: Some pairwise comparisons yield different results!"'
                stash(name: 'assets-files', includes: 'assets/*') 
            }
        }
    }
}
