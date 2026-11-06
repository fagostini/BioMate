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
                sh 'biomate --verbose blabber --format fastq --sample-sheet temp/SampleSheet.csv --seq-number 10 --taint --taint-count 5 --random-seed 42 --output temp --flowcell-id 20260310_LM43899_0385_A12GGASZR5 > temp/blabber.out 2> temp/blabber.err'
                sh 'biomate --verbose fastrewind --input-path temp --output-path temp/20260310_LM43899_0385_A12GGASZR5 --flowcell-id 20260310_LM43899_0385_A12GGASZR5 --threads 8 > temp/fastrewind.out 2> temp/fastrewind.err'
                sh 'ulimit -Sn 65535 2>/dev/null || true; if [ "$(ulimit -Sn)" -lt 4096 ]; then echo "Error: bcl-convert requires at least 4096 open files, but the soft nofile limit is only $(ulimit -Sn) (hard limit: $(ulimit -Hn)); increase the container nofile limit (e.g. docker --ulimit nofile=65535:65535) and retry"; exit 1; fi; assets/bcl-convert --output-directory temp/Demultiplexed --bcl-input-directory temp/20260310_LM43899_0385_A12GGASZR5 --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true --bcl-validate-sample-sheet-only true > temp/bcl-validate.out 2> temp/bcl-validate.err'
                sh 'ulimit -Sn 65535 2>/dev/null || true; if [ "$(ulimit -Sn)" -lt 4096 ]; then echo "Error: bcl-convert requires at least 4096 open files, but the soft nofile limit is only $(ulimit -Sn) (hard limit: $(ulimit -Hn)); increase the container nofile limit (e.g. docker --ulimit nofile=65535:65535) and retry"; exit 1; fi; assets/bcl-convert --output-directory temp/Demultiplexed --bcl-input-directory temp/20260310_LM43899_0385_A12GGASZR5 --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true --force > temp/bcl-convert.out 2> temp/bcl-convert.err'
                sh 'find temp/Demultiplexed -name "*.fastq.gz" | grep -v -e "Undetermined" -e "_I1_" -e "_I2_" | xargs -r zgrep -c ^@ || true > temp/fastq-counts.txt 2> temp/fastq-counts.err'
                sh 'find temp/Demultiplexed -name "*.fastq.gz" | grep "Undetermined" | grep -v -e "_I1_" -e "_I2_" | xargs -r zgrep -c ^@ || true > temp/undetermined-counts.txt 2> temp/undetermined-counts.err'
                sh 'echo "Verifying undetermined read counts (5 expected per indexed lane)..."; for spec in L001:0 L002:5 L003:5 L004:5 L005:5 L006:5 L007:5 L008:5; do lane=${spec%%:*}; expected=${spec##*:}; file=$(find temp/Demultiplexed -name "Undetermined_S0_${lane}_R1_*.fastq.gz" | head -n 1); if [ -n "$file" ]; then actual=$(zgrep -c ^@ "$file"); else actual=0; fi; if [ "$actual" != "$expected" ]; then echo "ERROR: lane ${lane} has ${actual} undetermined R1 reads, expected ${expected}"; exit 1; fi; echo "  lane ${lane}: ${actual}/${expected} undetermined R1 reads"; done'
                sh 'assets/compare_results.sh && echo "All files pairwise comparisons were successful!" || echo "WARNING: Some pairwise comparisons yield different results!"'
                stash(name: 'assets-files', includes: 'assets/*') 
            }
        }
    }
}
