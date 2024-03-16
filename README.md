## Code for AMAR-seq: automated multimodal sequencing of DNA methylation, chromatin accessibility, and RNA expression with single-cell resolution

## Running the Code
### Input files
1. The fq.gz file from AMAR-seq sequencing
2. Bismark reference file

### Processing steps

1. Specify the input sequence file and reference in nomve_process_bs.py
2. run python -u nomve_process_bs.py
3. The CpG and GpC profiles will be obtained separately.