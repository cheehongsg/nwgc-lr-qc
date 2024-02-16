# Coverage
Coverage computes basic coverage similar to picard CollectWgsMetrics but at the speed of mosdepth

## htslib static library
Coverage uses htslib static library for file parsing.

```
wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
tar xjvf htslib-1.19.1.tar.bz2
cd htslib-1.19.1
./configure
make
make install
```

## Compiling coverage
```
gcc -O2 -Wall coverage.c ./libhts.a -lz -lbz2 -llzma -lcurl -lpthread -lcrypto -lm -o coverage
```

## Command line
```
Usage: coverage [options] <aligned_bam> <referece_fasta>

Options:

       -t INT        number of threads [4]
       -q INT        have mapping quality >= INT [20]
       -b INT        have base quality >= INT [10]
       -c INT        max coverage at each locus <= INT [250]
       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [1]
       -r            row-based report; one field per row
```

Recommended ```-b 20``` for PacBio data and ```-b 10``` for ONT data.

## Examples
```bash
# uses additional 4 threads and min base quality of 20 for pacbio aligned data on GRCh38
# and report one mtric per column
coverage -b 20 -t 4 pacbio.aligned.bam human_GRCh38_no_alt_analysis_set.fasta

# uses additional 6 threads and min base quality of 10 for ont aligned data on GRCh38
# and report one metric per row
coverage -b 10 -t 6 -r ont.aligned.bam human_GRCh38_no_alt_analysis_set.fasta
```

## Output
```bash
###
# ont.aligned.bam is 86GB
coverage -t 4 ont.aligned.bam human_GRCh38_no_alt_analysis_set.fasta

# output:
# stdout has the metrics output
TOTAL_REF_BASES TOTAL_N_BASES   TOTAL_READS     TOTAL_BASES     TOTAL_EXCLUDED_MAPQ     TOTAL_EXCLUDED_BASEQ    TOTAL_EXCLUDED_CAPPED   TOTAL_EXCLUDED_BASES    COV_1X_BASESCOV_5X_BASES     COV_10X_BASES   COV_15X_BASES   COV_20X_BASES   COV_25X_BASES   COV_30X_BASES   COV_40X_BASES   MEAN_COV        FRACTION_1X_BASES       FRACTION_5X_BASES   FRACTION_10X_BASES       FRACTION_15X_BASES      FRACTION_20X_BASES      FRACTION_25X_BASES      FRACTION_30X_BASES      FRACTION_40X_BASES
3099922541      165046090       4869311 78070833404     4586842894      4465752587      736541112       9790105386      2864730897      2844875674      2835081584      2800484472   2574555611      1902073424      984822613       87317824        26.601063       0.976099        0.969334        0.965997        0.954209        0.877228        0.648093     0.335558        0.029752
# stderr has the execution and timing record
[main] Version: 0.0.1a-r1
[main] CMD: coverage -t 4 ont.aligned.bam human_GRCh38_no_alt_analysis_set.fasta
[main] Real time: 303.929 sec; CPU: 1322.831 sec


###
# alternatively, the same file (ont.aligned.bam, 86GB) with metric per row
coverage -t 6 -r ont.aligned.bam human_GRCh38_no_alt_analysis_set.fasta

# output:
# stdout has the metrics output
TOTAL_REF_BASES 3099922541
TOTAL_N_BASES   165046090
TOTAL_READS     4869311
TOTAL_BASES     78070833404
TOTAL_EXCLUDED_MAPQ     4586842894
TOTAL_EXCLUDED_BASEQ    4465752587
TOTAL_EXCLUDED_CAPPED   736541112
TOTAL_EXCLUDED_BASES    9790105386
COV_1X_BASES    2864730897
COV_5X_BASES    2844875674
COV_10X_BASES   2835081584
COV_15X_BASES   2800484472
COV_20X_BASES   2574555611
COV_25X_BASES   1902073424
COV_30X_BASES   984822613
COV_40X_BASES   87317824
MEAN_COV        26.601063
FRACTION_COV_1X_BASES   0.976099
FRACTION_COV_5X_BASES   0.969334
FRACTION_COV_10X_BASES  0.965997
FRACTION_COV_15X_BASES  0.954209
FRACTION_COV_20X_BASES  0.877228
FRACTION_COV_25X_BASES  0.648093
FRACTION_COV_30X_BASES  0.335558
FRACTION_COV_40X_BASES  0.029752

# stderr has the execution and timing record
[main] Version: 0.0.1a-r1
[main] CMD: coverage -t 6 -r ont.aligned.bam human_GRCh38_no_alt_analysis_set.fasta
[main] Real time: 211.501 sec; CPU: 1318.548 sec
```
