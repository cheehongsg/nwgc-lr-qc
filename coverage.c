#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1a-r1"
#endif

// TODO:
// verbose mode to report progress!!

/*********
 * Timer *
 *********/

int cov_verbose = 1;
double G_t_real;

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

/*********
 * coverage *
 *********/

static void print_usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Usage: coverage [options] <aligned_bam> <referece_fasta>\n\n");
    fprintf(fp, "Options:\n\n");
    fprintf(fp, "       -t INT        number of threads [%d]\n", 4);
    fprintf(fp, "       -q INT        have mapping quality >= INT [%d]\n", 20);
    fprintf(fp, "       -b INT        have base quality >= INT [%d]\n", 10);
    fprintf(fp, "       -c INT        max coverage at each locus <= INT [%d]\n", 250);
    fprintf(fp, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", cov_verbose);
    fprintf(fp, "       -r            row-based report; one field per row\n");
    fprintf(fp, "\n");
    return;
}

typedef struct coverage_t {
    uint64_t TOTAL_REF_BASES;
    uint64_t TOTAL_N_BASES;
    uint64_t TOTAL_READS;
    uint64_t READ_LENGTH;
    uint64_t TOTAL_BASES;
    uint64_t TOTAL_EXCLUDED_MAPQ;
    uint64_t TOTAL_EXCLUDED_BASEQ;
    uint64_t TOTAL_EXCLUDED_CAPPED;
    uint64_t TOTAL_EXCLUDED_BASES;
    uint64_t COV_1X_BASES;
    uint64_t COV_5X_BASES;
    uint64_t COV_10X_BASES;
    uint64_t COV_15X_BASES;
    uint64_t COV_20X_BASES;
    uint64_t COV_25X_BASES;
    uint64_t COV_30X_BASES;
    uint64_t COV_40X_BASES;
} coverage_t;

void rollupMetrics(bam1_t *bref, coverage_t *pMetrics, int32_t *starts, int32_t startsLen, uint32_t coverageCap) {
    int32_t depth = 0;
    for (int32_t c = 0; c < bref->core.l_qseq; ++c) {
        depth += starts[c];
        if (15==(bam_seqi(bam_get_seq(bref), c))) {
            // N base is omitted
            pMetrics->TOTAL_N_BASES++;
            pMetrics->TOTAL_EXCLUDED_BASES+=depth;
        } else {
            if (depth>coverageCap) {
                pMetrics->TOTAL_BASES+=coverageCap;
                uint32_t overage =  depth-coverageCap;
                pMetrics->TOTAL_EXCLUDED_BASES+=overage;
                pMetrics->TOTAL_EXCLUDED_CAPPED+=overage;
            } else {
                pMetrics->TOTAL_BASES+=depth;
            }
            if (depth>=40) {
                pMetrics->COV_40X_BASES++;
                pMetrics->COV_30X_BASES++;
                pMetrics->COV_25X_BASES++;
                pMetrics->COV_20X_BASES++;
                pMetrics->COV_15X_BASES++;
                pMetrics->COV_10X_BASES++;
                pMetrics->COV_5X_BASES++;
                pMetrics->COV_1X_BASES++;
            } else if (depth>=30) {
                pMetrics->COV_30X_BASES++;
                pMetrics->COV_25X_BASES++;
                pMetrics->COV_20X_BASES++;
                pMetrics->COV_15X_BASES++;
                pMetrics->COV_10X_BASES++;
                pMetrics->COV_5X_BASES++;
                pMetrics->COV_1X_BASES++;
            } else if (depth>=25) {
                pMetrics->COV_25X_BASES++;
                pMetrics->COV_20X_BASES++;
                pMetrics->COV_15X_BASES++;
                pMetrics->COV_10X_BASES++;
                pMetrics->COV_5X_BASES++;
                pMetrics->COV_1X_BASES++;
            } else if (depth>=20) {
                pMetrics->COV_20X_BASES++;
                pMetrics->COV_15X_BASES++;
                pMetrics->COV_10X_BASES++;
                pMetrics->COV_5X_BASES++;
                pMetrics->COV_1X_BASES++;
            } else if (depth>=15) {
                pMetrics->COV_15X_BASES++;
                pMetrics->COV_10X_BASES++;
                pMetrics->COV_5X_BASES++;
                pMetrics->COV_1X_BASES++;
            } else if (depth>=10) {
                pMetrics->COV_10X_BASES++;
                pMetrics->COV_5X_BASES++;
                pMetrics->COV_1X_BASES++;
            } else if (depth>=5) {
                pMetrics->COV_5X_BASES++;
                pMetrics->COV_1X_BASES++;
            } else if (depth>=1) {
                pMetrics->COV_1X_BASES++;
            }
        }
    }
    pMetrics->TOTAL_REF_BASES += bref->core.l_qseq;
}

int loadReferenceSequence(samFile *reffile, sam_hdr_t *ref_samhdr, bam1_t *bref, hts_pos_t tidlen, const char *tidname, coverage_t *pMetrics) {
    int ret = -1;
    
    while (0==(ret = sam_read1(reffile, ref_samhdr, bref))) {
        if (ret<0) {
            break;
        }
        // check that the sequence tally
        // or load the next
        if (0==strcmp(tidname, bam_get_qname(bref))) {
            if (tidlen == bref->core.l_qseq) {
                break;
            } else {
                fprintf(stderr, "[E::%s] Expecting %ld but counted %d for %s. Skipping to the next..\n", __func__, tidlen, bref->core.l_qseq, tidname);
                return -2;
            }
        } else {
            // we have to perform a quick update of the metrics
            pMetrics->TOTAL_REF_BASES += bref->core.l_qseq;
            for (int32_t c = 0; c < bref->core.l_qseq; ++c) {
                if (15==(bam_seqi(bam_get_seq(bref), c))) {
                    // N base is omitted
                    pMetrics->TOTAL_N_BASES++;
                }
            }

            if (cov_verbose >= 3) {
                fprintf(stderr, "[M::%s] Expecting %s, but seen %s (%d). Skipping to the next..\n", __func__, tidname, bam_get_qname(bref), bref->core.l_qseq);
            }
        }

    }
    return ret;
}

int loadRemainingReference(samFile *reffile, sam_hdr_t *ref_samhdr, bam1_t *bref, coverage_t *pMetrics) {
    int ret = -1;
    
    while (-1!=(ret = sam_read1(reffile, ref_samhdr, bref))) {
        if (ret>=0) {
            // we have to perform a quick update of the metrics
            pMetrics->TOTAL_REF_BASES += bref->core.l_qseq;
            for (int32_t c = 0; c < bref->core.l_qseq; ++c) {
                if (15==(bam_seqi(bam_get_seq(bref), c))) {
                    // N base is omitted
                    pMetrics->TOTAL_N_BASES++;
                }
            }
            if (cov_verbose >= 3) {
                fprintf(stderr, "[M::%s] Clearing remaining reference %s (%d).\n", __func__, bam_get_qname(bref), bref->core.l_qseq);
            }
        } else {
            // error!
            break;
        }
    }
    return ret;
}

// TODO: track timing ; progress indicator!
int main(int argc, char **argv) {
    int ret = EXIT_FAILURE;
    int c=0;
    int n_threads = 4;
    uint8_t minmapqual=20;
    uint8_t minbasequal=10;
    uint32_t coverageCap = 250;
    //uint32_t locusAccumulationCap = 100000;
    int rowBasedReport = 0;

    const char *inname = NULL;
    bam1_t *b = NULL;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    const char *tidname = NULL;
    int32_t tid=-1;
    hts_pos_t tidlen = 0;
    uint32_t *cigar = NULL;

    const char *refname = NULL;
    bam1_t *bref = NULL;
    samFile *reffile = NULL;
    sam_hdr_t *ref_samhdr = NULL;

    htsThreadPool tpool = {NULL, 0};

    int32_t i=0;
    uint8_t *quals = NULL;

    int32_t *starts = NULL;
    int32_t startsLen = 0;
    coverage_t metrics;
    memset(&metrics, 0, sizeof(metrics));

	while ((c = getopt(argc, argv, "rt:q:b:v:c:")) >= 0) {
		if (c == 't') n_threads = atoi(optarg), n_threads = n_threads > 1 ? n_threads : 1;
		else if (c == 'q') minmapqual = atoi(optarg), minmapqual = minmapqual >= 0 ? minmapqual : 20;
		else if (c == 'b') minbasequal = atoi(optarg), minbasequal = minbasequal >= 0 ? minbasequal : 10;
        else if (c == 'v') cov_verbose = atoi(optarg);
		else if (c == 'c') coverageCap = atoi(optarg), coverageCap = coverageCap >= 0 ? coverageCap : 250;
		else if (c == 'r') rowBasedReport = 1;
		else {
            print_usage(stderr);
            goto end;
		}
	}
	if (n_threads < 1) n_threads = 1;
	if (optind >= argc) {
        print_usage(stderr);
		goto end;
	}

    inname = argv[optind];
    if (optind + 1 < argc) {
        refname = argv[optind+1];
    } else {
        print_usage(stderr);
		goto end;
    }

    if (minbasequal<3) {
        if (cov_verbose >= 3) {
            fprintf(stderr, "[M::%s] Setting min base qual from %d to 3 for picard behaviour\n", __func__, minbasequal);
        }
        minbasequal = 3;
    }

	G_t_real = realtime();

    b = bam_init1();
    bref = bam_init1();

    if (!(infile = sam_open(inname, "r"))) {
        fprintf(stderr, "[E::%s] Could not open %s\n", __func__, inname);
        goto end;
    }
    if (infile->format.format != sam && infile->format.format != bam) {
        fprintf(stderr, "[E::%s] Alignment file - invalid format\n", __func__);
        goto end;
    }

    if (!(reffile = sam_open(refname, "r"))) {
        fprintf(stderr, "[E::%s] Could not open %s\n", refname, __func__);
        goto end;
    }
    if (reffile->format.format != fasta_format) {
        fprintf(stderr, "[E::%s] Reference fasta file - invalid format\n", __func__);
        goto end;
    }

    //create a pool of 4 threads
    if (!(tpool.pool = hts_tpool_init(n_threads))) {
        fprintf(stderr, "[E::%s] Failed to initialize the thread pool\n", __func__);
        goto end;
    }
    //share the pool with all the 3 files
    if (hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool) < 0 ||
        hts_set_opt(reffile, HTS_OPT_THREAD_POOL, &tpool) < 0) {
        fprintf(stderr, "[E::%s] Failed to set thread options\n", __func__);
        goto end;
    }

    if (!(in_samhdr = sam_hdr_read(infile))) {
        fprintf(stderr, "[E::%s] Failed to read header from file!\n", __func__);
        goto end;
    }

    if (!(ref_samhdr = sam_hdr_read(reffile))) {
        fprintf(stderr, "[E::%s] Failed to read header from file!\n", __func__);
        goto end;
    }

    while (sam_read1(infile, in_samhdr, b) >= 0) {

        // FIXME: track QCFAIL?
        if (BAM_FSECONDARY==(b->core.flag & BAM_FSECONDARY)) {
            continue;
        }

        // MAPQ filter
        if (b->core.qual<minmapqual) {
            cigar = bam_get_cigar(b);
            for (i = 0; i < b->core.n_cigar; ++i) {
                int op = bam_cigar_op(cigar[i]);
                int type = bam_cigar_type(op);
                if (3==type) {
                    // BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF
                    int olen = bam_cigar_oplen(cigar[i]);
                    metrics.TOTAL_EXCLUDED_MAPQ+=olen;
                    metrics.TOTAL_EXCLUDED_BASES+=olen;
                }
            }
            continue;
        }

        if (tid!=b->core.tid) {
            if (-1!=tid) {
                rollupMetrics(bref, &metrics, starts, startsLen, coverageCap);
            }
            tidname = sam_hdr_tid2name(in_samhdr, b->core.tid);
            tidlen = sam_hdr_tid2len(in_samhdr, b->core.tid);
            if (cov_verbose >= 3) {
                fprintf(stderr, "[M::%s] Processing %d - %s (%ld)\n", __func__, b->core.tid, tidname? tidname: "", tidlen);
            }
            tid=b->core.tid;

            // reset coverage depth for new reference
            if (0==startsLen) {
                if (cov_verbose >= 3) {
                    fprintf(stderr, "[M::%s] First coverage allocation %ld!\n", __func__, tidlen);
                }
                starts = calloc((tidlen+1), sizeof(int32_t));
                if (NULL==starts) {
                    fprintf(stderr, "[E::%s] Failed to allocate memory for coverage!\n", __func__);
                    goto end;
                }
                startsLen = tidlen;
            } else {
                if (tidlen>startsLen) {
                    // have to realloc
                    if (cov_verbose >= 3) {
                        fprintf(stderr, "[M::%s] Coverage re-allocation %ld!\n", __func__, tidlen);
                    }
                    starts = realloc(starts, (tidlen+1)*sizeof(int32_t));
                    if (NULL==starts) {
                        fprintf(stderr, "[E::%s] Failed to allocate memory for coverage!\n", __func__);
                        goto end;
                    }
                    startsLen = tidlen;
                    memset(starts, 0, (startsLen+1)*sizeof(int32_t));
                } else {
                    if (cov_verbose >= 3) {
                        fprintf(stderr, "[M::%s] Re-using coverage allocation %ld (vs %d)!\n", __func__, tidlen, startsLen);
                    }
                    // just initialize
                    memset(starts, 0, (startsLen+1)*sizeof(int32_t));
                }
            }
        
            // load the reference bases
            if (loadReferenceSequence(reffile, ref_samhdr, bref, tidlen, tidname, &metrics)<0) {
                fprintf(stderr, "[E::%s] Failed to read reference fasta entry for %d - %s (length: %ld)!\n", __func__, b->core.tid, tidname, tidlen);
                goto end;
            }
        }

        if (BAM_FSUPPLEMENTARY==(b->core.flag & BAM_FSUPPLEMENTARY)) {
            // TODO: have to process for split read
        } else {
            // consider primary reads
            metrics.TOTAL_READS++;
        }

        // FIXME: useOriginalQualities
        // rec.getOriginalBaseQualities() vs rec.getBaseQualities()

        // iterate the cigar for operation
        // own simplify version
        int refStart = b->core.pos;
        int refEnd = -1;
        int queryStart = 0;
        int queryEnd = -1;
        quals = bam_get_qual(b);
        cigar = bam_get_cigar(b);
        for (i = 0; i < b->core.n_cigar; ++i) {
            int op = bam_cigar_op(cigar[i]);
            int type = bam_cigar_type(op);
            int olen = bam_cigar_oplen(cigar[i]);
            if (3==type) {
                // BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF
                // moving query and reference
                // track query for base quality to decide usable
                // track reference to decide what is (starts) and (ends)
                refEnd = refStart + olen;
                queryEnd = queryStart + olen;

                int subRefPos = refStart;
                int subQueryPos = queryStart;
                while (subRefPos<refEnd) {
                    while (((quals[subQueryPos]<minbasequal) || (15==(bam_seqi(bam_get_seq(b), subQueryPos)))) && (subRefPos<refEnd)) {
                        subRefPos++;
                        subQueryPos++;
                        metrics.TOTAL_EXCLUDED_BASEQ++;
                        metrics.TOTAL_EXCLUDED_BASES++;
                    }
                    if (subRefPos<refEnd) {
                        int subStarts = subRefPos;
                        while ((quals[subQueryPos]>=minbasequal && (15!=(bam_seqi(bam_get_seq(b), subQueryPos)))) && (subRefPos<refEnd)) {
                            subRefPos++;
                            subQueryPos++;
                        }
                        starts[subStarts] += 1;
                        starts[subRefPos] += -1;
                    } else {
                        // no start found!
                    }
                }
                refStart = refEnd;
                queryStart = queryEnd;
            } else if (2==type) {
                // BAM_CDEL, BAM_CREF_SKIP
                // NOT moving query
                // move reference; but making stop and remembering possible start=pos+len?
                // we do not need to increase 'cov' as it is deletion
                refEnd = refStart + olen;
                
                refStart = refEnd;
            } else if (1==type) {
                // BAM_CINS, BAM_CSOFT_CLIP
                // move query, additional base(s) that are no present in the reference
                // NOT moving reference
                queryEnd = queryStart + olen;
                
                queryStart = queryEnd;
            } else {
                // BAM_CHARD_CLIP, BAM_CPAD, BAM_CBACK
                // do nothing
            }
        }
    }

    // TODO: handle the last entry of same tid!
    if (-1!=tid) {
        rollupMetrics(bref, &metrics, starts, startsLen, coverageCap);
    }
    // TODO: have to iterate the whole genome
    if (loadRemainingReference(reffile, ref_samhdr, bref, &metrics)<-1) {
        fprintf(stderr, "[E::%s] Failed to process remaining reference sequence(s)\n", __func__);
        goto end;
    }

    if (0!=rowBasedReport) {
        fprintf(stdout, "TOTAL_REF_BASES\t%lu\n", metrics.TOTAL_REF_BASES);
        fprintf(stdout, "TOTAL_N_BASES\t%lu\n", metrics.TOTAL_N_BASES);
        fprintf(stdout, "TOTAL_READS\t%lu\n", metrics.TOTAL_READS);
        fprintf(stdout, "TOTAL_BASES\t%lu\n", metrics.TOTAL_BASES);
        fprintf(stdout, "TOTAL_EXCLUDED_MAPQ\t%lu\n", metrics.TOTAL_EXCLUDED_MAPQ);
        fprintf(stdout, "TOTAL_EXCLUDED_BASEQ\t%lu\n", metrics.TOTAL_EXCLUDED_BASEQ);
        fprintf(stdout, "TOTAL_EXCLUDED_CAPPED\t%lu\n", metrics.TOTAL_EXCLUDED_CAPPED);
        fprintf(stdout, "TOTAL_EXCLUDED_BASES\t%lu\n", metrics.TOTAL_EXCLUDED_BASES);
        fprintf(stdout, "COV_1X_BASES\t%lu\n", metrics.COV_1X_BASES);
        fprintf(stdout, "COV_5X_BASES\t%lu\n", metrics.COV_5X_BASES);
        fprintf(stdout, "COV_10X_BASES\t%lu\n", metrics.COV_10X_BASES);
        fprintf(stdout, "COV_15X_BASES\t%lu\n", metrics.COV_15X_BASES);
        fprintf(stdout, "COV_20X_BASES\t%lu\n", metrics.COV_20X_BASES);
        fprintf(stdout, "COV_25X_BASES\t%lu\n", metrics.COV_25X_BASES);
        fprintf(stdout, "COV_30X_BASES\t%lu\n", metrics.COV_30X_BASES);
        fprintf(stdout, "COV_40X_BASES\t%lu\n", metrics.COV_40X_BASES);

        uint64_t effectiveGenome = metrics.TOTAL_REF_BASES - metrics.TOTAL_N_BASES;
        fprintf(stdout, "MEAN_COV\t%.6f\n", metrics.TOTAL_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_1X_BASES\t%.6f\n", metrics.COV_1X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_5X_BASES\t%.6f\n", metrics.COV_5X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_10X_BASES\t%.6f\n", metrics.COV_10X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_15X_BASES\t%.6f\n", metrics.COV_15X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_20X_BASES\t%.6f\n", metrics.COV_20X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_25X_BASES\t%.6f\n", metrics.COV_25X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_30X_BASES\t%.6f\n", metrics.COV_30X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "FRACTION_COV_40X_BASES\t%.6f\n", metrics.COV_40X_BASES * 1.0f / effectiveGenome);
    } else {
        fprintf(stdout, "TOTAL_REF_BASES\tTOTAL_N_BASES\tTOTAL_READS\tTOTAL_BASES");
        fprintf(stdout, "\tTOTAL_EXCLUDED_MAPQ\tTOTAL_EXCLUDED_BASEQ\tTOTAL_EXCLUDED_CAPPED\tTOTAL_EXCLUDED_BASES");
        fprintf(stdout, "\tCOV_1X_BASES\tCOV_5X_BASES\tCOV_10X_BASES\tCOV_15X_BASES");
        fprintf(stdout, "\tCOV_20X_BASES\tCOV_25X_BASES\tCOV_30X_BASES\tCOV_40X_BASES");
        fprintf(stdout, "\tMEAN_COV");
        fprintf(stdout, "\tFRACTION_1X_BASES\tFRACTION_5X_BASES\tFRACTION_10X_BASES\tFRACTION_15X_BASES");
        fprintf(stdout, "\tFRACTION_20X_BASES\tFRACTION_25X_BASES\tFRACTION_30X_BASES\tFRACTION_40X_BASES");
        fprintf(stdout, "\n");

        fprintf(stdout, "%lu\t%lu\t%lu\t%lu", metrics.TOTAL_REF_BASES, metrics.TOTAL_N_BASES, metrics.TOTAL_READS, metrics.TOTAL_BASES);
        fprintf(stdout, "\t%lu\t%lu\t%lu\t%lu", metrics.TOTAL_EXCLUDED_MAPQ, metrics.TOTAL_EXCLUDED_BASEQ, metrics.TOTAL_EXCLUDED_CAPPED, metrics.TOTAL_EXCLUDED_BASES);
        fprintf(stdout, "\t%lu\t%lu\t%lu\t%lu", metrics.COV_1X_BASES, metrics.COV_5X_BASES, metrics.COV_10X_BASES, metrics.COV_15X_BASES);
        fprintf(stdout, "\t%lu\t%lu\t%lu\t%lu", metrics.COV_20X_BASES, metrics.COV_25X_BASES, metrics.COV_30X_BASES, metrics.COV_40X_BASES);
        uint64_t effectiveGenome = metrics.TOTAL_REF_BASES - metrics.TOTAL_N_BASES;
        fprintf(stdout, "\t%.6f", metrics.TOTAL_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "\t%.6f\t%.6f", metrics.COV_1X_BASES * 1.0f / effectiveGenome, metrics.COV_5X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "\t%.6f\t%.6f", metrics.COV_10X_BASES * 1.0f / effectiveGenome, metrics.COV_15X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "\t%.6f\t%.6f", metrics.COV_20X_BASES * 1.0f / effectiveGenome, metrics.COV_25X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "\t%.6f\t%.6f", metrics.COV_30X_BASES * 1.0f / effectiveGenome, metrics.COV_40X_BASES * 1.0f / effectiveGenome);
        fprintf(stdout, "\n");
    }

    ret = EXIT_SUCCESS;

end:
    if (b) {
        bam_destroy1(b);
    }

    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }

    if (bref) {
        bam_destroy1(bref);
    }

    if (ref_samhdr) {
        sam_hdr_destroy(ref_samhdr);
    }
    if (reffile) {
        sam_close(reffile);
    }

    if (starts) {
        free(starts);
        starts = NULL;
    }

	if (EXIT_SUCCESS == ret) {
        int i;
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - G_t_real, cputime());
	}
    return ret;
}
