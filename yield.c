#include "htslib/sam.h"

#include <sys/resource.h>
#include <sys/time.h>

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1a-r1"
#endif

/*********
 * Timer *
 *********/

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
 * yield *
 *********/

static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: yield <bamfile>\n\
Report yield metrics\n");
    return;
}

typedef struct yield_t {
    uint64_t TOTAL_READS;
    uint64_t PF_READS;
    uint64_t READ_LENGTH;
    uint64_t TOTAL_BASES;
    uint64_t PF_BASES;
    uint64_t Q10_BASES;
    uint64_t PF_Q10_BASES;
    uint64_t Q20_BASES;
    uint64_t PF_Q20_BASES;
    uint64_t Q30_BASES;
    uint64_t PF_Q30_BASES;
    uint64_t Q40_BASES;
    uint64_t PF_Q40_BASES;
    uint64_t Q20_EQUIVALENT_YIELD;
    uint64_t PF_Q20_EQUIVALENT_YIELD;
} yield_t;

// TODO: track timing ; progress indicator!
int main(int argc, char **argv) {
    int ret = EXIT_FAILURE;
    const char *inname = NULL;
    bam1_t *b = NULL;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;

    uint16_t isPfRead=0;
    int32_t c=0;
    uint8_t *quals = NULL;
    uint8_t qual=0;

    yield_t metrics;
    memset(&metrics, 0, sizeof(metrics));

    if (argc != 2) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];

	G_t_real = realtime();

    b = bam_init1();

    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }

    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    while (sam_read1(infile, in_samhdr, b) >= 0) {

        if (BAM_FSECONDARY==(b->core.flag & BAM_FSECONDARY)) {
            continue;
        }

        if (BAM_FSUPPLEMENTARY==(b->core.flag & BAM_FSUPPLEMENTARY)) {
            continue;
        }

        metrics.TOTAL_READS++;
        metrics.TOTAL_BASES += b->core.l_qseq;

        isPfRead = (b->core.flag & BAM_FQCFAIL);
        if (0==isPfRead) {
            metrics.PF_READS++;
            metrics.PF_BASES += b->core.l_qseq;
        }

        // FIXME: useOriginalQualities
        // rec.getOriginalBaseQualities() vs rec.getBaseQualities()

        // infile->format.format == fastq_format
        quals = bam_get_qual(b);
        for (c = 0; c < b->core.l_qseq; ++c) {
            qual=quals[c];

            metrics.Q20_EQUIVALENT_YIELD += qual;
            if (0==isPfRead) {
                metrics.PF_Q20_EQUIVALENT_YIELD += qual;
            }
            if (qual >= 40) {
                metrics.Q10_BASES++;
                metrics.Q20_BASES++;
                metrics.Q30_BASES++;
                metrics.Q40_BASES++;
                if (0==isPfRead) {
                    metrics.PF_Q10_BASES++;
                    metrics.PF_Q20_BASES++;
                    metrics.PF_Q30_BASES++;
                    metrics.PF_Q40_BASES++;
                }
            } else if (qual >= 30) {
                metrics.Q10_BASES++;
                metrics.Q20_BASES++;
                metrics.Q30_BASES++;
                if (0==isPfRead) {
                    metrics.PF_Q10_BASES++;
                    metrics.PF_Q20_BASES++;
                    metrics.PF_Q30_BASES++;
                }
            } else if (qual >= 20) {
                metrics.Q10_BASES++;
                metrics.Q20_BASES++;
                if (0==isPfRead) {
                    metrics.PF_Q10_BASES++;
                    metrics.PF_Q20_BASES++;
                }
            } else if (qual >= 10) {
                metrics.Q10_BASES++;
                if (0==isPfRead) {
                    metrics.PF_Q10_BASES++;
                }
            }
        }
    }

    metrics.Q20_EQUIVALENT_YIELD = metrics.Q20_EQUIVALENT_YIELD / 20;
    metrics.PF_Q20_EQUIVALENT_YIELD = metrics.PF_Q20_EQUIVALENT_YIELD / 20;
    metrics.READ_LENGTH = (metrics.TOTAL_READS == 0) ? 0 : (int) (metrics.TOTAL_BASES / metrics.TOTAL_READS);

    fprintf(stdout, "TOTAL_READS\tPF_READS\tREAD_LENGTH\tTOTAL_BASES\tPF_BASES\tQ10_BASES\tPF_Q10_BASES\tQ20_BASES\tPF_Q20_BASES\tQ30_BASES\tPF_Q30_BASES\tQ40_BASES\tPF_Q40_BASES\tQ20_EQUIVALENT_YIELD\tPF_Q20_EQUIVALENT_YIELD\n");
    fprintf(stdout, "%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n", 
        metrics.TOTAL_READS, metrics.PF_READS, metrics.READ_LENGTH, metrics.TOTAL_BASES, metrics.PF_BASES,
        metrics.Q10_BASES, metrics.PF_Q10_BASES,
        metrics.Q20_BASES, metrics.PF_Q20_BASES, metrics.Q30_BASES, metrics.PF_Q30_BASES,
        metrics.Q40_BASES, metrics.PF_Q40_BASES,
        metrics.Q20_EQUIVALENT_YIELD, metrics.PF_Q20_EQUIVALENT_YIELD);

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
