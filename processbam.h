#ifndef __PROCESSBAM
#define __PROCESSBAM

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "common.h"

/* Sample this many fragments to calculate avg/median/std per library */
#define SAMPLEFRAG 1000000 

/* Maximum sequence/quality length */
#define MAX_SEQ 1000

typedef struct _bam_info
{
	htsFile* bam_file; /* file pointer to the BAM file */
  	char* sample_name; /* name of the sample, parsed from SM in the BAM header */
	int num_libraries; /* number of libraries, counted from the RG tags in the BAM header */
	struct library_properties** libraries; /* each library_properties struct holds statistical/other info */
} bam_info;

/* Function Prototypes */
void load_bam( bam_info* in_bam, char* path);

void print_bam( bam_info* in_bam);


/* BAM Utility functions */
void get_sample_name( bam_info* in_bam, char* header_text);



#endif
