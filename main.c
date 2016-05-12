#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/khash.h"
#include "htslib/htslib/hts.h"

//argv[1] = bam file argv[2] = str file
int main(int argc, char** argv) {

	hts_itr_t *iter = NULL;
	hts_idx_t *idx = NULL;
	htsFile *infile = NULL;
	// Structure for one alignment
	bam1_t *b = NULL;
	bam_hdr_t *header = NULL;

	FILE *strFile = fopen(argv[2], "r");
	assert(strFile != NULL);

	infile = hts_open(argv[1], "r");
	assert(infile != NULL);

	idx = sam_index_load(infile,  argv[1]);
	assert(idx != NULL);

	header = bam_hdr_read((infile->fp).bgzf);
	assert(header != NULL);

	// printf("Text of header: %s\n", header->text);


	// iter  = sam_itr_querys(idx, header, ".");
	// assert(iter != NULL);

	b = bam_init1();

	// iterate over str file, get only chr-start-end positions of the str
	int strLine[3], strCount = 0;
	char line[100];
	char buffer[100];

	while(fgets(line, sizeof line, strFile) != NULL) {
		int count = 0;
		sscanf(line, "%d %d %d", strLine, strLine + 1, strLine + 2);
		// printf("strCount:%d %d %d %d\n", strCount, strLine[0], strLine[1], strLine[2]);
		strCount++;

		//search the STRs in bam file
		iter  = sam_itr_queryi(idx, strLine[0], strLine[1], strLine[2]);
		assert(iter != NULL);

		// if (strLine[0] == 20)
		// 	continue;

		while(sam_itr_next(infile, iter, b) >= 0) {
			count++;
			printf("*******%d\n", count);
			printf("STR: %d %d %d\n", strLine[0], strLine[1], strLine[2]);
			printf("QNAME: %s\n", bam_get_qname(b));
			printf("TID: %d\n", (b)->core.tid);
			printf("POS: %d\n", (b)->core.pos);
			printf("SEQ LENGTH: %d\n", (b)->core.l_qseq);
			printf("END: %d\n", (b)->core.pos + (b)->core.l_qseq);

		}
		// printf("For str %d: %d reads\n", strCount, count);
	}

	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(infile);
	fclose(strFile);

	return 0;
}
