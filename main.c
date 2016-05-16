#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/khash.h"
#include "htslib/htslib/hts.h"

#define MAX_STR_SIZE 1000

// Structure for pair reads information: sequences, quality, positions/names/tids, cigar, flag.
typedef struct {
	unsigned char *qname;
    int32_t tid, mtid;
    int32_t pos, mpos;
    uint32_t flag:16, mflag:16, cigar:16, mcigar:16;
    int32_t l_qseq, ml_qseq;
    uint8_t *seq, *mseq;
    uint8_t *qual, *mqual;
} bam1_pair;

/* Decode 4-bit encoded bases to their corresponding characters */
char base_as_char( int base_as_int)
{
	if( base_as_int == 1)
	{
		return 'A';
	}
	else if( base_as_int == 2)
	{
		return 'C';
	}
	else if( base_as_int == 4)
	{
		return 'G';
	}
	else if( base_as_int == 8)
	{
		return 'T';
	}
	else if( base_as_int == 15)
	{
		return 'N';
	}
	return 0;
}

char *seqToString(int32_t l_qseq, uint8_t *seq) {
	char sequence[MAX_STR_SIZE];
	char *read = (char *) malloc(MAX_STR_SIZE * sizeof(char));
	int i;
	strncpy(sequence, seq, l_qseq);
	sequence[l_qseq] = '\0';

	for(i = 0; i < strlen(sequence); i++) {
		char next_char = base_as_char(bam_seqi(sequence, i));
		read[i] = next_char;
	}
	read[i] = '\0';
	return read;
}

// char *qualToString(uint8_t *qual, int32_t *l_qseq){
// 	char qual_str[MAX_STR_SIZE];
// 	strncpy(qual_str, qual, l_qseq);
// 	qual[l_qseq] = '\0';

// 	for(int i = 0; i < strlen(qual_str); i++) {
// 		qual_str[i] = qual_str[i] + 33;
// 	}
// 	return qual_str;
// }

char *serialize_bam1_pair(bam1_pair pair) {
	char* str = (char *) malloc(sizeof(char) * 10000);
	char *sequence = seqToString(pair.l_qseq, pair.seq);
	sprintf(str, "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s", pair.qname, pair.flag, pair.tid, pair.pos, pair.l_qseq, strlen(sequence), sequence, pair.mflag, pair.mtid, pair.mpos, pair.ml_qseq, pair.mseq);
	free(sequence);
	return str;
}

void print(FILE *out, bam1_pair pairs[], int length) {
	for(int i = 0; i < length; i++) {
		char *line = serialize_bam1_pair(pairs[i]);
		fprintf(out,"%s\n", line);
		free(line);
	}
}

//argv[1] = bam file argv[2] = str file
int main(int argc, char** argv) {

	hts_itr_t *iter = NULL;
	hts_idx_t *idx, *mate_idx = NULL;
	htsFile *infile, *mate_file = NULL;
	// Structure for one alignment
	bam1_t *b = NULL;
	bam_hdr_t *header = NULL;

	FILE *strFile = fopen(argv[2], "r");
	assert(strFile != NULL);

	// file for printing pair structure
	FILE *out = fopen("out.txt", "w");
	assert(out != NULL);

	infile = hts_open(argv[1], "r");
	assert(infile != NULL);

	mate_file = hts_open(argv[1], "r");
	assert(mate_file != NULL);

	idx = sam_index_load(infile,  argv[1]);
	assert(idx != NULL);

	mate_idx = sam_index_load(mate_file,  argv[1]);
	assert(mate_idx != NULL);

	header = bam_hdr_read((infile->fp).bgzf);
	assert(header != NULL);

	// printf("Text of header: %s\n", header->text);

	b = bam_init1();

	// iterate over str file, get only chr-start-end positions of the str
	int strLine[3], strCount = 0;
	char line[100];
	char buffer[100];

	while(fgets(line, sizeof line, strFile) != NULL) {
		// printf("STR COUNT: %d\n", strCount);
		int count = 0, pair_idx = 0;
		bam1_pair pairs[10000];
		sscanf(line, "%d %d %d", strLine, strLine + 1, strLine + 2);
		// printf("strCount:%d %d %d %d\n", strCount, strLine[0], strLine[1], strLine[2]);
		strCount++;

		//search the STRs in bam file
		iter  = sam_itr_queryi(idx, strLine[0], strLine[1], strLine[2]);
		assert(iter != NULL);

		while(sam_itr_next(infile, iter, b) >= 0) {
			bam1_t *m = bam_init1();
			count++;

			//Each segment properly aligned(bit 2) and unmapped bit is not set(bit 3)
			if(((b)->core.flag&BAM_FPROPER_PAIR) && !((b)->core.flag&BAM_FUNMAP)) {

				// printf("READ QNAME: %s ", bam_get_qname(b));
				// printf("READ FLAG: %d \n", (b)->core.flag);
				// printf("READ POS: %d\n", (b)->core.pos);
				// printf("MPOS: %d\n", (b)->core.mpos);

				hts_itr_t *mate_iter  = sam_itr_queryi(mate_idx, (b)->core.mtid, (b)->core.mpos, (b)->core.mpos + 1);
				assert(mate_iter != NULL);

				int count_mates = 0;
				while(sam_itr_next(mate_file, mate_iter, m) > 0) {
					if (strcmp(bam_get_qname(b), bam_get_qname(m)) != 0) continue;
					// printf("MATE QNAME: %s ", bam_get_qname(m));
					// printf("MATE FLAG: %d \n", (m)->core.flag);
					// printf("MATE POS: %d\n", (m)->core.pos);
					count_mates++;
				}

				if(count_mates == 1) {

					bam1_pair pair;
					pair.qname = (unsigned char *) bam_get_qname(b);
					pair.flag = (b)->core.flag;
					pair.mflag = (m)->core.flag;
					pair.tid = (b)->core.tid;
					pair.mtid = (m)->core.tid;
					pair.pos = (b)->core.pos;
					pair.mpos = (m)->core.pos;
					pair.l_qseq = (b)->core.l_qseq;
					pair.ml_qseq = (m)->core.l_qseq;
					pair.seq = bam_get_seq(b);
					pair.mseq = bam_get_seq(m);
					pair.qual = bam_get_qual(b);
					pair.mqual = bam_get_qual(m);
					pairs[pair_idx++] = pair;

				}
			}
		}
		if(pair_idx != 0) {
			printf("How many pairs for an STR: %d\n", pair_idx);
			print(out, pairs, pair_idx);
			return 0;
		}
	}

	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(infile);
	sam_close(mate_file);
	fclose(strFile);
	fclose(out);

	return 0;
}
