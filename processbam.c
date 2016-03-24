#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

/* tardis headers */
#include "processbam.h"

void load_bam( bam_info* in_bam, char* path)
{
	/* Variables */
	htsFile* bam_file;
	bam_hdr_t* bam_header;
	hts_idx_t* bam_index;

	fprintf( stderr, "Processing BAM file %s.\n", path);


	/* Open the BAM file for reading. htslib automatically detects the format
		of the file, so appending "b" after "r" in mode is redundant. */
	bam_file = safe_hts_open( path, "r");
	bam_index = sam_index_load( bam_file, path);

	if (bam_index == NULL){
	  fprintf(stderr, "BAM index not found.\n");
	  exit (EXIT_COMMON);
	}
	/* Read in BAM header information */
	bam_header = bam_hdr_read( ( bam_file->fp).bgzf);

	get_sample_name( in_bam, bam_header->text); 
	in_bam->bam_file = bam_file;
	in_bam->bam_index = bam_index;
	in_bam->bam_header = bam_header;

}

void read_alignment( bam_info* in_bam, parameters *params)
{
	bam1_core_t bam_alignment_core;
	bam1_t*	bam_alignment;	
	hts_idx_t* bam_index;
	int return_value;
	int i;
	int j;

	char sequence[MAX_SEQ];
	char read[MAX_SEQ];
	char qual[MAX_SEQ];
	char read2[MAX_SEQ];

	char next_char;
	const uint32_t *cigar;
	int n_cigar;
	htsFile *bam_file;
	
	bam_hdr_t* bam_header;
	char rand_loc[MAX_SEQ];	
	int result;
	char md[MAX_SEQ];
	int chrom_id;
	int start; int end;
	int maps_to_test = params->maps_to_test;
	char map_chr[MAX_SEQ];
	int map_tid;
	int map_loc;
	char *ref_seq;
	char ref_seq2[MAX_SEQ];
	int loc_len;
	int clipped;
	int cigar_add_len;

	bam_file = in_bam->bam_file;
	bam_header = in_bam->bam_header;
	bam_index = in_bam->bam_index;

	bam_alignment = bam_init1();
	j=0;
	
	while (j<maps_to_test){
      	  
	  //	  return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);

	  //	  if (rand() % 181 != 0)
	  // continue;

	  
	  if (bam_index == NULL){
	    fprintf(stderr, "BAM index not found.\n");
	    exit (EXIT_COMMON);
	  }

	  if (bam_header == NULL){
	    fprintf(stderr, "BAM header not found.\n");
	    exit (EXIT_COMMON);
	  }

	  
	  chrom_id = rand() % (params->num_chrom);
	  
	  start = (rand() % (params->chrom_lengths[chrom_id])) - 1000;
	  end = start + 1000;

	  sprintf(rand_loc, "%s:%d-%d", params->chrom_names[chrom_id], start, end);

	  hts_itr_t *iter=sam_itr_queryi(bam_index, chrom_id, start, end);

	  if (iter == NULL) { // region invalid or reference name not found

	    hts_itr_destroy(iter);
	    continue;
	  }

	  result = sam_itr_next(bam_file, iter, bam_alignment);

	  if (result < 0){
	    //printf("res %s\n", rand_loc);
	    hts_itr_destroy(iter);
	    continue;
	  }

	  bam_alignment_core = bam_alignment->core;
	  
	  if (bam_alignment_core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) // skip secondary and supplementary alignments
	    {
	      hts_itr_destroy(iter);	  
	      continue;
	    }
	  if (bam_aux_get(bam_alignment, "MD") == NULL){
	      hts_itr_destroy(iter);	  
	      continue;
	  }

	  cigar = bam_get_cigar(bam_alignment);
	  n_cigar = bam_alignment_core.n_cigar;
	  strncpy( sequence, bam_get_seq( bam_alignment), bam_alignment_core.l_qseq);
	  sequence[bam_alignment_core.l_qseq] = '\0';
	  strncpy( qual, bam_get_qual( bam_alignment), bam_alignment_core.l_qseq);
	  qual[bam_alignment_core.l_qseq] = '\0';
	  qual_to_ascii( qual);
	  
	  for( i = 0; i < strlen( sequence); i++)
	    {
	      next_char = base_as_char( bam_seqi( sequence, i));
	      read[i] = next_char;
	    }
	  read[i] = '\0';
	  strcpy(read2, read);
	  //fprintf(stdout, "%s\n", bam_get_qname(bam_alignment));
	  //fprintf(stdout, "%s\n", read);
	  //fprintf(stdout, "%s\n",  qual);
	  //fprintf(stdout, "n_cigar: %d\n", n_cigar);

	  clipped=0;
	  cigar_add_len = 0;

	  fprintf(stdout, "\nCIGAR: ");
	  for (i=0; i<n_cigar; i++){
	    if (bam_cigar_opchr(cigar[i]) == 'S' || bam_cigar_opchr(cigar[i]) == 'H'){
	      hts_itr_destroy(iter);
	      clipped=1;
	      break;
	    }
	    else if (bam_cigar_opchr(cigar[i]) == 'D')
	      cigar_add_len += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_opchr(cigar[i]) == 'I')
	      cigar_add_len -= bam_cigar_oplen(cigar[i]);
	    //fprintf(stdout, "%d\t%c\t%d\t", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]), bam_cigar_type(cigar[i]));
	    fprintf(stdout, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
	  }
	
	  fprintf(stdout, "\n");
	  
	  if (clipped) continue;

	  strcpy(md, bam_aux_get(bam_alignment, "MD"));
	  fprintf(stdout, "MD: %s\n", md);

	  map_tid = bam_alignment_core.tid;
	  map_loc = bam_alignment_core.pos;
	  strcpy(map_chr, bam_header->target_name[map_tid]);
	  
	  ref_seq = faidx_fetch_seq(params->ref_fai, map_chr, map_loc, map_loc-1+strlen(read)+cigar_add_len, &loc_len);
	  strcpy(ref_seq2, ref_seq);
	  
	  //fprintf(stdout, "%s\t%d\t%d\n%s\n%s\n", map_chr, map_loc, loc_len, read, ref_seq);
	  //fprintf(stdout, "pre\n%s\n%s\n", read, ref_seq);

	  //	  return_value = check_map(read, ref_seq, cigar, md);
	  
	  apply_cigar_md(ref_seq, read, md+1, n_cigar, cigar);

	  //	  fprintf(stdout, "\npos\n%s\n%s\n", read, ref_seq);

	  if (strcmp(read, ref_seq)){
	    fprintf(stdout, "\npre\n%s\n%s\n", read2, ref_seq2);
	    fprintf(stdout, "\npos\n%s\n%s\n", read, ref_seq);
	    return;
	  }

	  free(ref_seq);

	 
	  hts_itr_destroy(iter);
		  
	  j++;
	}
}


void get_sample_name( bam_info* in_bam, char* header_text)
{
	/* Delimit the BAM header text with tabs and newlines */

        char *tmp_header = NULL;
	set_str( &( tmp_header), header_text);
	char* p = strtok( tmp_header, "\t\n");
	char sample_name_buffer[1024];

	while( p != NULL)
	{
		/* If the current token has "SM" as the first two characters,
			we have found our Sample Name */
		if( p[0] == 'S' && p[1] == 'M')
		{
			/* Get the Sample Name */
			strncpy( sample_name_buffer, p + 3, strlen( p) - 3);

			/* Add the NULL terminator */
			sample_name_buffer[strlen( p) - 3] = '\0';

			/* Exit loop */
			break;
		}
		p = strtok( NULL, "\t\n");
	}

	set_str( &( in_bam->sample_name), sample_name_buffer);
	free( tmp_header);
}


		
void print_bam( bam_info* in_bam)
{
	/*printf( "Number of Chromosomes: %d\n", in_bam->num_chrom);

	int i;
	for( i = 0; i < in_bam->num_chrom; i++)
	{
		printf( "Chromosome Name: %s\n", ( in_bam->chrom_names)[i]);
		printf( "Length of the Chromosome: %d\n", ( in_bam->chrom_lengths)[i]);
	}*/
}

