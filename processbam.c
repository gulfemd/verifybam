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
	bam1_core_t bam_alignment_core;
	bam1_t*	bam_alignment;
	int** fragment_size;
	int** second_pass_fragments;
	int* fragments_sampled;
	int* second_test_pass;
	int* fragment_size_total;
	float* variance;
	int diff;
	int return_value;
	int i;
	int j;

	fprintf( stderr, "Processing BAM file %s.\n", path);


	/* Open the BAM file for reading. htslib automatically detects the format
		of the file, so appending "b" after "r" in mode is redundant. */
	bam_file = safe_hts_open( path, "r");

	/* Read in BAM header information */
	bam_header = bam_hdr_read( ( bam_file->fp).bgzf);


	get_sample_name( in_bam, bam_header->text);


	/* Initial read */	
	bam_alignment = bam_init1();
	return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);


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

