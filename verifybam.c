#include "verifybam.h"

int main( int argc, char** argv)
{
	bam_info* in_bam;
	parameters* params;
	int return_value;
	char username[MAX_SEQ];
	int i;
	int j;
	FILE *refFile;

	/* Set program parameters */
	init_params( &params);
	
	/* Seed random number generator */
	srand(time(NULL));

	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params);
	if( return_value == 0)
	{
		exit( EXIT_SUCCESS);
	}
	else if( return_value != 1)
	{
		exit( return_value);
	}

	if ( VERIFYBAM_DEBUG)
        {
 	        print_params( params);
	}

	load_chrom_properties(params);

	/* read reference genome */
	/*
	refFile = safe_fopen(params->ref_genome, "r");
	readSingleFasta(refFile);
	*/
	
	params->ref_fai = fai_load(params->ref_genome);
	
	/* Read BAM files and calculate the median/avg/std of fragment sizes per library */
	in_bam = ( bam_info*) malloc( sizeof( bam_info));
	in_bam->sample_name = NULL;
	load_bam( in_bam, params->bam_file);
	
	/* Initial read */	
	read_alignment(in_bam, params);


	/* BAM is loaded, min/max/avg/std are calculated. Now, extract FASTQs of discordants, OEAs, and orphans */


	return EXIT_SUCCESS;
}
