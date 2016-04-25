#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include "verifybam.h"
#include "cmdline.h"


int parse_command_line( int argc, char** argv, parameters* params)
{
	int index;
	int o;
	static struct option long_options[] = 
	{
		{"input"  , required_argument,   0, 'i'},
		{"ref"    , required_argument,   0, 'f'},
		{"maps"   , required_argument,   0, 'm'},
		{"help"   , no_argument,         0, 'h'},
		{"threads", required_argument,   0, 't'},	
		{"version", no_argument,         0, 'v'},
		{0        , 0,                   0,  0 }
	};
  
	if( argc == 1)
	{
		print_help();
		return 0;
	}
  
	while( ( o = getopt_long( argc, argv, "hv:i:f:t:m:", long_options, &index)) != -1)
	{
		switch( o)
		{
			case 'i':
			  set_str( &( params->bam_file), optarg);
			break;
	  
			case 'f':
				set_str( &( params->ref_genome), optarg);
			break;

			case 'm':
			        params->maps_to_test = atoi (optarg);
			break;
			
			case 't':
				params->threads = atoi( optarg);
			break;

			case 'h':
				print_help();
				return 0;
			break;

			case 'v':
				fprintf( stderr, "\nVERIFYBAM: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
				fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", VERIFYBAM_VERSION, VERIFYBAM_UPDATE, BUILD_DATE);
				fprintf( stderr, "It is bigger on the inside!\n\n");
				return 0;
			break; 
		}
	}
  
	/* TODO: check parameter validity */
  

	if( params->bam_file == NULL)
	{
	        fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter an input BAM file using the --input option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --ref   is invoked */
	if( params->ref_genome == NULL)
	{
		fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter reference genome file (FASTA) using the --ref option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if threads>0 */
	if( params->threads <= 0)
	{
		fprintf( stderr, "[VERIFYBAM CMDLINE WARNING] Invalid number of threads was entered (%d). Resetted to 1.\n", params->threads);
		params->threads = 1;
	}

	return 1;

}

void print_help( void)
{  
	fprintf( stdout, "\nVERIFYBAM: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", VERIFYBAM_VERSION, VERIFYBAM_UPDATE, BUILD_DATE);	
	fprintf( stdout, "\t--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.\n");
	fprintf( stdout, "\t--input [BAM files]        : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.\n");
	fprintf( stdout, "\t--out   [output prefix]    : Prefix for the output file names.\n");
	fprintf( stdout, "\t--ref   [reference genome] : Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--gaps  [gaps file]        : Assembly gap coordinates in BED3 format.\n");
	fprintf( stdout, "\t--dups  [dups file]        : Segmental duplication coordinates in BED3 format.\n");
	fprintf( stdout, "\t--reps  [reps file]        : RepeatMasker annotation coordinates in BED6 format. See manual for details.\n");
	fprintf( stdout, "\t--mei   [\"Alu:L1Hs:SVA\"]   : List of mobile element names.\n");
	/*
	fprintf( stdout, "\t--xx                       : Sample is male.\n");
	fprintf( stdout, "\t--xy                       : Sample is female.\n");
	*/
	fprintf( stdout, "\t--vh                       : Run VariationHunter/CommonLAW (read pair + read depth).\n");
	/* not  yet implemented, hide the parameters
	fprintf( stdout, "\t--ns                       : Run NovelSeq (read pair + assembly).\n");
	fprintf( stdout, "\t--sr                       : Run SPLITREAD (split read).\n");
	fprintf( stdout, "\t--all                      : Run all three algorithms above [DEFAULT].\n");
	*/
	fprintf( stdout, "\t--skip-fastq               : Skip FASTQ dump for discordants. Use this only if you are regenerating the calls.\n");
	fprintf( stdout, "\t--skip-sort                : Skip FASTQ sort for discordants. Use this only if you are regenerating the calls.\n");
	fprintf( stdout, "\t--skip-remap               : Skip FASTQ remapping for discordants. Use this only if you are regenerating the calls.\n");
	fprintf( stdout, "\t--version                  : Print version and exit.\n");
	fprintf( stdout, "\t--help                     : Print this help screen and exit.\n\n");
	fprintf( stdout, "It is bigger on the inside!\n\n");
}

