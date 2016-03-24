#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "common.h"

// Track memory usage
long long memUsage = 0;

void init_params( parameters** params)
{
	int i;
	/* initialize parameters */
	*params = ( parameters*) malloc( sizeof( parameters));
	( *params)->ref_genome = NULL;
	( *params)->bam_file = NULL;
	( *params)->maps_to_test = 100;
	( *params)->threads = 1;
}

void load_chrom_properties(parameters* params)
{
	FILE* fai_file;
	int ln_count=0,i,c;
	char filename[255],first_arg[255],sec_arg[255];
	int return_value;

	sprintf( filename, "%s.fai", ( params)->ref_genome);
	fai_file= safe_fopen( filename, "r");


	//count the number of chromosomes by counting the non-empty lines in the .fai file
	do{
		c=fgetc(fai_file);
		if(c=='\n') ln_count++;

	}while (c!=EOF);
	params->num_chrom=ln_count;

	/* Reset the file pointer to the start of the file */
	rewind(fai_file);

	//We need the names and lengths of each chromosome
	params->chrom_lengths = ( int*) malloc( params->num_chrom * sizeof( int));
	params->chrom_names = ( char**) malloc( params->num_chrom * sizeof( char*));
	for( i = 0; i < params->num_chrom; i++)
		{
		  return_value = fscanf( fai_file, "%[^\t]\t%[^\t]%*[^\n]",first_arg,sec_arg);
		  params->chrom_names[i]=NULL;
		  set_str( &(params->chrom_names[i]), first_arg);
		  params->chrom_lengths[i]=atoi(sec_arg);
		}
	fclose(fai_file);
}

void print_params( parameters* params)
{
	int i;

	printf( "BAM input: %s\n", params->bam_file);
	
	printf( "ref_genome: %s\n", params->ref_genome);

}

void print_error( char* msg)
{
	/* print error message and exit */
	fprintf( stderr, "\n%s\n", msg);
	fprintf( stderr, "Invoke parameter -h for help.\n");
	exit( EXIT_COMMON);
}


FILE* safe_fopen( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];

	file = fopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
		
	}
	return file;
}

gzFile safe_fopen_gz( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
        gzFile file;
	char err[500];

	file = gzopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);		
	}
	return file;
}

htsFile* safe_hts_open( char* path, char* mode)
{
	htsFile* bam_file;
	char err[500];
	
	bam_file = hts_open( path, mode);
	if( !bam_file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}

	return bam_file;
}

int is_proper( int flag)
{
        if ( (flag & BAM_FPAIRED) != 0 && (flag & BAM_FSECONDARY) == 0 && (flag & BAM_FSUPPLEMENTARY) == 0 && (flag & BAM_FDUP) == 0 && (flag & BAM_FQCFAIL) == 0)
	         return 1;

	return 0;
}

int is_concordant( bam1_core_t bam_alignment_core, int min, int max)
{
	int flag = bam_alignment_core.flag;

	if( ( flag & BAM_FPAIRED) == 0) 
	{
	        /* Read is single-end. Skip this by calling it concordant */
	        return 1;
	}

	if( ( flag & BAM_FPROPER_PAIR) == 0) 
	{
		/* Not proper pair */
		return 0;
	}

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return 0;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return 0;
	}

	if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
	{
		/* -- orientation = inversion */
		return 0;
	}

	if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
	{
		/* ++ orientation = inversion */
		return 0;
	}

	if( bam_alignment_core.tid != bam_alignment_core.mtid) 
	{
		/* On different chromosomes */
		return 0;
	}

	if( bam_alignment_core.pos <= bam_alignment_core.mpos) // c.a.
	{
		/* Read is placed BEFORE its mate */
		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			/* -+ orientation = tandem duplication */
			return 0;
		}
	}
	else
	{
		/* Read is placed AFTER its mate */
		if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			/* +- orientation = tandem duplication */
			return 0;
		}
	}

	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs(bam_alignment_core.isize) < min || abs(bam_alignment_core.isize) > max) // c.a.
	{
		/* Deletion or Insertion */
		return 0;
	}

	/* All passed. Read is concordant */
	return 1;
}

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
}

/* Return the complement of a base */
char complement_char( char base)
{
	switch( base)
	{
		case 'A': 
			return 'T';
			break;
		case 'C': 
			return 'G';
			break;
		case 'G': 
			return 'C';
			break;
		case 'T': 
			return 'A';
			break;
		default: 
			return 'N';
			break;
	}
	return 'X';
}

/* Add 33 to the integer value of the qual characters to convert them to ASCII */
void qual_to_ascii( char* qual)
{
	int i;
	for( i = 0; i < strlen( qual); i++)
	{
		qual[i] = qual[i] + 33;
	}
}

/* Even safer than strncpy as it dynamically allocates space for the string if
 there hasn't been already */
void set_str( char** target, char* source)
{
	if( *target != NULL)
	{
		free( ( *target));
	}
	
	if (source != NULL)
	{
	        ( *target) = ( char*) malloc( sizeof( char) * ( strlen( source) + 1));
		strncpy( ( *target), source, ( strlen( source) + 1));
	}
	else
	{
	        ( *target) = NULL;
	}
}


/* Reverse a given string */
void reverse_string( char* str)
{
	int i;
	char swap;
	int len = strlen( str);

	for( i = 0; i < len / 2; i++)
	{
		swap = str[i];
		str[i] = str[len - i - 1];
		str[len - i - 1] = swap;
	}
}

int compare_size_int( const void* p, const void* q)
{
    int i = *( const int*) p;
    int j = *( const int*) q;

	if( i < j)
	{
		return -1;
	}
	else if( i == j)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void* getMem( size_t size)
{
	void* ret;

	ret = malloc( size);
	if( ret == NULL)
	{
		fprintf( stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory.\n", getMemUsage(), ( float) ( size / 1048576.0));
		exit( 0);
	}

	memUsage = memUsage + size;
	return ret;
}

void freeMem( void* ptr, size_t size)
{
	memUsage = memUsage - size;
	free( ptr);
}

double getMemUsage()
{
	return memUsage / 1048576.0;
}

void del_char(char *ref, int start, int len){
  int ref_len = strlen(ref);
  int i;
  for (i=start;i<start+len && i<ref_len;i++)
    ref[i]=ref[i+len];
  while (i < ref_len - len){
    ref[i]=ref[i+len];
    i++;
  }
  ref[i]=0;
}

void ins_char(char *ref, char *read, int start, int len){
  int ref_len = strlen(ref);
  int i;

  for (i=ref_len+len; i>=start; i--)
    ref[i]=ref[i-len];

  ref[ref_len+len]=0;
  ref_len = ref_len+len;

  for (i=start;i<start+len && i<ref_len;i++){   
    ref[i]= read[i]; //'.';
  }

}

void apply_cigar_md(char *ref, char *read, char *md, int n_cigar, const uint32_t *cigar){
  //void applymd(char *ref, char *md){
  /* assuming Z is removed */
  int i,j,k;
  char buf[1000];
  int oplen;
  int refptr=0;
  int edit_loc;
  int delcnt;
  int thisdel;
  int inserted;
  int skipk;

  /* TODO: when I , add to refptr*/

  edit_loc = 0;
  for (i=0; i<n_cigar; i++){
    if (bam_cigar_opchr(cigar[i]) == 'M')
      edit_loc += bam_cigar_oplen(cigar[i]);
    else if (bam_cigar_opchr(cigar[i]) == 'D'){
      del_char(ref, edit_loc, bam_cigar_oplen(cigar[i]));
      //edit_loc -= bam_cigar_oplen(cigar[i]);
    }
    else if (bam_cigar_opchr(cigar[i]) == 'I'){
      ins_char(ref, read, edit_loc, bam_cigar_oplen(cigar[i]));
      edit_loc += bam_cigar_oplen(cigar[i]);
    }
    //fprintf(stdout, "%d\t%c\t", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(c igar[i]));
    //fprintf(stdout, "%d\t%d\t%c\t%d\t", bam_cigar_op(cigar[i]), bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]), bam_cigar_type(cigar[i]));
  }
  


  j=0;
  buf[0]=0;
  i=0;
  delcnt=1;
  inserted=0;
  skipk=0;
  //printf("md: %s\n", md);
  while (i<strlen(md)){
    //printf("md[%d]: %c\n", i, md[i]);
    if (isdigit(md[i])){
      //printf("digit %c\n", md[i]);
      buf[j++]=md[i];
    }
    else {
      buf[j]=0;
      j=0;
      oplen=atoi(buf);
      refptr+=oplen;
      printf("buf: %s : %d\n", buf, oplen);
    }

    if (md[i]=='^'){ // del. skip
      //while (isalpha(md[i])) {i++; refptr++;}      
      thisdel = 0;
      for (k=0; k<n_cigar; k++){
	if (bam_cigar_opchr(cigar[k]) == 'D'){
	  thisdel++;
	  if (thisdel == delcnt){
	    //while (isalpha(md[i])) {i++; refptr++;}
	    //while (isalpha(md[i])) {i++; }
	    
	    i+=bam_cigar_oplen(cigar[k]); 
	    //refptr+=bam_cigar_oplen(cigar[k]); 
	  }
	  delcnt++;
	}
      }
    }
    else if (isalpha(md[i])){
      inserted = 0;
      edit_loc = 0;
      for (k=skipk; k<n_cigar; k++){
	if (bam_cigar_opchr(cigar[k]) == 'I' && edit_loc < refptr+bam_cigar_oplen(cigar[k])){
	  inserted += bam_cigar_oplen(cigar[k]);
	  skipk=k+1;
	}
	edit_loc += bam_cigar_oplen(cigar[k]);
      }
      refptr += inserted;
      printf("changing %d:%c %d:%c\n", refptr, read[refptr], i, md[i]);
      read[refptr++]=md[i];
    }
    i++; 
  }
}
