
/* ******************* */
/* SeqCount.c:  */

/* A program that takes a FASTA file of sequences and produces a file */
/* summarizing the counts for each codon or for each nucleotide.      */

/* Author: John Novembre (novembre@berkeley.edu)             */
/* Date: Feb 28th, 2006                                               */

/* ******************* */


#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#define MAXLINE 500
struct sequence {
  char name[MAXLINE];
  char *seq;
  int cod_count[64];
  int codonlength;
  double f_A,f_T,f_C,f_G;
};

void read_file(char *inputfile, int *num_seq,struct sequence *seqs);

int main(int argc, char **argv) {
  

  FILE *out, *outfreq,*outgc,*outgccnt;
  char outfilename[MAXLINE];
  int num_seq;
  struct sequence *seqs;
  struct sequence totals;
  char codon[64][4];
  int i,j,k,l;
  int ctr;
  int n_cod;
  int file_format;
  int length; // in nucleotides
  int total_length;
  int countwhat;  /* 1=codons 2=nucleotides */
  char inputfile[40];
  char input[40];
  char ch;
  if(argc==4){
    if(!strcmp(argv[1],"-c")){
      countwhat=1;
    }else if(!strcmp(argv[1],"-n")){
      countwhat=2;
    }else{
      fprintf(stderr,"%s option not recognized.\n",argv[1]);
      fprintf(stderr,"Standard options are: -c for counting codons and -n for counting nucleotides.\n",argv[1]);
      exit(-1);
    }
    sprintf(inputfile,"%s",argv[2]);
    num_seq=atoi(argv[3]);
    if(num_seq<1) {
      fprintf(stderr,"Error reading number of sequences.  Read in %d.\n",num_seq);
      exit(-1);
    }
  }else{
    do{
      printf("Enter the FASTA-format input filename:\n");
      scanf("%s",input);	  
      if(!fopen(input,"r")){
	printf("Cannot find the file %s.\n",input);	    
	i=1;
      }else{i=0;}
      strcpy(inputfile,input);
    }while(i);    
    getchar(); /* to clean the newline */
    countwhat=0;
    do{
      printf("Enter 'C' to count codons, 'N' to count nucleotides:\n");
      scanf("%c",&ch);
      if(toupper(ch)=='C') countwhat=1;
      if(toupper(ch)=='N') countwhat=2;
    }while(countwhat==0);
    do{
      printf("Enter the number of sequences:\n");
      scanf("%d",&num_seq);
    }while((num_seq<1));
  }

  seqs=(struct sequence *)calloc((size_t)num_seq,sizeof(struct sequence));
  if(seqs==NULL) {
    fprintf(stderr,"Memory allocation failed.");
    exit(-1);
  }

  read_file(inputfile,&num_seq,seqs);

  if(countwhat==1){
    ctr=0;
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) {
	for(k=0;k<4;k++) {
	  codon[ctr][0]=i+1;
	  codon[ctr][1]=j+1;
	  codon[ctr][2]=k+1;
	  codon[ctr][3]='\0';
	  ctr++;
	}
      } 
    }
    
    
    for(i=0;i<64;i++) {
      for(j=0;j<strlen(codon[i]);j++) {
	switch(codon[i][j]) {
	case 1 : {codon[i][j]='T'; break;}
	case 2 : {codon[i][j]='C'; break;}
	case 3 : {codon[i][j]='A'; break;}
	case 4 : {codon[i][j]='G'; break;}
	default: {codon[i][j]='\0';break;}
	}
      }
    }
    
    printf("Initializing codon count array...\n");
    for(i=0;i<num_seq;i++){
      for(k=0;k<64;k++) {
	seqs[i].cod_count[k]=0;
      }
      seqs[i].f_A=seqs[i].f_T=seqs[i].f_C=seqs[i].f_G=0;
      seqs[i].codonlength=0;
    }
    

    for(i=0;i<64;i++){
    	totals.cod_count[i]=0;
    }

    totals.f_A=totals.f_T=totals.f_C=totals.f_G=0.0;
    totals.codonlength=0;
    
    printf("Counting codons...\n"); 
    for(i=0;i<num_seq;i++){
      
      // Convert to uppercase
      for(j=0;seqs[i].seq[j]!='\0';j++){
	seqs[i].seq[j]=toupper(seqs[i].seq[j]);
      }
      for(j=0;seqs[i].seq[3*j]!='\0';j++) {
	for(k=0;(strncmp(seqs[i].seq+3*j,codon[k],3))&&(k<65);k++){
	  
	};

	seqs[i].cod_count[k]++;
	seqs[i].codonlength++;
	totals.cod_count[k]++;
	totals.codonlength++;
	
      }    
    }
    
    strcpy(outfilename,inputfile);
    strcat(outfilename,".codcnt");

    if((out=fopen(outfilename,"w"))==NULL){
      fprintf(stderr, "Couldn't create output file\n");
      exit(-1);
    }

    strcpy(outfilename,inputfile);
    strcat(outfilename,".codfreq");

    if((outfreq=fopen(outfilename,"w"))==NULL){
      fprintf(stderr, "Couldn't create output file\n");
      exit(-1);
    }


    fprintf(out,"%d\n64\n", num_seq+1);
    fprintf(outfreq,"%d\n64\n", num_seq+1);
    
      
    for(i=0;i<64;i++) {
      fprintf(out,"%s ",codon[i]);
      fprintf(outfreq,"%s ",codon[i]);
    }
    fprintf(out, "\n");
    fprintf(outfreq, "\n");

    for(i=0;i<num_seq;i++) {
      fprintf(out,"%s>",seqs[i].name);
      fprintf(outfreq,"%s> ",seqs[i].name);
      for(j=0;j<64;j++) {
	fprintf(out,"%d ",seqs[i].cod_count[j]);
	fprintf(outfreq,"%.3f ",(double)seqs[i].cod_count[j]/(double)seqs[i].codonlength);
      }
      fprintf(out,"\n");
      fprintf(outfreq,"\n");
    }
    

    fprintf(out,"Totals> ");
    fprintf(outfreq,"Totals> ");
    for(i=0;i<64;i++){
    	fprintf(out,"%d ",totals.cod_count[i]);
    	fprintf(outfreq,"%f ",(double)totals.cod_count[i]/(double)totals.codonlength);
    }
    
    fclose(out);
    fclose(outfreq);
    
  }else if(countwhat==2){
    printf("Finding nucleotide content...\n");
    total_length=0;
    totals.f_A=0.0;
    totals.f_T=0.0;
    totals.f_C=0.0;
    totals.f_G=0.0;
    
    strcpy(outfilename,inputfile);
    strcat(outfilename,".acgtfreq");
   
    if((outgc=fopen(outfilename,"w"))==NULL){ 
      fprintf(stderr, "Couldn't create output nucleotide composition frequency file\n"); 
      exit(-1);  
    } 

    strcpy(outfilename,inputfile);
    strcat(outfilename,".acgtcnt");
   
    if((outgccnt=fopen(outfilename,"w"))==NULL){ 
      fprintf(stderr, "Couldn't create output nucleotide composition count file\n"); 
      exit(-1);  
    } 

    
    fprintf(outgc,"Nucleotide f_A f_C f_G f_T\n");
    fprintf(outgccnt,"Nucleotide n_A n_C n_G n_T n_Total\n");
    
    for(i=0;i<num_seq;i++) {
      length=0;
      for(j=0;seqs[i].seq[j]!='\0';j++,length++) {
	
	switch(toupper(seqs[i].seq[j])) {
	case 'A' : {seqs[i].f_A++; totals.f_A++; break;}
	case 'T' : {seqs[i].f_T++; totals.f_T++; break;}
	case 'C' : {seqs[i].f_C++; totals.f_C++; break;}
	case 'G' : {seqs[i].f_G++; totals.f_G++; break;}
	};
      }
      fprintf(outgccnt,"%s> %.0f %.0f %.0f %.0f %d\n",seqs[i].name,seqs[i].f_A, seqs[i].f_C, seqs[i].f_G, seqs[i].f_T,length);

      seqs[i].f_A/=length;
      seqs[i].f_T/=length;
      seqs[i].f_C/=length;
      seqs[i].f_G/=length;

      fprintf(outgc,"%s> %f %f %f %f\n",seqs[i].name,seqs[i].f_A, seqs[i].f_C, seqs[i].f_G, seqs[i].f_T);

      total_length+=length;  
    }

    fprintf(outgccnt,"Totals> %.0f %.0f %.0f %.0f %d\n",totals.f_A,totals.f_C,totals.f_G,totals.f_T,total_length);

    totals.f_A/=total_length;
    totals.f_C/=total_length;
    totals.f_G/=total_length;
    totals.f_T/=total_length;

    fprintf(outgc,"Totals> %f %f %f %f\n",totals.f_A,totals.f_C,totals.f_G,totals.f_T);
	
    fclose(outgc);    
    fclose(outgccnt);
  }

  free(seqs);
  exit(0);
}


void read_file(char *inputfile, int *num_seq,struct sequence *seqs) {

  
  int ch;
  FILE *in;
  int seqlength;
  long seqstart;
  int i,j;

  if((in =fopen(inputfile,"r"))==NULL){
    fprintf(stderr, "Couldn't open %s \n",inputfile);
    exit(-1);
  }

  printf("Reading input file...\n");
  /* Clear the first '>' */
  ch=getc(in);

  for(i=0;i<*num_seq;i++) {
    
    if(!fgets(seqs[i].name,MAXLINE,in)){
      fprintf(stderr, "Specified number of sequences appears to be too high.  Only reading in %d sequences",i);
      *num_seq=i;
      i=*num_seq;
    };
    if(strchr(seqs[i].name,'\n')){
      *strchr(seqs[i].name,'\n')='\0';
    }

    seqstart=ftell(in);
    /* Find length of sequence */
    for(j=0;((ch=getc(in))!='>')&&(ch!=EOF);j++) {
      if((ch=='\n')||(ch==' ')){ j--; continue;}
    }


    seqlength=j;

    /* Declare memory */
    seqs[i].seq=(char *)calloc((size_t)seqlength,sizeof(char));
    if(seqs[i].seq==NULL){
      fprintf(stderr,"Memory allocation failure.  Sequence %s may be too long at %d bases.",seqs[i].name,seqlength);
    }
    
    fseek(in,seqstart,SEEK_SET);


    /* Read in sequence     */
    for(j=0;((ch=getc(in))!='>')&&(ch!=EOF);j++) {
      if((ch=='\n')||(ch==' ')){ j--; continue;}
      else seqs[i].seq[j]=ch;
    }
    seqs[i].seq[j]='\0';
  }
  
}



