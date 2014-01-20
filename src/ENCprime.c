/* ********************/
/* ENCprime.c:  */

/* A program to calculate a handful of codon usage bias statistics.  */
/* These include: Nc, Nc', Akashi's scaled chi-square, a sum of  */
/* chi-squares across codons, and Karlin and Mrazek's B*(a).  */

/* Author: John Novembre (novembre@berkeley.edu)    */
/* Date: July 17th, 2002   */

/* ********************/


/* Things that would be nice to do: */
/*   --Output Nc,Ncp statistics for each amino acid */


#include <string.h>
#include <ctype.h>  /* for toupper()...  */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __MWERKS__
	#ifdef macintosh
		#include <console.h>
	#endif
#endif


#define MAX_CODONS 64
#define MAX_AA 45
#define MAX_REDUN 8
#define MAX_NAME_LENGTH 500
#define MAX_SEQ_LENGTH 500000
#define MIN_AA_N 20
#define NUM_CODES 25
#define MAX_CODE_LENGTH 800
#define LINE_LENGTH 200
/*  non-taxa specific information about each codon   */
struct codon{

  /*  read in from input data files  */
  char seq[4];
  struct a_acid *aa;

  /*  computed by the program  */
  int gc_bp; /*  number of g/c basepairs in the codon  */
  int var_sites[3]; /*  an array marking the variable sites as 1,   */
                    /*  non-varying sites as 0.  */
};

/*  non-taxa specific info about each amino acid  */
struct a_acid{
  char abbrev;
  int redun;    /*  redundancy family  */
  struct codon *cods[MAX_REDUN]; /*  array of the codons for one amino acid  */
};

/*  taxa specific codon usage data and statistics  */
struct codon_datum {
  int n; /*  number of times observed in sample  */

  struct codon *cod; /*  pointer to the codon info  */

  double ex_num;  /*  numerator of expected frequency  */
  double ex_denom; /*  denominator of expected frequency  */
  double ex; /*  expected frequency  */

  double p; /*  actual frequency  */
  double bias; /*  bias value  */

};

/*  taxa specific stats calculated for each amino acid  */
struct aa_stat {
  struct a_acid *aa;
  struct codon_datum *cods[MAX_REDUN];
  int n; /*  number of times the aa is observed  */
  double chisq; /*  chi-sqaure values  */
  double p; /*  p values for each model  */
  int m; /*  true if null  model is rejected, false if not rejected...  */
  int bad_ex; /*  if true, expected proportion is too small   */
  double F; /*  "homozygosity"  of codons for this amino acid (for ENC calc)  */
  double F_p; /*  "homozyosity" calculated from chi-square  */
  double Bg;  /*  Mrazek and Karlin measure  */
};

struct sequence{

  /*  read in from input data files  */
  char name[MAX_NAME_LENGTH];
  char seq[MAX_SEQ_LENGTH];
  int n_cod_obs;
  double f_A,f_T,f_C,f_G;
  double f_A3,f_T3,f_C3,f_G3;

  struct codon_datum cod_data[MAX_CODONS];  /* array of codon data  */

  /*  computed by program  */

  struct aa_stat a_stat[MAX_AA];  /* array of stats for each amino acid   */

  double sumchi;
  int df;
  double p;
  double scaledchi; /*  Akashi's scaled chi-square  */
  double N_c;  /*  Wright's effective number of codons for the sequence  */
  double N_c_p; /*  Nc prime  */

  double B_KM; /*  Mrazek and Karlin measure  */
};


int verbosity=2;
struct sequence *seq;
struct codon cod[MAX_CODONS];
struct a_acid aa[MAX_AA];
int n_seq=0;
int n_cod=0;
int n_aa=0;
int n_a_data=0;
double alpha=0.05;
double eps=1e-3;
int quick_mode=0;
int expected_given=0;
char outputfile[40];


/*  From Eric Anderson  */
double IncGammaCF(double a, double x);
double IncGammaSer(double a, double x);
double GammaP(double a, double x);
double LogGamma(double x);

char **codes;

void init_codes();
void get_code(char *filename);
void get_codon_count_data(char *filename);
void get_acgt(char *filename);
void calc_cod_data();
void calc_F();
void test_null();
void calc_n_aa();
void calc_sumb();
void calc_sumchi();
void calc_Nc();
void calc_p_gc_pref();
void output2screen();
void output2file(char *outputfile);

void init_codes();
void get_code(char *filename);
void get_codon_count_data(char *filename);
void get_acgt(char *filename);
void calc_cod_data();
void calc_F();
void test_null();
void calc_n_aa();
void calc_sumb();
void calc_sumchi();
void calc_Nc();
void calc_p_gc_pref();
void output2screen();
void output2file(char *outputfile);


int main (int argc, char **argv){

  /*  Game plan:  */

  /*  get genetic code  */
  /*  read in codon count */
  /*  read in nucleotide comp (or generate it)  */
  /*  calc cod_data  */
  /*  calc aa_data (X^2, F, F_p)  */
  /*  calc Nc  */
  /*  output2 file  */


  char codefile[LINE_LENGTH];
  char datafile[LINE_LENGTH];
  char acgtfile[LINE_LENGTH];
  char input[LINE_LENGTH];
  char junk[LINE_LENGTH];
  int i,j,k;
  char ch;
  FILE *in;
  codes=(char**)calloc((size_t)NUM_CODES,sizeof(char*));
  for(i=0;i<NUM_CODES;i++) {
  	codes[i]=(char*)calloc((size_t)MAX_CODE_LENGTH,sizeof(char));
  }

  /* Read in the defaults */
  if(!(in=fopen("ENCprime.defaults","r"))){
    printf("Could not find ENCprime.defaults file.\n");
    sprintf(datafile,"../data/s_pombe.dna.cnt");
    sprintf(acgtfile,"../data/s_pombe.dna.expected");
    sprintf(codefile,"1");
    sprintf(outputfile,"example");
    verbosity=2;
    quick_mode=0;
  } else{
    fgets(input,LINE_LENGTH,in);
    strtok(input,":");
    sscanf(strtok(NULL,":"),"%s",datafile);
    fgets(input,LINE_LENGTH,in);
    strtok(input,":");
    sscanf(strtok(NULL,":"),"%s",acgtfile);
    fgets(input,LINE_LENGTH,in);
    strtok(input,":");
    sscanf(strtok(NULL,":"),"%s",codefile);
    fgets(input,LINE_LENGTH,in);
    strtok(input,":");
    sscanf(strtok(NULL,":"),"%s",outputfile);
    fgets(input,LINE_LENGTH,in);
    strtok(input,":");
    verbosity=atoi(strtok(NULL,":"));
    fgets(input,LINE_LENGTH,in);
    strtok(input,":");
    quick_mode=atoi(strtok(NULL,":"));
  }

#ifdef __MWERKS__
	#ifdef macintosh
		argc = ccommand(&argv);
	#endif
#endif


  printf("--ENCprime Program--\n");

  if(argc==1){

    while(1){
      /*printf("Menu-driven action baby.\n");*/
      printf("====Main Menu====\n");
      printf("(1) Codon counts file: \t\t\t%s\n",datafile);
      printf("(2) Nucleotide composition file: \t%s\n",acgtfile);
      printf("(3) Genetic Code: \t\t\t%s \t(Number between 1 and %d indicates a Genbank code, otherwise input filename).\n",codefile,NUM_CODES);
      printf("(4) Output file name:\t\t\t%s\n",outputfile);
      printf("(5) Verbosity: \t\t\t\t%d \t(0=silent, 1=show input, 2=show input and pause)\n",verbosity);
      printf("(6) Data Explorer: \t\t\t%d \t(0=Enter data explorer, 1=Exit immediately after analysis)\n",quick_mode);

      do{
	printf("\nType 1-6 and return to make a change to the options...\n");
	printf("Type 7 and return to begin analysis, or type 8 to exit.\n");
	do
	  ch=getchar();
	while(ch=='\n');
	j=atoi(&ch);
      }while((j<1)||(j>8));
      if(j==7) break;
      switch(j){
      case 8: { printf("Exiting ENCprime...\n");exit(0);}
      case 1: {
	do{
	  printf("Enter the codon counts filename or 0 to cancel:\n");
	  scanf("%s",input);
	  if((strlen(input)==1)&&(atoi(input)==0)) break;
	  if(!fopen(input,"r")){
	    printf("Cannot find the file %s.\n",input);
	    i=1;
	  }else{i=0;sprintf(datafile,"%s",input);}

	}while(i);

	break;
      }
      case 2: {
	do{
	  printf("Enter the nucleotide composition filename or 0 to cancel:\n");
	  scanf("%s",input);
	  if((strlen(input)==1)&&(atoi(input)==0)) break;
	  if(!fopen(input,"r")){
	    printf("Cannot find the file %s.\n",input);
	    i=1;
	  }else{i=0;sprintf(acgtfile,"%s",input);}

	}while(i);

	break;
      }
      case 3: {
	do{
	  printf("Enter the genetic code's Genbank id number or the filename or 0 to cancel:\n");
	  scanf("%s",input);


	  k=strlen(input);
	  if(k<=2){
	    i=atoi(input);
	    if(i==0) break;
	    if((i<1)||(i>NUM_CODES)||(i==7)||(i==8)||((i>16)&&(i<21))){
	      printf("No code %d.  Try again.\n",i);
	      i=1;
	    }  else {
	      sprintf(codefile,"%s",input);
	      i=0;
	    }
	  } else if(!fopen(input,"r")){
	    printf("Cannot find the file %s.\n",input);
	    i=1;
	  } else{i=0;sprintf(codefile,"%s",input);}

	}while(i);

	break;
      }
      case 4: {
	printf("Enter the output file's prefix:\n");
	scanf("%s",input);
	if((strlen(input)==1)&&(atoi(input)==0)) break;
	sprintf(outputfile,"%s",input);
	break;
      }
      case 5:{
	do{
	  printf("Enter the verbosity (0,1,2):\n");
	  scanf("%d",&verbosity);
	}while((verbosity<0)||(verbosity>2));
	break;
      }
      case 6:{
	do{
	  printf("Skip the data explorer? (1=yes, 0=no)\n");
	  scanf("%d",&quick_mode);
	}while((quick_mode<0)||(quick_mode>1));
	break;
      }

      }

    }
    getchar();
  } else if(argc<6) {
    printf("Usage: %s <cnt file> <nuc comp file> <gen code file> <output file> <verbosity> <-q for quick>\n",argv[0]);
    exit(-1);
  } else {
    sprintf(datafile,"%s",argv[1]);
    sprintf(acgtfile,"%s",argv[2]);
    sprintf(codefile,"%s",argv[3]);
    sprintf(outputfile,"%s",argv[4]);
    verbosity=atoi(argv[5]);
    if(argc==7) {
      if(!strcmp(argv[6],"-q")) {
	quick_mode=1;
      }
    }
  }

  printf("Codon usage file: %s\n",datafile);
  printf("Null nucleotide usage file: %s\n",acgtfile);
  printf("Genetic code file: %s\n",codefile);
  printf("Output data file: %s\n",outputfile);
  printf("Verbosity: %d\n",verbosity);

  init_codes();
  get_code(codefile);
  get_codon_count_data(datafile);
  get_acgt(acgtfile);

  calc_cod_data();

  calc_n_aa();
  test_null();
  calc_F();
  calc_Nc();  /*  now calculate the N_c for the taxa  */
  calc_sumchi();

  output2screen();

  output2file(outputfile);
  free(seq);
  for(i=0;i<NUM_CODES;i++) {
  	free(codes[i]);
  }
  free(codes);
}


void get_code(char *filename) {
  int ch, pos;
  int codeID;
  FILE *in;
  int i, j, k,l;

  /* Init aa abbreviations*/
  for(i=0;i<MAX_AA;i++){
    aa[i].abbrev='\0';
  }

  /* If filename is 1-11? choose the appropriate genetic code */

  i=strlen(filename);
  if(i<=2){
    codeID=atoi(filename);
    if((codeID<1)||(codeID>NUM_CODES)||(codeID==7)||(codeID==8)||((codeID>16)&&(codeID<21))){
      fprintf(stderr,"No code %d ",codeID);
      exit(-1);
    }
  } else{
    /* Open the file and copy it to a string */
    if(!(in=fopen(filename,"r"))){
      fprintf(stderr,"Error opening genetic code file %s\n", filename);
      exit(-1);
    }

    printf("Reading in genetic code...\n");

    codeID=0;
    i=0;
    do{
      i++;
      codes[codeID][i]=fgetc(in);
    }while((ch!=EOF)&&(i<MAX_CODE_LENGTH));

    printf("Closing genetic code data file...\n");
    fclose(in);

  }

  pos=0;
  /*  read data  */
  /*  skip to second line  */
  do pos++; while(codes[codeID][pos]!='=');pos++;
  do pos++; while(codes[codeID][pos]==' ');

  for(k=0;k<64;k++){
    ch=codes[codeID][pos];
    i=0;
    /* advance until you find a spot that either matchs c or is null  */
    while((aa[i].abbrev!=ch)&&(aa[i].abbrev!='\0')){
      i++;
    };

    n_aa=i+1;
    if(n_aa>MAX_AA) {
      fprintf(stderr,"More amino acids (%d) in genetic code than program allows(%d).\n",n_aa,MAX_AA);
      exit(-1);
    }

    aa[i].abbrev=ch;
    cod[k].aa=&aa[i];
    pos++;
  }

  /*  skip to second line  */
  do pos++; while(codes[codeID][pos]!='=');

  for (i=0;i<3;i++){
    /*  skip to next line   */
    do pos++; while(codes[codeID][pos]!='=');pos++;
    do pos++; while(codes[codeID][pos]==' ');
    /*  read data  */
    for(k=0;k<64;k++){
      cod[k].seq[i]=codes[codeID][pos];
      pos++;

    }
  }

  /*  add null character  */
  for(k=0;k<64;k++)
    cod[k].seq[3]='\0';

  /*  calculate redundancy of each amino acid  */
  for(i=0;i<n_aa;i++){
    aa[i].redun=0;
    l=0;
    for(k=0;k<64;k++){
      if(cod[k].aa==&aa[i]) {
	aa[i].redun++;
	aa[i].cods[l]=&cod[k];
	l++;
      }
    }
  }

  /*  Calculate variable sites of each codon   */
  for(k=0;k<MAX_CODONS;k++){

    /* initialize vars  */
    for(i=0;i<3;i++){
      cod[k].var_sites[i]=0;
    };

    /*  calculate vars  */
    for(i=0;i<MAX_CODONS;i++){

      /*  if same amino acid  */
      /*  1.  up redundancy counter  */
      /*  2.  run through sites on sequence and mark variable ones  */
      if(cod[k].aa->abbrev==cod[i].aa->abbrev){
	if(strcmp(cod[k].seq,cod[i].seq)) {
	  for(j=0;j<3;j++){
	    if(toupper(cod[k].seq[j])!=toupper(cod[i].seq[j])) {
	      cod[k].var_sites[j]=1;
	    }
	  }
	}
      }
    }
  }

  if(verbosity>0) {
    printf("\n----Genetic code----\n");
    for(k=0;k<64;k++){
      printf("%s %c %d %d %d %d / ", cod[k].seq, cod[k].aa->abbrev,cod[k].aa->redun, cod[k].var_sites[0], cod[k].var_sites[1], cod[k].var_sites[2]);
      if(((k+1)%5)==0) printf("\n");
    }
    printf("\n----------------------------------------\n");
    if(verbosity>1) {
      printf("Press enter to continue...\n");
      do
	ch=getchar();
      while(ch!='\n');
    }
  }

}


void get_codon_count_data(char *filename){

  FILE *in;
  char cod_seq[4];
  char ch;
  int i, j,k,l;
  if(!(in=fopen(filename,"r"))){
    printf("Error opening codon usage data file %s\n", filename);
    ch=getchar();
    exit(-1);
  }

  /*  Read number of taxa and number of codons from file  */
  printf("Reading codon usage data file... \n");
  fscanf(in, "%d %d", &n_seq,&n_cod);

  /* Declare memory for n_seq sequences */
  seq=(struct sequence *) calloc((size_t)n_seq,sizeof(struct sequence));
  if(seq==NULL) {
    fprintf(stderr,"Memory allocation failed.");
    exit(-1);
  }

  /*  read codon sequences  */
  for (i=0;i<n_cod;i++){

    fscanf(in,"%s ",cod_seq);
    cod_seq[3]='\0';
    /*  assign codon pointer to appropriate codon   */
    k=0;

    while(strcmp(cod_seq,cod[k].seq)!=0){
      k++;
      if(k>64) {
	k=-1;
	fprintf(stderr,"The codon listed as %s in the codon usage data file has no match \nin the genetic code file.  Check both files and adjust accordingly.\n",cod_seq);

	exit(-1);
      };
    };
    for(j=0;j<n_seq;j++) {
      seq[j].cod_data[i].cod=&cod[k];
    }

  };

    printf("\n");
    /*  establish pointers from aa_stats to amino acids  */
  for(i=0;i<n_aa;i++){
    for(j=0;j<n_seq;j++) {
      seq[j].a_stat[i].aa=&aa[i];
    }

  }

  /*  Connect seq.a_stat to its respective codons  */
  for(j=0;j<n_seq;j++) {
    for(k=0;k<n_aa;k++) {
      l=0;
      for(i=0;i<n_cod;i++) {
	if(seq[j].a_stat[k].aa==seq[j].cod_data[i].cod->aa){
	  seq[j].a_stat[k].cods[l]=&seq[j].cod_data[i];
	  l++;
	}
      }
    }
  }

  /*  read sequence name and codon occurence data  */
  for (j=0;j<n_seq;j++){


    /*  grab a character  */
    ch=fgetc(in);

    /*  if it's line feed, tab, carriage return or a space, keep grabbing  */
    while((ch==13)||(ch==9)||(ch==10)||(ch==' ')){
      ch=fgetc(in);
    }

    /*  then assign to the name until you hit a colon (or LF or CR)  */
    for(k=0;k<MAX_NAME_LENGTH;k++) {

      if((ch==13)||(ch==10)||(ch=='>')) {
      	seq[j].name[k]='\0';
	break;
      } else {
	seq[j].name[k]=ch;
	ch=fgetc(in);
      }
    }


    /*  read the data and count total codons  */
    seq[j].n_cod_obs=0;
    for (i=0;i<n_cod;i++){
      fscanf(in, "%d",&seq[j].cod_data[i].n);
      seq[j].n_cod_obs+=seq[j].cod_data[i].n;
    }

  }

  printf("Closing codon usage data file...\n");

  fclose(in);

  if(verbosity>0) {
    printf("----Data from codon usage data file----\n");
    printf("\nNum sequences: %d   Num codons: %d\n", n_seq,n_cod);

    for(i=0;i<n_seq;i++){
      printf("\n%s (Length: %d) ",seq[i].name,seq[i].n_cod_obs*3);
      for(k=0;k<n_cod;k++){
	printf("%d ", seq[i].cod_data[k].n);
      };
    }
    printf("\n");
    for(i=0; i<n_cod;i++){
      printf("%s ",seq[0].cod_data[i].cod->seq);
    }

    printf("\n--------------------------------------\n\n");
    if(verbosity>1) {
      printf("Press enter to continue...\n");
      do
	ch=getchar();
      while(ch!='\n');
    }

  }

}


void get_acgt(char *filename) {

  FILE *out;
  FILE *in;
  int i,j,k,l;

  char ch;

  /*  open a filestream    */
  if(!(in=fopen(filename,"r"))) {
    fprintf(stderr,"Error opening nucleotide composition file %s\n", filename);
    ch=getchar();
    exit(-1);
  }

  printf( "Reading nucleotide composition data file...\n");
  /*  grab a character  */
  ch=fgetc(in);
  if(ch=='C'){
    expected_given=1;
    printf("Reading file as expected frequencies for each codon...\n");
  } else if (ch=='N'){
    expected_given=0;
    printf("Reading file as background nucleotide compositions...\n");
  }
  /*  then read until you hit a LF or CR)  */
  for(k=0;k<MAX_CODE_LENGTH;k++) {
    if((ch==13)||(ch==10)) {
      break;
    } else {
      ch=fgetc(in);
    }
  }

  for(i=0;i<n_seq;i++){

    /* Read until you hit the colon */
    for(k=0;k<MAX_NAME_LENGTH;k++) {
      if(ch=='>') {
	break;
      } else {
	ch=fgetc(in);
      }
    }


    if(expected_given){
      for(k=0;k<n_cod;k++) {
	fscanf(in," %lf",&seq[i].cod_data[k].ex);
      }
    }else{
      fscanf(in," %lf %lf %lf %lf",&seq[i].f_A,&seq[i].f_C,&seq[i].f_G,&seq[i].f_T);
    }

    /*  then read until you hit a LF or CR)  */
    for(k=0;k<MAX_CODE_LENGTH;k++) {
      if((ch==13)||(ch==10)) {
	break;
      } else {
	ch=fgetc(in);
      }
    }

  }
  printf("Closing nucleotide composition data file...\n");

  fclose(in);
  if (verbosity>0) {
    printf("\n----Data from nucleotide composition data file----\n");

    for(i=0;i<n_seq;i++){
      if(expected_given){
	printf("%s ",seq[i].name);
	for(k=0;k<n_cod;k++){
	  printf("%.3f ",seq[i].cod_data[k].ex);
	}
	printf("\n");
      }else{
	printf("%s %f %f %f %f\n",seq[i].name, seq[i].f_A, seq[i].f_C, seq[i].f_G,seq[i].f_T);
      }
    }
    printf("------------------------------\n");
  }
  if(verbosity>1) {
    printf("Press enter to continue...\n");
    do
      ch=getchar();
    while(ch!='\n');
  }

}

void calc_cod_data(){

  int h,i,j,k,l,n;
  int p_denom;
  double sort[MAX_REDUN];
  double b,sum;
  int incr;
  int o;

  /*  first calculate the numerator of the expected frequency  */

  /*  cycle through the taxa  */
  if(expected_given){
      /* Check input frequencies... */

    printf("Checking expected frequencies sum to 1 for each amino acid.\n");

    for(j=0;j<n_seq;j++) {
      for(i=0;i<n_aa;i++){
	sum=0;
	for(k=0;k<aa[i].redun;k++)
	  sum+=seq[j].a_stat[i].cods[k]->ex;
	/* Catch egregioius violations */
	if((sum>1.05)||(sum<0.95)){
	  fprintf(stderr,"Expected frequencies for the codons of amino acid %c of sequence %s sum to %f not 1.\n",aa[i].abbrev,seq[j].name,sum);
	  for(k=0;k<aa[i].redun;k++){
	    fprintf(stderr,"%s %f\n",aa[i].cods[k]->seq,seq[j].a_stat[i].cods[k]->ex);
	  }
	  exit(-1);
	}else{
	  /* and for the rest, just renormalize */
	  for(k=0;k<aa[i].redun;k++){
	    seq[j].a_stat[i].cods[k]->ex/=sum;
	  }
	}
      }

    }
  }else {
    /*  calculate expected frequencies...  */
    printf("Calculating expected frequencies...\n");


    for(j=0;j<n_seq;j++){


      /*  cycle through the codons  */
      for(i=0;i<n_cod;i++) {
	/*  initialize ex_numerator  */
	seq[j].cod_data[i].ex_num=1;
	/*  and cycle through the sites  */
	for(k=0;k<3;k++){
	  /*  if variable site  */
	  if(seq[j].cod_data[i].cod->var_sites[k]){
	    /*  if A, numerator = numerator * f_A  */
	    if(toupper(seq[j].cod_data[i].cod->seq[k])=='A') seq[j].cod_data[i].ex_num*=seq[j].f_A;
	    if(toupper(seq[j].cod_data[i].cod->seq[k])=='T') seq[j].cod_data[i].ex_num*=seq[j].f_T;
	    if(toupper(seq[j].cod_data[i].cod->seq[k])=='C') seq[j].cod_data[i].ex_num*=seq[j].f_C;
	    if(toupper(seq[j].cod_data[i].cod->seq[k])=='G') seq[j].cod_data[i].ex_num*=seq[j].f_G;
	  }
	}
      }

      /*  cycle through the codons again, this time summing up the ex_nums...  */
      for(i=0;i<n_cod;i++) {
	/*  initialize the denominator  */
	seq[j].cod_data[i].ex_denom=0;
	for(k=0;k<n_cod;k++){
	  /*  if two codons belong to the same amino acid  */
	  if(toupper(seq[j].cod_data[i].cod->aa->abbrev)==toupper(seq[j].cod_data[k].cod->aa->abbrev)) {
	    /*  sum their numerators...  */
	    seq[j].cod_data[i].ex_denom=seq[j].cod_data[i].ex_denom+seq[j].cod_data[k].ex_num;
	  }
	}
      }

      for(i=0;i<n_cod;i++){
	if(seq[j].cod_data[i].ex_denom==0.0) seq[j].cod_data[i].ex_denom=1;
	seq[j].cod_data[i].ex=seq[j].cod_data[i].ex_num/seq[j].cod_data[i].ex_denom;
      }
    }/*  end cycling through sequences   */

  }

  printf("Calculating observed frequencies...\n");

  /*  now calculate p  */
  for(j=0;j<n_seq;j++){
    for(i=0;i<n_cod;i++){
      p_denom=0;
      for(k=0;k<n_cod;k++){
	/*  if two codons belong to the same amino acid  */
	if(toupper(seq[j].cod_data[i].cod->aa->abbrev)==toupper(seq[j].cod_data[k].cod->aa->abbrev)) {
	  /*  sum their numerators...  */
	  p_denom=p_denom+seq[j].cod_data[k].n;
	}
      } /*  compare to the next codon  */

      if(p_denom==0) {
	seq[j].cod_data[i].p=0;
      }	else {
	/*  calc the p value  */
	  seq[j].cod_data[i].p=(double)seq[j].cod_data[i].n / p_denom;
	}

      /*  calc the bias value (observed-expected)   */
      seq[j].cod_data[i].bias=seq[j].cod_data[i].p-seq[j].cod_data[i].ex;

    } /*  cycle to the next codon  */
  }
}


void calc_n_aa(){

  int i,j,k;
  printf("Counting number of each amino acid observed...\n");

  for(k=0;k<n_seq;k++) {
    for (i=0;i<n_aa;i++) {
      seq[k].a_stat[i].n=0;
      for(j=0;j<aa[i].redun;j++) {
	seq[k].a_stat[i].n+=seq[k].a_stat[i].cods[j]->n;
      }
    }
  }
}


void test_null() {

  int i,j,k,l;

  printf("Calculating Chi-square...\n");

  for(j=0;j<n_seq;j++){

    for(i=0;i<n_aa;i++){

      seq[j].a_stat[i].chisq=0.0;

      seq[j].a_stat[i].Bg=0.0;

      if(aa[i].redun<=1) {
	seq[j].a_stat[i].chisq=-2.0;
	continue;
      }

      /*  add up n*bias^2 / ex to get chisq values  */

      for(k=0;k<n_cod;k++){

	/* if(seq[j].cod_data[k].ex==0) seq[j].cod_data[k].ex=.01;  */
	if(aa[i].abbrev==seq[j].cod_data[k].cod->aa->abbrev) {

	  /* option to catch extremely low expected values	   */
  	  if(seq[j].cod_data[k].ex<eps*eps) {
	    seq[j].a_stat[i].chisq=-1;
	    seq[j].a_stat[i].Bg=0;
	    printf("A very small expectation (less than %lf).  Sequence: %s Codon: %s, Expectation: %lf \n This is not necessarily a problem, but it usually only occurs if an error reading the input files has occurred.  Check to make sure nucleotide composition file is being read in correctly by setting verbosity to 1.  Also, check to make sure the order of sequences in the same in both the nucleotide content file and the codon counts file.\n",eps*eps,seq[j].name,seq[j].cod_data[k].cod->seq,seq[j].cod_data[k].ex);
	    seq[j].a_stat[i].bad_ex=1;
	  }else{
	    seq[j].a_stat[i].chisq+=(pow(seq[j].cod_data[k].bias,2.0)/seq[j].cod_data[k].ex)*seq[j].a_stat[i].n;
	    if(seq[j].cod_data[k].p!=seq[j].cod_data[k].ex) {
	      seq[j].a_stat[i].Bg+=fabs((seq[j].cod_data[k].p/seq[j].cod_data[k].ex)-1.0);

	    }else {
	      seq[j].a_stat[i].Bg=0;
	    }
	  };
	}

      } /*  end codon cycle  */

      /*  calculate p value from gamma distribution (X^2)  */
      if(!seq[j].a_stat[i].bad_ex)
      	seq[j].a_stat[i].p=1-GammaP(0.5*(aa[i].redun-1), 0.5*(seq[j].a_stat[i].chisq));
      else
      	seq[j].a_stat[i].p=1;

      /* if p < 0.05, then reject null  */
      if(seq[j].a_stat[i].p<alpha) {
		seq[j].a_stat[i].m=1;
      } else {
      	seq[j].a_stat[i].m=0;
      }

      /*  Divide Bg value by redudancy of amino acid  */
      seq[j].a_stat[i].Bg/=(double)aa[i].redun;

    } /*  end aa cycle  */
  } /*  end sequence cycle    */
}


void calc_F(){

  int i,j,k,l,n;
  double sum_p2;
  /*  see Frank Wright paper 1990 Gene 87:23-29   */

  for(j=0;j<n_seq;j++) {
    /*  cycle through all amino acids and calc F for each  */
    for(k=0;k<n_aa;k++) {

      n=seq[j].a_stat[k].n;
      seq[j].a_stat[k].F=0;
      seq[j].a_stat[k].F_p=0;

      sum_p2=0;

      for(i=0;i<aa[k].redun;i++){
	sum_p2+=pow((seq[j].a_stat[k].cods[i]->p),2.0);
      }

      /*  if enough observations that chi-sq is legit (i.e. >5)  */
      /*  calculate F'  */
      /*  else set F'= 1/redun... as a conservative assumption  */
      if(seq[j].a_stat[k].n>5){
	seq[j].a_stat[k].F_p=(seq[j].a_stat[k].chisq+(double)n-(double)aa[k].redun)/((double)aa[k].redun*((double)n-1.0));
	seq[j].a_stat[k].F=((double)n*sum_p2-1.0)/((double)n-1.0);
      }
      else {
	seq[j].a_stat[k].F=-1;
	seq[j].a_stat[k].F_p=-1;
      }
    }
  }
}


void calc_sumchi() {
  int j,i;
  /*  Calculate the Sum of Chi-Squares measure.  In addition, calculate B*(a) of Karlin and Mrazek.    */
  /*  Here B*(a) is denoted as B_KM.  */

  for(j=0;j<n_seq;j++) {
    seq[j].sumchi=0;
    seq[j].B_KM=0;

    for(i=0;i<n_aa;i++){
      if(seq[j].a_stat[i].chisq>-1) {
	seq[j].sumchi+=seq[j].a_stat[i].chisq;
	seq[j].df+=seq[j].a_stat[i].aa->redun-1;
      }
      seq[j].B_KM+=((double)seq[j].a_stat[i].n/(double)seq[j].n_cod_obs)*seq[j].a_stat[i].Bg;
    }

    seq[j].p=1-GammaP(0.5*(seq[j].df), 0.5*(seq[j].sumchi));
    seq[j].scaledchi=seq[j].sumchi/seq[j].n_cod_obs;
  }
}


void calc_Nc() {

  double F_tot[MAX_REDUN+1];
  double F_tot_p[MAX_REDUN+1];
  double F_bar[MAX_REDUN+1];
  double F_bar_p[MAX_REDUN+1];
  int n_fold[MAX_REDUN+1];
  int n_fold_p[MAX_REDUN+1];
  int n_fold_o[MAX_REDUN+1];
  int i,j,k,l;

  printf("Calculating Nc and Nc'...\n");
  for(k=0;k<n_seq;k++){

    seq[k].N_c=0;
    seq[k].N_c_p=0;

    /*  for each redundancy class, get average F and add to N_c  */
    for(i=1;i<(MAX_REDUN+1);i++) {
      n_fold[i]=0;
      n_fold_p[i]=0;
      n_fold_o[i]=0;
      F_tot[i]=0;
      F_tot_p[i]=0;

      /*  scan all amino acids to see if they are included in redundancy class  */
      for(j=0;j<n_aa;j++) {

	/*  if in class, not termination codon, and have been observed  */
	/*  add their contribution to the mean of the class  */
	if((aa[j].redun==i)&&(aa[j].abbrev!='*')/*&&(seq[k].a_stat[j].n>1)*/) {

	  /*  increment the number of aa acids in the redun class  */
	  n_fold_o[i]++;

	  /*  only if F was calculated above  */
	  if(seq[k].a_stat[j].F!=-1) {
	    F_tot[i]+=seq[k].a_stat[j].F;
	    n_fold[i]++;
	  }

	  /*  only if F_p was calculated above  */
	  if(seq[k].a_stat[j].F_p!=-1) {
	    F_tot_p[i]+=seq[k].a_stat[j].F_p;
	    n_fold_p[i]++;
	  }

	}
      }

    } /*  end redudancy class  */


    for(i=1;i<MAX_REDUN+1;i++) {

      /*  get F_bar  */
      if(n_fold[i])
        F_bar[i]=F_tot[i]/n_fold[i];
      else if(i==3 && n_fold[2]&&n_fold[4])
        /* Bug found by Anders Fuglsang  */
        /*F_bar[i]=(n_fold[2]/F_tot[2]+n_fold[4]/F_tot[4])*0.5;*/
        /* Bug fix by JN 2/28/06 */
        F_bar[i]=(F_tot[2]/n_fold[2]+F_tot[4]/n_fold[4])*0.5;
      else /*  assume equal usage to be conservative     */
        F_bar[i]=1/i;

      if(n_fold_p[i])
        F_bar_p[i]=F_tot_p[i]/n_fold_p[i];
      else if(i==3 && n_fold_p[2]&&n_fold_p[4])
        /* Bug found by Anders Fuglsang */
      	/*F_bar_p[i]=(n_fold_p[2]/F_tot_p[2]+n_fold_p[4]/F_tot_p[4])*0.5;*/
        /* Bug fix by JN 2/28/06 */
        F_bar_p[i]=(F_tot_p[2]/n_fold_p[2]+F_tot_p[4]/n_fold_p[4])*0.5;
      else /*  assume equal usage to be conservative  */
      	F_bar_p[i]=1/i;

      /*  Force F_bar to be bounded (1/k,1)  */

      /*  To prevent divide by zero's if F_bar[i]<eps set to 1/redun  */
      if(F_bar[i]<1/(double)i) F_bar[i]=1/(double)i;
      if(F_bar_p[i]<1/(double)i) F_bar_p[i]=1/(double)i;
      /*       if(k==5) printf("bar: %d %f\n",i,F_bar[i]);  */
      if(n_fold_o[i]){
	seq[k].N_c+=n_fold_o[i]/F_bar[i];
	seq[k].N_c_p+=n_fold_o[i]/F_bar_p[i];
      }
    }
  }
}


void output2screen(){

  int i,j,k;
  char c;

  double scaledchi_m=0,N_c_m=0,N_c_p_m=0;
  double scaledchi_v=0,N_c_v=0,N_c_p_v=0;
  double scaledchi_CV,N_c_CV,N_c_p_CV;
  double scaledchi_d,N_c_d,N_c_p_d;

  for(j=0;j<n_seq;j++){

    scaledchi_d=seq[j].scaledchi-scaledchi_m;
    N_c_d=seq[j].N_c-N_c_m;
    N_c_p_d=seq[j].N_c_p-N_c_p_m;

    scaledchi_m+=scaledchi_d/(j+1);
    N_c_m+=N_c_d/(j+1);
    N_c_p_m+=N_c_p_d/(j+1);

    scaledchi_v+=scaledchi_d*(seq[j].scaledchi-scaledchi_m);
    N_c_v+=N_c_d*(seq[j].N_c-N_c_m);
    N_c_p_v+=N_c_p_d*(seq[j].N_c_p-N_c_p_m);

  }
  scaledchi_v/=n_seq;
  N_c_v/=n_seq;
  N_c_p_v/=n_seq;

  scaledchi_CV=sqrt(scaledchi_v)/scaledchi_m;
  N_c_CV=sqrt(N_c_v)/N_c_m;
  N_c_p_CV=sqrt(N_c_p_v)/N_c_p_m;

  /* printf("\n%s@%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",outputfile,N_c_m,N_c_p_m,scaledchi_m,N_c_v,N_c_p_v,scaledchi_v,N_c_CV,N_c_p_CV,scaledchi_CV);  */

  if(!quick_mode) {
    printf("\n\n----Data Explorer----\n");

    while(1) {
      do{
	printf("\nWhich sequence?  Enter a # between 1 and %d, 0 to quit: \n",n_seq);
	scanf("%d",&j);

      }while((j<0)||(j>n_seq) ) ;
      if(j==0) return;
      j--;
      printf("=================================================\n");
      printf("Name: %s ",seq[j].name);
      if(!expected_given)
	printf(" Nuc content:%lf %lf %lf %lf",seq[j].f_A, seq[j].f_C, seq[j].f_G,seq[j].f_T);

      printf("\n\tAkashi scaled Sum of X^2 values: %lf", seq[j].scaledchi);

      printf("\n\tKarlin and Mrazek B*(a): %lf",seq[j].B_KM);
      printf("\n\tNc:  %lf \n\tNcp: %lf\n", seq[j].N_c, seq[j].N_c_p);
      printf("\tNum cods observed: %d\n",seq[j].n_cod_obs);

      /*  clean the newline off the buffer  */
      c=getchar();

      while(1){
	do{
	  printf("Which amino acid? Enter the one-letter abbrev (0 for new sequence): \n");
	  c=getchar();
	  while(getchar()!='\n') continue;
	  for(i=0;i<(n_aa);i++){
	    if(aa[i].abbrev==toupper(c)) {
	      k=i;
	      break;
	    }
	    k=-1;
	  }
	}while((k==-1)&&(c!='0'));

	if(c=='0') break;

	printf("Amino acid: %c \t Redun: %d  Obs: %d  Chisq: %lf \n", aa[k].abbrev, aa[k].redun, seq[j].a_stat[k].n,seq[j].a_stat[k].chisq);
	printf("   F: %lf F_p: %lf\n",seq[j].a_stat[k].F,seq[j].a_stat[k].F_p);
	for(i=0;i<aa[k].redun;i++){
	  printf("   %s %c ", aa[k].cods[i]->seq, aa[k].abbrev);
	  printf("%4d ", seq[j].a_stat[k].cods[i]->n);
	  printf("f_obs: %.4lf f_ex: %.4lf diff: %.4lf\n", seq[j].a_stat[k].cods[i]->p, seq[j].a_stat[k].cods[i]->ex,seq[j].a_stat[k].cods[i]->bias);

	}

      }
    }
  }
  printf("Finished...\n");
}


void output2file(char *filename) {

  FILE *out;
  int i,j,k,l;
  char ch;
  /* strcat(filename,".stats");*/
  if(!(out=fopen(filename,"w"))){
    fprintf(stderr,"Error opening output file %s\n",filename);
  	ch=getchar();
  }

  fprintf(out,"Name Nc Ncp ScaledChi SumChi df p B_KM n_codons\n");
  for(j=0;j<n_seq;j++) {
    fprintf(out, "%s: %.4f %.4f %.4f %.4f %d %.4f %.4f %d\n", seq[j].name, seq[j].N_c,seq[j].N_c_p,seq[j].scaledchi,seq[j].sumchi,seq[j].df,seq[j].p,seq[j].B_KM,seq[j].n_cod_obs);
  }

  fclose(out);

}


void init_codes(){

  /* These represent genetic codes in Genbank format, and with corresponding ID numbers */
  sprintf(codes[1],"%s",
"   AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = ---M---------------M---------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
");
  sprintf(codes[2],"%s",
 "  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG \
  Starts = --------------------------------MMMM---------------M------------ \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
");
  sprintf(codes[3],"%s",
 "  AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = ----------------------------------MM---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[4],"%s",
  " AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = --MM---------------M------------MMMM---------------M------------ \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[5],"%s",
  " AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG \
  Starts = ---M----------------------------MMMM---------------M------------ \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[6],"%s",
  " AAs  = FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  /* Note: codes 7 and 8 do not exist */
  sprintf(codes[9],"%s",
  " AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------M------------ \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[10],"%s",
  " AAs  = FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[11],"%s",
  " AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = ---M---------------M------------MMMM---------------M------------ \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[12],"%s",
  " AAs  = FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = -------------------M---------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[13],"%s",
  " AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[14],"%s",
  " AAs  = FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[15],"%s",
  " AAs  = FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[16],"%s",
 "  AAs  = FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
");
  /* Note: codes 17-20 do not exist */
  sprintf(codes[21],"%s",
  " AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------M------------ \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[22],"%s",
  " AAs  = FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = -----------------------------------M---------------------------- \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");
  sprintf(codes[23],"%s",
  " AAs  = FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG \
  Starts = --------------------------------M--M---------------M------------ \
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG \
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG \
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG \
 ");


}






/* The following code was provided kindly by Eric C. Anderson */


/*************
*  FUNCTION:  double GammaP
*  PURPOSE:   Compute the incomplete Gamma function.  P(a,x)


*  INPUT:     double a;
			  double x;

   	Requires the functions IncGammaCF and IncGammaSer

   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
double
GammaP(double a, double x) {
	if(x<0.0 || a <= 0.0)  {
		printf("\n\nx<0.0 or a<=0.0 in GammaP\n\nExiting to system...");
		exit(1);
	}
	if(x <= a + 1.0) {
		return(IncGammaSer(a,x));
	}
	else {
		return(1.0 - IncGammaCF(a,x));
	}
}



/*************
*  FUNCTION:  double IncGammaSer
*  PURPOSE:   Compute the incomplete Gamma function by its series representation.
			  Note that this computes the part from 0 to x (called P(a,x) )

			  This is the preferred method when x < (a+1)

*  INPUT:     double x

	This function is based on the series representation described
	in "Numerical Recipes in C" (Cambridge University Press),

   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
#define MAXIT 100
#define PRECIS 1.0e-7

double
IncGammaSer(double a, double x)
{
  double old_funct, new_funct;  /*  this is what we ultimately want to return;  */
  double term;  /*  the next term to be added to the series  */
	double i;
	double numer;  /*  the normalizing constant---Gamma(a).  Also turns out to be the numerator of the terms that  */
	/*  get summed.  */
	double denom;  /*  the denominator of the terms that get summed  */

	/*  set value of the log of the normalizing constant  */
	numer = exp(LogGamma(a));
	denom = a * numer;

	term = numer/denom;
	old_funct = term;

	for(i=1;i<MAXIT;i++)  {
	  term *= 1.0/(a + i) * x;  /*  the terms may be computed recursively  */
		new_funct = old_funct + term;

		/*  check for convergence  */
		if( fabs(new_funct - old_funct) < PRECIS)  {
			break;
		}
		else {
			old_funct = new_funct;
		}
	}

	if(i >= MAXIT - .001) {
	  /*  done here we warn if it failed to converge.   */
		printf("\n\nFailed to converge in IncGammaSer\n\n");
		exit(1);
	}

	/*  otherwise, we return the result  */
	/*  first we have to add the coefficients and normalize it:  */
	new_funct *= exp(-x) * pow(x,a) * 1.0/numer;
	return(new_funct);
}
#undef MAXIT
#undef PRECIS

/*************
*  FUNCTION:  double IncGammaCF
*  PURPOSE:   Compute the incomplete Gamma function by its continued
			  fraction representation.  Note that this computes the
			  part from x to infinity (called Q(a,x) )

			  This is the preferred method for x > (a + 1)

*  INPUT:     double x

	This function is based on the continued fraction representation described
	in "Numerical Recipes in C" (Cambridge University Press), and evaluated using
	Lentz's method for continued fractions.  Basically following the pseudocode on
	page 171.

   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
#define SMIDGEN 1.0e-30
#define MAXIT 100
#define PRECIS 2.0e-7

double
IncGammaCF(double a, double x)
{
  double funct;  /*  this is what we ultimately want to return;  */
	double A,B,C,D;
	double Delta;
	double j;
	double log_normo;  /*  the log of the normalizing constant---Gamma(a)  */

	/*  set value of the log of the normalizing constant  */
	log_normo = LogGamma(a);

	/*  here for the zero subscript part:  */
	B = 0.0;  /*  b_0 is really zero, so we just make it tiny  */
	funct = SMIDGEN;
	C = funct;
	D = 0.0;

	/*  then for the 1 subscript part things are sort of different still  */
	A = 1.0;
	B = x + 1.0 - a;
	D = B + A * D;
	if(fabs(D) < SMIDGEN) D = SMIDGEN;
	C = B +  A/C;
	if(fabs(C) < SMIDGEN) C = SMIDGEN;
	D = 1.0/D;
	Delta = C * D;
	funct = Delta * funct;
	/*  don't even bother checking for convergence after this first "iteration"  */


	/*  now we iterate through the 2,3... and so forth subscript parts until converged  */
	/*  so each loop corresponds to the j+1-th subscript of the continued fraction  */
	for(j=1.0;j<MAXIT;j++)  {
	  /*  now, we define the new values for A and B on the j-th level of the continued fraction  */
		A = -j * (j - a);
		B = x + 1.0 + (2 * j) - a;

		D = B + A * D;
		if(fabs(D) < SMIDGEN) D = SMIDGEN;
		C = B + A/C;
		if(fabs(C) < SMIDGEN) C = SMIDGEN;
		D = 1.0/D;
		Delta = C * D;
		funct = Delta * funct;

		/*  here we check for convergence.  If it has, we return the appropriate result  */
		if(fabs(Delta - 1.0) < PRECIS) {
		  /*  add the coefficients and the normalizing constant  */
			funct *= exp(-x - log_normo) * pow(x,a);
			return(funct);
		}
	}

	/*  done here we warn if it failed to converge.  Should never get here  */
	printf("\n\nFailed to converge in IncGammaCF\n\n");
	exit(1);

	return(-55.55);  /*  just put this in so that it doesn't give a compile warning  */
	/*  about having no return value.  */

}
#undef SMIDGEN
#undef MAXIT
#undef PRECIS

/*************
*  FUNCTION:  double LogGamma
*  PURPOSE:   Compute log of the Gamma Function

*  INPUT:     double x

	This function is based on the six term series of Lanczos as described
	in "Numerical Recipes in C" (Cambridge University Press).  I use the
	choice of gamma = 5 and N = 6.

*************/
double
LogGamma(double x)
{
  /*  declare the coefficients in the series as a static double array  */
	static double Coeff[6] = {76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	double Prelim;  /*  for the part before the series  */
	double Series;
	int i;
	double denom;

	/*  We start off by computing Gamma(x+1) by the series:  */

	/*  compute the preliminary part  */
	Prelim = (x+.5) * log(x+5.5) - (x+5.5);

	/*  add the log of sqrt(2\pi) to that:  */
	Prelim += 0.91893853320467274;

	/*  now compute the series.  Start with the constant term c_0  */
	Series = 1.000000000190015;

	denom = x;
	for(i=0;i<6;i++)  {
		Series += Coeff[i]/(++denom);
	}

	/*  now we just have to return the right thing.  Since we just computed  */
	/*  Recall Gamma(x+1) = x * Gamma(x), so we have to subtract log(x) from this.  */
	/*  We can do this with only one call to log by dividing Series by x.  */
	return(Prelim + log(Series/x));


}
