
#include <stdio.h>
#include <signal.h>

//#include "R.h"
#include "bootdist.h"
#include "stdlib.h"


#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
#endif /* OSX_CARBON */

#ifdef WIN32
#include <windows.h>
/* for console code (clear screen, text color settings) */
CONSOLE_SCREEN_BUFFER_INFO      savecsbi;
boolean savecsbi_valid = false;
HANDLE  hConsoleOutput;

#endif /* WIN32 */

#ifndef OLDC
static void crash_handler(int signum);

#endif

#if defined(OSX_CARBON) && defined(__MWERKS__)
boolean fixedpath = false;
#endif


FILE *outfile, *infile, *intree, *intree2, *outtree, *weightfile, *catfile, *ancfile, *mixfile, *factfile;
long spp, words, bits, nonodes, endsite, outgrno, nextree, which;
boolean ibmpc, ansi, tranvsp, interleaved, printdata, outgropt, treeprint, dotdiff, transvp;
naym *nayme;                     /* names of species */
steptr weight, category, alias, location, ally;
sequence y;


#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   seqboot_inputnumbers(void);
void   inputoptions(void);
char **matrix_char_new(long rows, long cols);
void   matrix_char_delete(char **mat, long rows);
double **matrix_double_new(long rows, long cols);
void   matrix_double_delete(double **mat, long rows);
void   seqboot_inputdata(void);
void   allocrest(void);
void   freerest(void);
void   allocnew(void);
void   freenew(void);
void   allocnewer(long newergroups, long newersites);
void   doinput(int argc, Char *argv[]);
void   writebootmatrix(void);
void   bootweights(void);
void   writebootdist(void);
void   bootwrite(void);
void   freenewer(void);
/* function prototypes */
#endif

/*** Config vars ***/
/* Mutually exclusive booleans for boostrap type */
boolean bootstrap;

/* Bootstrap sample frequency */
double fracsample = 0.5; /* ...or user-defined sample freq, [0..inf) */

boolean weights = false;/* Read weights file */

boolean firstrep; /* TODO Must this be global? */
longer seed;

/* Filehandles and paths */
/* Usual suspects declared in phylip.c/h */
Char infilename[FNMLNGTH], outfilename[FNMLNGTH];
Char *tmp;
long sites, loci, maxalleles, groups, 
  reps, maxnewsites;

steptr oldweight, where, how_many, mixdata, ancdata;

/* Original dataset */
/* [0..spp-1][0..sites-1] */
Char **nodep   = NULL;           /* molecular or morph data */
double **nodef = NULL;         /* gene freqs */

Char *factor = NULL;  /* factor[sites] - direct read-in of factors file */
long *factorr = NULL; /* [0..sites-1] => nondecreasing [1..groups] */

long *alleles = NULL;

/* Mapping with read-in weights eliminated
 * Allocated once in allocnew() */
long newsites;
long newgroups;
long *newwhere   = NULL;    /* Map [0..newgroups-1] => [1..newsites] */
long *newhowmany = NULL;    /* Number of chars for each [0..newgroups-1] */

/* Mapping with bootstrapped weights applied */
/* (re)allocated by allocnewer() */
long newersites, newergroups;
long *newerfactor  = NULL;  /* Map [0..newersites-1] => [1..newergroups] */
long *newerwhere   = NULL;  /* Map [0..newergroups-1] => [1..newersites] */
long *newerhowmany = NULL;  /* Number of chars for each [0..newergroups-1] */
long **charorder   = NULL;  /* Permutation [0..spp-1][0..newergroups-1] */
long **sppord      = NULL;  /* Permutation [0..newergroups-1][0..spp-1] */

void getoptions()
{
  /* interactively set options */
  long inseed, inseed0, loopcount;
  bootstrap = true;
  fracsample = 1.0;

  weights = false;
  interleaved = true;
  loopcount = 0;
  for (;;) {
    //cleerhome();
    printf("\nAlgorithms revised from PHYLIP, version %s\n\n",VERSION);
		
#ifdef WIN32
#endif
    fflush(stdout);
	//scanf("%ld%*[^\n]", &reps);
    //getchar();
	//printf("%d", reps);
	reps = atoi(tmp);
    if (reps == (int)(atoi(tmp))) {
		if(reps <= 0){
			printf("\n\nBAD NUMBER: We recommand a positive number > 100.\n\t(Set to default 1000)\n");
			reps = 1000;
		}else if(reps < 100){
			printf("\n\nWARNING: We recommand a positive number > 100.\n");
		}
	}else{
			printf("\n\nBAD NUMBER: Input uncorrect. \n\t(Set to default 1000)\n");
			reps = 1000;
	}
    break;

    countup(&loopcount, 100);
  }
  fracsample = 1.0;
  
  initseed(&inseed, &inseed0, seed);

}  /* getoptions */


void seqboot_inputnumbers()
{
  /* read numbers of species and of sites */
  fscanf(infile, "%ld%ld", &spp, &sites);
  loci = sites;
  maxalleles = 1;

}  /* seqboot_inputnumbers */


void inputoptions()
{
  /* input the information on the options */
  long maxfactsize, i;

  for (i = 1; i <= (sites); i++)
    factorr[i - 1] = i;

  for (i = 0; i < (sites); i++)
    oldweight[i] = 1;

  for (i = 0; i < (loci); i++)
    how_many[i] = 0;
  for (i = 0; i < (loci); i++)
    where[i] = 0;
  for (i = 1; i <= (sites); i++) {
    how_many[factorr[i - 1] - 1]++;
    if (where[factorr[i - 1] - 1] == 0)
      where[factorr[i - 1] - 1] = i;
  }
  groups = factorr[sites - 1];
  newgroups = 0;
  newsites = 0;
  maxfactsize = 0;
  for(i = 0 ; i < loci ; i++){
    if(how_many[i] > maxfactsize){
      maxfactsize = how_many[i];
    }
  }
  maxnewsites = groups * maxfactsize;
  allocnew();
  for (i = 0; i < groups; i++) {
    if (oldweight[where[i] - 1] > 0) {
      newgroups++;
      newsites += how_many[i];
      newwhere[newgroups - 1] = where[i];
      newhowmany[newgroups - 1] = how_many[i];
    }
  }
}  /* inputoptions */


char **matrix_char_new(long rows, long cols)
{
  char **mat;
  long i;

  assert(rows > 0); assert(cols > 0);

  mat = (char **)Malloc(rows*sizeof(char *));
  for (i = 0; i < rows; i++)
    mat[i] = (char *)Malloc(cols*sizeof(char));

  return mat;
}


void matrix_char_delete(char **mat, long rows)
{
  long i;
  
  assert(mat != NULL);
  for (i = 0; i < rows; i++)
    free(mat[i]);
  free(mat);
}


double **matrix_double_new(long rows, long cols)
{
  double **mat;
  long i;

  assert(rows > 0); assert(cols > 0);
  
  mat = (double **)Malloc(rows*sizeof(double *));
  for (i = 0; i < rows; i++)
    mat[i] = (double *)Malloc(cols*sizeof(double));

  return mat;
}


void matrix_double_delete(double **mat, long rows)
{
  long i;
  
  assert(mat != NULL);
  for (i = 0; i < rows; i++)
    free(mat[i]);
  free(mat);
}


void seqboot_inputdata()
{
  /* input the names and sequences for each species */
  long i, j, k, l, m, n, basesread, basesnew=0;
  Char charstate;
  boolean allread, done;

  nodep = matrix_char_new(spp, sites);
  
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  
  interleaved = false;

  basesread = 0;
  allread = false;
  while (!allread) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile)) 
      scan_eoln(infile);
    i = 1;
    while (i <= spp) {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i-1);
      j = interleaved ? basesread : 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) ||eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ')
            continue;
          uppercase(&charstate);
          j++;
          if (charstate == '.')
            charstate = nodep[0][j-1];
          nodep[i-1][j-1] = charstate;
        }
        if (interleaved)
          continue;
        if (j < sites) 
          scan_eoln(infile);
        else if (j == sites)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j;
      scan_eoln(infile);
      if ((interleaved && j != basesnew) || ((!interleaved) && j != sites)){
        printf("\n\nERROR: sequences out of alignment at site %ld", j+1);
        printf(" of species %ld\n\n", i);
        exxit(-1);}
      i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == sites);
    } else
      allread = (i > spp);
  }
  return;

  m = (sites - 1) / 60 + 1;
  
  for (i = 1; i <= m; i++) {
    for (j = 0; j < spp; j++) {
      for (k = 0; k < nmlngth; k++)
        putc(nayme[j][k], outfile);
      fprintf(outfile, "   ");
      l = i * 60;
      if (l > sites)
        l = sites;
      n = (i - 1) * 60;
      for (k = n; k < l; k++) {
        if (j + 1 > 1 && nodep[j][k] == nodep[0][k])
          charstate = '.';
        else
          charstate = nodep[j][k];
        putc(charstate, outfile);
        if ((k + 1) % 10 == 0 && (k + 1) % 60 != 0)
          putc(' ', outfile);        
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* seqboot_inputdata */


void allocrest()
{ /* allocate memory for bookkeeping arrays */

  oldweight = (steptr)Malloc(sites*sizeof(long));
  weight = (steptr)Malloc(sites*sizeof(long));
  where = (steptr)Malloc(loci*sizeof(long));
  how_many = (steptr)Malloc(loci*sizeof(long));
  factor = (Char *)Malloc(sites*sizeof(Char));
  factorr = (steptr)Malloc(sites*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
}  /* allocrest */

void freerest()
{
  /* Free bookkeeping arrays */
  free(oldweight);
  free(weight);
  free(where);
  free(how_many);
  free(factor);
  free(factorr);
  free(nayme);
}


void allocnew(void)
{ /* allocate memory for arrays that depend on the lenght of the 
     output sequence*/
  /* Only call this function once */
  assert(newwhere == NULL && newhowmany == NULL); 

  newwhere = (steptr)Malloc(loci*sizeof(long));
  newhowmany = (steptr)Malloc(loci*sizeof(long));
} /* allocnew */


void freenew(void)
{ /* free arrays allocated by allocnew() */
  /* Only call this function once */
  assert(newwhere != NULL);
  assert(newhowmany != NULL);

  free(newwhere);
  free(newhowmany);
}


void allocnewer(long newergroups, long newersites)
{ /* allocate memory for arrays that depend on the length of the bootstrapped
     output sequence */
  /* Assumes that spp remains constant */
  static long curnewergroups = 0;
  static long curnewersites  = 0;

  long i;

  if (newerwhere != NULL) {
    if (newergroups > curnewergroups) {
      free(newerwhere);
      free(newerhowmany);
      for (i = 0; i < spp; i++)
        free(charorder[i]);
      newerwhere = NULL;
    }
    if (newersites > curnewersites) {
      free(newerfactor);
      newerfactor = NULL;
    }
  }

  if (charorder == NULL)
    charorder = (steptr *)Malloc(spp*sizeof(steptr));

  /* Malloc() will fail if either is 0, so add a dummy element */
  if (newergroups == 0)
    newergroups++;
  if (newersites == 0)
    newersites++;
  
  if (newerwhere == NULL) {
    newerwhere = (steptr)Malloc(newergroups*sizeof(long));
    newerhowmany = (steptr)Malloc(newergroups*sizeof(long));
    for (i = 0; i < spp; i++)
      charorder[i] = (steptr)Malloc(newergroups*sizeof(long));
    curnewergroups = newergroups;
  }
  if (newerfactor == NULL) {
    newerfactor = (steptr)Malloc(newersites*sizeof(long));
    curnewersites = newersites;
  }
}


void freenewer()
{
  /* Free memory allocated by allocnewer() */
  /* spp must be the same as when allocnewer was called */
  long i;

  if (newerwhere) {
    free(newerwhere);
    free(newerhowmany);
    free(newerfactor);
    for (i = 0; i < spp; i++)
      free(charorder[i]);
    free(charorder);
  }
}


void doinput(int argc, Char *argv[])
{ /* reads the input data */
  tmp = argv[3]; 
  getoptions();

  seqboot_inputnumbers();
  allocrest();
  openfile(&outfile,argv[2],"output data file","w",argv[0],outfilename);
  //openfile(&outfile,OUTFILE,"output data file","w",argv[0],outfilename);
  inputoptions();
  seqboot_inputdata();
}  /* doinput */


void bootweights()
{ /* sets up weights by resampling data */
  long i, j, k, blocks;
  long grp = 0, site = 0;
  
  for (i = 0; i < (newgroups); i++)
    weight[i] = 0;
  if (bootstrap) {
    blocks = fracsample * newgroups;
    for (i = 1; i <= (blocks); i++) {
      j = (long)(newgroups * randum(seed)) + 1;
	  weight[j - 1]++;
      j++;
      if (j > newgroups)
        j = 1;
    }

  }
  /* Count number of replicated groups */
  newergroups = 0;
  newersites  = 0;
  for (i = 0; i < newgroups; i++) {
    newergroups += weight[i];
    newersites  += newhowmany[i] * weight[i];
  }

  if (newergroups < 1) {
    fprintf(stdout, "ERROR: sampling frequency or number of sites is too small\n");
    exxit(-1);
  }
  
  /* reallocate "newer" arrays, sized by output groups:
   * newerfactor, newerwhere, newerhowmany, and charorder */
  allocnewer(newergroups, newersites);
  
  /* Replicate each group i weight[i] times */
  grp = 0;
  site = 0;
  for (i = 0; i < newgroups; i++) {
    for (j = 0; j < weight[i]; j++) {
      for (k = 0; k < newhowmany[i]; k++) {
        newerfactor[site] = grp + 1;
        site++;
      }
      newerwhere[grp] = newwhere[i];
      newerhowmany[grp] = newhowmany[i];
      grp++;
    }
  }
}  /* bootweights */


void writebootmatrix()
{
  /* write out one set of bootstrapped sequences */
  long i, j, k, l, m, n;  //, n2;
  
  Char charstate;

  sppord = (long **)Malloc(newergroups*sizeof(long *));
  for (i = 0; i < (newergroups); i++)
    sppord[i] = (long *)Malloc(spp*sizeof(long));
  for (j = 1; j <= spp; j++)
    sppord[0][j - 1] = j;
  for (i = 1; i < newergroups; i++) {
    for (j = 1; j <= (spp); j++)
      sppord[i][j - 1] = sppord[i - 1][j - 1];
  }
  fprintf(outfile, "%5ld %5ld\n", spp, newersites);
  
  l = 1;
  
  long n2 = 0; 
  
  m = interleaved ? 60 : newergroups;
  do {
    if (m > newergroups)
      m = newergroups;
    for (j = 0; j < spp; j++) {
      n = 0;
      for (k = 0; k < nmlngth-n2; k++)
        fprintf(outfile, " "); 
      fprintf(outfile, " "); 
      
      for (k = l - 1; k < m; k++) {
        for (n2 = -1; n2 <= (newerhowmany[charorder[j][k]] - 2); n2++) {
          n++;
          charstate = nodep[sppord[charorder[j][k]][j] - 1]
                           [newerwhere[charorder[j][k]] + n2];
          putc(charstate, outfile);
          if (n % 10 == 0 && n % 60 != 0)
            putc(' ', outfile);
          
        }
      }
      putc('\n', outfile);
    }
  } while (interleaved && l <= newersites);
  for (i = 0; i < (newergroups); i++)
    free(sppord[i]);
  free(sppord);
}  /* writebootmatrix */


void writebootdist()
{
  /* write out one set of bootstrapped sequences */
  long i, j, k, l, m, n, n2;
  double doa, dob, dab;
  int na, nb, nab;
  
  sppord = (long **)Malloc(newergroups*sizeof(long *));
  for (i = 0; i < (newergroups); i++)
    sppord[i] = (long *)Malloc(spp*sizeof(long));
  for (j = 1; j <= spp; j++)
    sppord[0][j - 1] = j;
  for (i = 1; i < newergroups; i++) {
    for (j = 1; j <= (spp); j++)
      sppord[i][j - 1] = sppord[i - 1][j - 1];
  }
  
  //fprintf(outfile, "%d\n", spp);
  fprintf(outfile, "%ld\n", spp);
  
  l = 1;
  m = interleaved ? 60 : newergroups;
  
  int bootmatrix[spp][newersites];
  double distance[spp][spp];
  
  do {
    if (m > newergroups)
      m = newergroups;
    for (j = 0; j < spp; j++) {

      for (k = l - 1; k < m; k++) {		
        for (n2 = -1; n2 <= (newerhowmany[charorder[j][k]] - 2); n2++) {
		  //sppord: sample count
          n++;
          bootmatrix[j][k] = (int)nodep[sppord[charorder[j][k]][j] - 1]
                           [newerwhere[charorder[j][k]] + n2] - 48;
		  // k: site count
		  // bootmatrix: each site value after bootstrapping
        }
      }
    }
  } while (interleaved && l <= newersites);
  for (i = 0; i < (newergroups); i++)
    free(sppord[i]);
  free(sppord);
  
  for(i =0;i < spp; i++){
	  n = 0;
      if ((l == 1)) {
        n2 = nmlngth;
        //for (k = 0; k < n2; k++)
          //putc(nayme[i][k], outfile);  // j: species count; k: site count
		  //printf("%d\t%d\t\n", j,k); 
	  }
      for (k = 0; k < nmlngth-n2; k++)
        fprintf(outfile, " "); 
	  fprintf(outfile, "%s", nayme[i]);   //nayme: names of species
	  //printf("   Processing Row %d of Matrix %d", spp, newersites);
	  for(j =0; j<=i; j++){
		  if(j < i){
			  na = 0;
			  nb = 0;
			  nab = 0;
			  for(k = 0;k < newersites; k++){
				  if(bootmatrix[i][k] == 1 && bootmatrix[j][k]!=2)  na++;
				  if(bootmatrix[j][k] == 1 && bootmatrix[i][k]!=2)  nb++;
				  if(bootmatrix[i][k] == 1 && bootmatrix[j][k]==1)  nab++;
			  }
			  //printf("i = %d, j = %d, na = %d, nb = %d, nab = %d\n", i,j,na,nb,nab);
			  
              if((int)nab == 0){
				  distance[i][j] = 5.0;
				  fprintf(outfile, "\t%.7lf", distance[i][j]);
			  }else{
				  doa = -log(nab*1.0/na);
				  dob = -log(nab*1.0/nb);
				  dab = sqrt(doa*dob);
				  
				  doa = (doa <= 0) ? 0.0 : doa;
				  dob = (dob <= 0) ? 0.0 : dob;
				  dab = (dab <= 0) ? 0.0 : dab;
				  
				  distance[i][j] = dab;
				  fprintf(outfile, "\t%.7f", distance[i][j]);
				  //printf("doa = %.10f, dob = %.10f, dab = %.10f\n", doa,dob,dab);
		  
			  }
		  }else{
			  distance[i][j] = 0.0;
			  fprintf(outfile, "\t%.1lf", distance[i][j]);
		  }		  
	  }
	  fprintf(outfile, "\n");
  }
}  /* writebootdist */

void bootwrite()
{ /* does bootstrapping, distance computing and writes out data sets */
  long i, j, rr, repdiv10;

  repdiv10 = reps / 10;
  if (repdiv10 < 1)
    repdiv10 = 1;
  putchar('\n');
  firstrep = true;
  
  printf("Bootstrapping and distance computing. Please wait patiently...\n\n\n");
  
  for (rr = 1; rr <= (reps); rr++) {
    bootweights();
    for (i = 0; i < spp; i++)
      for (j = 0; j < newergroups; j++)
        charorder[i][j] = j;
    
    writebootdist();
	
    if (reps < 10 || rr % repdiv10 == 0) {
      printf("Progress status: completed bootstrapping and distance computing number: %4ld\n", rr);
#ifdef WIN32

#endif
      firstrep = false;
    }
  }
  printf("\nOutput written to file \"%s\"\n\n", outfilename);
}  /* bootwrite */

int main(int argc, Char *argv[])
{  /* Read in sequences or frequencies and bootstrap them */
#ifdef MAC
  argc = 1;                /* macsetup("bootdist","");                */
  argv[0] = "bootdist";
 
  
#endif
  init(argc,argv);
  //  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&infile, argv[1], "input file", "r", argv[0], infilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinput(argc, argv);
  bootwrite();
  
  freenewer();
  freenew();
  freerest();

  if (nodep)
    matrix_char_delete(nodep, spp);
  if (nodef)
    matrix_double_delete(nodef, spp);

  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  printf("Finish.\n\n");
  return 0;
}


////////////////////////  functions from phylip.c ///////////////////////


static void crash_handler(int sig_num)
{ /* when we crash, lets print out something usefull */
  printf("ERROR:  ");
  switch(sig_num) {
#ifdef SIGSEGV
    case SIGSEGV:
      puts("This program has caused a Segmentation fault.");
      break;
#endif /* SIGSEGV */
#ifdef SIGFPE
    case SIGFPE:
      puts("This program has caused a Floating Point Exception");
      break;
#endif  /* SIGFPE */
#ifdef SIGILL
    case SIGILL:
      puts("This program has attempted an illegal instruction");
      break;
#endif  /* SIGILL */
#ifdef SIGPIPE 
    case SIGPIPE:
      puts("This program tried to write to a broken pipe");
      break;
#endif  /* SIGPIPE */
#ifdef SIGBUS
    case SIGBUS:
      puts("This program had a bus error");
      break;
#endif /* SIGBUS */
  }   
  if (sig_num == SIGSEGV) {
    puts("       This may have been caused by an incorrectly formatted input file");
    puts("       or input tree file.  You should check those files carefully.");
  }
  abort();
}


void init(int argc, char** argv) 
{
#ifdef SIGSEGV
  signal(SIGSEGV, crash_handler);
#endif /* SIGSEGV */
#ifdef SIGFPE
  signal(SIGFPE, crash_handler);
#endif /* SIGFPE */
#ifdef SIGILL
  signal(SIGILL, crash_handler);
#endif /* SIGILL */
#ifdef SIGPIPE
  signal(SIGPIPE, crash_handler);
#endif /* SIGPIPE */
#ifdef SIGBUS
  signal(SIGBUS, crash_handler);
#endif /* SIGBUS */

  /* Set default terminal characteristics */
  ibmpc = IBMCRT;
  ansi = ANSICRT;

  /* Clear the screen */
  //cleerhome();

}

void scan_eoln(FILE *f) 
{ /* Eat everything up to EOF or newline, including newline */
  char ch;
  while (!eoff(f) && !eoln(f)) 
    gettc(f);
  if (!eoff(f)) 
    ch = gettc(f);
}

boolean eoff(FILE *f)
{ /* Return true iff next getc() is EOF */
    int ch;
    if (feof(f)) 
      return true;
    ch = getc(f);
    if (ch == EOF) {
      ungetc(ch, f);
      return true;
    }
    ungetc(ch, f);
    return false;
}  /*eoff*/


boolean eoln(FILE *f)
{ /* Return true iff next getc() is EOL or EOF */
    register int ch;
    ch = getc(f);
    if (ch == EOF)
      return true;
    ungetc(ch, f);
    return ((ch == '\n') || (ch == '\r'));
}  /*eoln*/


int filexists(char *filename)
{ /* Return true iff file already exists */
  FILE *fp;
  fp = fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}  /*filexists*/


const char* get_command_name (const char *vektor)
{ /* returns the name of the program from vektor without the whole path */
  char *last_slash;

  /* Point to the last slash... */
  last_slash = strrchr (vektor, DELIMITER);

  if (last_slash)
    /* If there was a last slash, return the character after it */
    return last_slash + 1;
  else
    /* If not, return the vector */
    return vektor;

}  /* get_command_name */


void EOF_error()
{ /* Print a message and exit when EOF is reached prematurely. */
  puts("\n\nERROR: Unexpected end-of-file.\n");
  exxit(-1);
}  /* EOF_error */


void getstryng(char *fname)
{ /* read in a file name from stdin and take off newline if any */
  char *end;

  fflush(stdout);
  fname = fgets(fname, FNMLNGTH, stdin);
  if ( fname == NULL )
    EOF_error();

  if ( (end = strpbrk(fname, "\n\r")) != NULL)
    *end = '\0';
    
} /* getstryng */


void countup(long *loopcount, long maxcount)
{ /* count how many times this loop has tried to read data, bail out
     if exceeds maxcount */

  (*loopcount)++;
  if ((*loopcount) >= maxcount) {
    printf("\nERROR: Made %ld attempts to read input in loop. Aborting run.\n",
            *loopcount);
    exxit(-1);
  }
} /* countup */


void openfile(FILE **fp,const char *filename,const char *filedesc,
              const char *mode,const char *application, char *perm)
{ /* open a file, testing whether it exists etc. */
  FILE *of;
  char file[FNMLNGTH];
  char filemode[3];
  const char *progname_without_path;
  long loopcount, loopcount2;
#if defined(OSX_CARBON) && defined(__MWERKS__)
  ProcessSerialNumber myProcess;
  FSRef bundleLocation;
  unsigned char bundlePath[FNMLNGTH];
  
  if(!fixedpath){
    /* change path to the bundle location instead of root directory */
    GetCurrentProcess(&myProcess);
    GetProcessBundleLocation(&myProcess, &bundleLocation);
    FSRefMakePath(&bundleLocation, bundlePath, FNMLNGTH);
    chdir((const char*)bundlePath);
    chdir(".."); /* get out of the .app directory */
    
    fixedpath = true;
  }
#endif

  progname_without_path = get_command_name(application);

  strcpy(file, filename);
  strcpy(filemode, mode);
  loopcount = 0;
  while (1){
#if ! OVERWRITE_FILES
#endif /* ! OVERWRITE_FILES */
    of = fopen(file,filemode);
    if (of)
      break;
    else {
      switch (filemode[0]){

      case 'r':
        printf("%s: can't find %s \"%s\"\n", progname_without_path,
            filedesc, file);
        file[0] = '\0';
        loopcount2 = 0;
        while ( file[0] =='\0' ) {
          printf("Please enter a new file name> ");
          fflush(stdout);
          countup(&loopcount2, 10);
          getstryng(file);
        }
        break;

      case 'w':
      case 'a':
        printf("%s: can't write %s \"%s\"\n", progname_without_path,
            filedesc, file);
        file[0] = '\0';
        loopcount2 = 0;
        while (file[0] =='\0') {
          printf("Please enter a new file name> ");
          fflush(stdout);
          countup(&loopcount2, 10);
          getstryng(file);
        }
        continue;

      default:
     printf("There is some error in the call of openfile. Unknown mode.\n");
        exxit(-1);
      }
    }
    countup(&loopcount, 20);
  }
  *fp = of;
  if (perm != NULL)
    strcpy(perm,file);
} /* openfile */


void cleerhome()
{ /* home cursor and clear screen, if possible */
#ifdef WIN32
  if(ibmpc || ansi){
  } else {
    printf("\n\n");
  }
#else
  printf("%s", ((ibmpc || ansi) ? ("\033[2J\033[H") : "\n\n"));
#endif
} /* cleerhome */


double randum(longer seed)
{ /* random number generator -- slow but machine independent
  This is a multiplicative congruential 32-bit generator
  x(t+1) = 1664525 * x(t) mod 2^32, one that passes the
  Coveyou-Macpherson and Lehmer tests, see Knuth ACP vol. 2
  We here implement it representing each integer in base-64
  notation -- i.e. as an array of 6 six-bit chunks   */
  
  long i, j, k, sum;
  longer mult, newseed;
  double x;

  mult[0] = 13;   /* these four statements set the multiplier */
  mult[1] = 24;   /* -- they are its "digits" in a base-64    */
  mult[2] = 22;   /*    notation: 1664525 = 6*64^3+22*64^2    */
  mult[3] = 6;    /*                         +24*64+13        */
  for (i = 0; i <= 5; i++)
    newseed[i] = 0;
  for (i = 0; i <= 5; i++) {  /* do the multiplication piecewise */
    sum = newseed[i];
    k = i;
    if (i > 3)
      k = 3;
    for (j = 0; j <= k; j++)
      sum += mult[j] * seed[i - j];
    newseed[i] = sum;
    for (j = i; j <= 4; j++) {
      newseed[j + 1] += newseed[j] / 64;
      newseed[j] &= 63;
    }
  }
  memcpy(seed, newseed, sizeof(longer));        /* new seed replaces old one */
  seed[5] &= 3;          /* from the new seed, get a floating point fraction */
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}  /* randum */


void randumize(longer seed, long *enterorder)
{ /* randomize input order of species -- randomly permute array enterorder */
  long i, j, k;
  
  for (i = 0; i < spp; i++) {
    j = (long)(randum(seed) * (i+1));
    k = enterorder[j];
    enterorder[j] = enterorder[i];
    enterorder[i] = k;
  }
} /* randumize */


double normrand(longer seed)
{/* standardized Normal random variate */
  double x;
  x = randum(seed)+randum(seed)+randum(seed)+randum(seed)
       + randum(seed)+randum(seed)+randum(seed)+randum(seed)
       + randum(seed)+randum(seed)+randum(seed)+randum(seed)-6.0;
  return(x);
} /* normrand */ 


long readlong(const char *prompt)
{ /* read a long */
  long res, loopcount;
  char string[100];

  loopcount = 0;
  do {
    printf("%s", prompt);
    fflush(stdout);
    getstryng(string);
    if (sscanf(string,"%ld",&res) == 1)
      break;
    countup(&loopcount, 10);
   } while (1);
  return res;
}  /* readlong */


void uppercase(Char *ch)
{ /* convert ch to upper case */
  *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */



void initseed(long *inseed, long *inseed0, longer seed)
{ /* input random number seed */
  long i;

  assert(inseed);
  assert(inseed0);
  
  *inseed = rand()%500 * 2 + 1;
  
  *inseed0 = *inseed;
  for (i = 0; i <= 5; i++)
    seed[i] = 0;
  i = 0;
  do {
    seed[i] = *inseed & 63;
    *inseed /= 64;
    i++;
  } while (*inseed != 0);
}  /*initseed*/


void initjumble(long *inseed, long *inseed0, longer seed, long *njumble)
{ /* input number of jumblings for jumble option */
  long loopcount;
  initseed(inseed, inseed0, seed);
  loopcount = 0;
  for (;;) {
    printf("Number of times to jumble?\n");
    fflush(stdout);
    if (scanf("%ld%*[^\n]", njumble) == 1) {
      getchar();
      if (*njumble >= 1)
        break;
    }
    countup(&loopcount, 10);
  }
}  /*initjumble*/


void initoutgroup(long *outgrno, long spp)
{ /* input outgroup number */
  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("Type number of the outgroup:\n");
    fflush(stdout);
    if (scanf("%ld%*[^\n]", outgrno) == 1) {
      getchar();
      if (*outgrno >= 1 && *outgrno <= spp)
        break;
      else {
        printf("BAD OUTGROUP NUMBER: %ld\n", *outgrno);
        printf("  Must be in range 1 - %ld\n", spp);
      }
    }
    countup(&loopcount, 10);
  }
}  /*initoutgroup*/


void initthreshold(double *threshold)
{ /* input threshold for threshold parsimony option */
  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("What will be the threshold value?\n");
    fflush(stdout);
    if (scanf("%lf%*[^\n]", threshold) == 1) {
      getchar();
      if (*threshold >= 1.0)
        break;
      else
        printf("BAD THRESHOLD VALUE:  it must be greater than 1\n");
    }
    countup(&loopcount, 10);
  }
  *threshold = (long)(*threshold * 10.0 + 0.5) / 10.0;
}  /*initthreshold*/


void initcatn(long *categs)
{ /* initialize category number for rate categories */
  long loopcount;

  loopcount = 0;
  *categs = 0;
  for (;;) {
    printf("Number of categories (1-%d)?\n", maxcategs);
    fflush(stdout);
    if (scanf("%ld%*[^\n]", categs) == 1) {
      getchar();
      if (*categs > maxcategs || *categs < 1)
        continue;
      else
        break;
    }
    countup(&loopcount, 10);
  }
}  /*initcatn*/


void initcategs(long categs, double *rate)
{ /* initialize category rates for HMM rates */
  long i, loopcount, scanned;
  char line[100], rest[100];
  boolean done;

  loopcount = 0;
  for (;;){
    printf("Rate for each category? (use a space to separate)\n");
    fflush(stdout);
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++){
      scanned = sscanf(line,"%lf %[^\n]", &rate[i], rest);
      if ((scanned < 2 && i < (categs - 1))
          || (scanned < 1 && i == (categs - 1)))
      {
        printf("Please enter exactly %ld values.\n", categs);
        done = false;
        break;
      }
      strcpy(line, rest);
    }
    if (done)
      break;
    countup(&loopcount, 100);
  }
}  /*initcategs*/


void initprobcat(long categs, double *probsum, double *probcat)
{ /* input probabilities of rate categores for HMM rates */
  long i, loopcount, scanned;
  boolean done;
  char line[100], rest[100];

  loopcount = 0;
  do {
    printf("Probability for each category?");
    printf(" (use a space to separate)\n");
    fflush(stdout);
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++) {
      scanned = sscanf(line, "%lf %[^\n]", &probcat[i], rest);
      if ((scanned < 2 && i < (categs - 1)) ||
          (scanned < 1 && i == (categs - 1))) {
        done = false;
        printf("Please enter exactly %ld values.\n", categs);
        break;
      }
      strcpy(line, rest);
    }
    if (!done)
      continue;
    *probsum = 0.0;
    for (i = 0; i < categs; i++)
      *probsum += probcat[i];
    if (fabs(1.0 - (*probsum)) > 0.001) {
      done = false;
      printf("Probabilities must add up to");
      printf(" 1.0, plus or minus 0.001.\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /*initprobcat*/


void lgr(long m, double b, raterootarray lgroot)
{ /* For use by initgammacat.  Get roots of m-th Generalized Laguerre
     polynomial, given roots of (m-1)-th, these are to be
     stored in lgroot[m][] */
  long i;
  double upper, lower, x, y;
  boolean dwn;   /* is function declining in this interval? */

  if (m == 1) {
    lgroot[1][1] = 1.0+b;
  } else {
    dwn = true;
    for (i=1; i<=m; i++) {
      if (i < m) {
        if (i == 1)
          lower = 0.0;
        else
          lower = lgroot[m-1][i-1];
        upper = lgroot[m-1][i];
      } else {                 /* i == m, must search above */
        lower = lgroot[m-1][i-1];
        x = lgroot[m-1][m-1];
        do {
          x = 2.0*x;
          y = glaguerre(m, b, x);
        } while ((dwn && (y > 0.0)) || ((!dwn) && (y < 0.0)));
        upper = x;
      }
      while (upper-lower > 0.000000001) {
        x = (upper+lower)/2.0;
        if (glaguerre(m, b, x) > 0.0) {
          if (dwn)
            lower = x;
          else
            upper = x;
        } else {
          if (dwn)
            upper = x;
          else
            lower = x;
        }        
      }
      lgroot[m][i] = (lower+upper)/2.0;
      dwn = !dwn;                /* switch for next one */
    }
  }
} /* lgr */


double logfac (long n)
{ /* log(n!) values were calculated with Mathematica
     with a precision of 30 digits */
  long i;
  double x;

  switch (n)
    {
    case 0:
      return 0.;
    case 1:
      return 0.;
    case 2:
      return 0.693147180559945309417232121458;
    case 3:
      return 1.791759469228055000812477358381;
    case 4:
      return 3.1780538303479456196469416013;
    case 5:
      return 4.78749174278204599424770093452;
    case 6:
      return 6.5792512120101009950601782929;
    case 7:
      return 8.52516136106541430016553103635;
    case 8:
      return 10.60460290274525022841722740072;
    case 9:
      return 12.80182748008146961120771787457;
    case 10:
      return 15.10441257307551529522570932925;
    case 11:
      return 17.50230784587388583928765290722;
    case 12:
      return 19.98721449566188614951736238706;
    default:
      x = 19.98721449566188614951736238706;
      for (i = 13; i <= n; i++)
        x += log(i);
      return x;
    }
} /* logfac */


double glaguerre(long m, double b, double x)
{ /* Generalized Laguerre polynomial computed recursively.
     For use by initgammacat */
  long i;
  double gln, glnm1, glnp1; /* L_n, L_(n-1), L_(n+1) */

  if (m == 0)
    return 1.0;
  else {
    if (m == 1)
      return 1.0 + b - x;
    else {
      gln = 1.0+b-x;
      glnm1 = 1.0;
      for (i=2; i <= m; i++) {
        glnp1 = ((2*(i-1)+b+1.0-x)*gln - (i-1+b)*glnm1)/i;
        glnm1 = gln;
        gln = glnp1;
      }
      return gln;
    }
  }
} /* glaguerre */


void initlaguerrecat(long categs, double alpha, double *rate, double *probcat)
{ /* calculate rates and probabilities to approximate Gamma distribution
     of rates with "categs" categories and shape parameter "alpha" using
     rates and weights from Generalized Laguerre quadrature */
  long i;
  raterootarray lgroot; /* roots of GLaguerre polynomials */
  double f, x, xi, y;

  alpha = alpha - 1.0;
  lgroot[1][1] = 1.0+alpha;
  for (i = 2; i <= categs; i++)
    lgr(i, alpha, lgroot);                   /* get roots for L^(a)_n */
  /* here get weights */
  /* Gamma weights are (1+a)(1+a/2) ... (1+a/n)*x_i/((n+1)^2 [L_{n+1}^a(x_i)]^2)  */
  f = 1;
  for (i = 1; i <= categs; i++)
    f *= (1.0+alpha/i);
  for (i = 1; i <= categs; i++) {
    xi = lgroot[categs][i];
    y = glaguerre(categs+1, alpha, xi);
    x = f*xi/((categs+1)*(categs+1)*y*y);
    rate[i-1] = xi/(1.0+alpha);
    probcat[i-1] = x;
  }
} /* initlaguerrecat */


double hermite(long n, double x)
{ /* calculates hermite polynomial with degree n and parameter x */
  /* seems to be unprecise for n>13 -> root finder does not converge*/
  double h1 = 1.;
  double h2 = 2. * x;
  double xx = 2. * x;
  long i;

  for (i = 1; i < n; i++) {
    xx = 2. * x * h2 - 2. * (i) * h1;
    h1 = h2;
    h2 = xx;
  }
  return xx;
} /* hermite */


void root_hermite(long n, double *hroot)
{ /* find roots of Hermite polynmials */
  long z;
  long ii;
  long start;

  if (n % 2 == 0) {
    start = n/2;
    z = 1;
  } else {
    start = n/2 + 1;
    z=2;
    hroot[start-1] = 0.0;
  }
  for (ii = start; ii < n; ii++) {         /* search only upwards*/
    hroot[ii] = halfroot(hermite, n, hroot[ii-1]+EPSILON, 1./n);
    hroot[start - z] = -hroot[ii];
    z++;
  }
} /* root_hermite */


double halfroot(double (*func)(long m, double x), long n, double startx,
                double delta)
{ /* searches from the bound (startx) only in one direction
     (by positive or negative delta, which results in
     other-bound=startx+delta)
     delta should be small.
     (*func) is a function with two arguments  */
  double xl;
  double xu;
  double xm = 0;
  double fu;
  double fl;
  double fm = 100000.;
  double gradient;
  boolean dwn = false; 

  /* decide if we search above or below startx and escapes to trace back
     to the starting point that most often will be
     the root from the previous calculation */
  if (delta < 0) {
    xu = startx;
    xl = xu + delta;
  } else {
    xl = startx;
    xu = xl + delta;
  }
  delta = fabs(delta);
  fu = (*func)(n, xu);
  fl = (*func)(n, xl);
  gradient = (fl-fu)/(xl-xu);
  while(fabs(fm) > EPSILON) {        /* is root outside of our bracket?*/
    if ((fu<0.0 && fl<0.0) || (fu>0.0 && fl > 0.0)) {
      xu += delta;
      fu = (*func)(n, xu);
      fl = (*func)(n, xl);
      gradient = (fl-fu)/(xl-xu);
      dwn = (gradient < 0.0) ? true : false;
    } else {
      xm = xl - fl / gradient;
      fm = (*func)(n, xm);
      if (dwn) {
        if (fm > 0.) {
          xl = xm;
          fl = fm;
        } else {
          xu = xm;
          fu = fm;
        }
      } else {
        if (fm > 0.) {
          xu = xm;
          fu = fm;
        } else {
          xl = xm;
          fl = fm;
        }
      }
      gradient = (fl-fu)/(xl-xu);
    }
  }
  return xm;
} /* halfroot */


void hermite_weight(long n, double * hroot, double * weights)
{
  /* calculate the weights for the hermite polynomial at the roots
     using formula from Abramowitz and Stegun chapter 25.4.46 p.890 */
  long i;
  double hr2;
  double numerator;

  numerator = exp(0.6931471805599 * ( n-1.) + logfac(n)) / (n*n);
  for (i = 0; i < n; i++) {
    hr2 = hermite(n-1, hroot[i]);
    weights[i] = numerator / (hr2*hr2);
  }
} /* hermiteweight */


void inithermitcat(long categs, double alpha, double *rate, double *probcat)
{ /* calculates rates and probabilities */
  long i;
  double *hroot;
  double std;

  std = SQRT2 /sqrt(alpha);
  hroot = (double *) Malloc((categs+1) * sizeof(double));
  root_hermite(categs, hroot);         /* calculate roots */
  hermite_weight(categs, hroot, probcat);  /* set weights */
  for (i=0; i<categs; i++) {           /* set rates */
    rate[i] = 1.0 + std * hroot[i];
    probcat[i] = probcat[i];
  }
  free(hroot);
} /* inithermitcat */


void initgammacat (long categs, double alpha, double *rate, double *probcat)
{ /* calculate rates and probabilities to approximate Gamma distribution
   of rates with "categs" categories and shape parameter "alpha" using
   rates and weights from Generalized Laguerre quadrature or from
   Hermite quadrature */

  if (alpha >= 100.0)
    inithermitcat(categs, alpha, rate, probcat);
  else
    initlaguerrecat(categs, alpha, rate, probcat);
} /* initgammacat */


void inithowmany(long *howmanny, long howoften)
{/* input how many cycles */
  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("How many cycles of %4ld trees?\n", howoften);
    fflush(stdout);
    if (scanf("%ld%*[^\n]", howmanny) == 1) {
      getchar();
      if (*howmanny >= 1)
        break;
    }
    countup(&loopcount, 10);
  }
}  /*inithowmany*/


void inithowoften(long *howoften)
{ /* input how many trees per cycle */
  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("How many trees per cycle?\n");
    fflush(stdout);
    if (scanf("%ld%*[^\n]", howoften) == 1) {
      getchar();
      if (*howoften >= 1)
        break;
    }
    countup(&loopcount, 10);
  }
}  /*inithowoften*/


void initlambda(double *lambda)
{ /* input patch length parameter for autocorrelated HMM rates */

  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("Mean block length of sites having the same rate (greater than 1)?\n");
    fflush(stdout);
    if (scanf("%lf%*[^\n]", lambda) == 1) {
      getchar();
      if (*lambda > 1.0)
        break;
    }
    countup(&loopcount, 10);
  }
  *lambda = 1.0 / *lambda;
} /* initlambda */ 



void initfreqs(double *freqa, double *freqc, double *freqg, double *freqt)
{ /* input frequencies of the four bases */
  char input[100];
  long scanned, loopcount;

  printf("Base frequencies for A, C, G, T/U (use blanks to separate)?\n");
  loopcount = 0;
  do {
    fflush(stdout);
    getstryng(input);
    scanned = sscanf(input,"%lf%lf%lf%lf%*[^\n]", freqa, freqc, freqg, freqt);
    if (scanned == 4)
      break;
    else
      printf("Please enter exactly 4 values.\n");
    countup(&loopcount, 100);
  } while (1);
}  /* initfreqs */


void initratio(double *ttratio)
{ /* input transition/transversion ratio */
  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("Transition/transversion ratio?\n");
    fflush(stdout);
    if (scanf("%lf%*[^\n]", ttratio) == 1) {
      getchar();
      if (*ttratio >= 0.0)
        break;
      else
        printf("Transition/transversion ratio cannot be negative.\n");
    }
    countup(&loopcount, 10);
  }
}  /* initratio */


void initpower(double *power)
{
  for (;;) {
    printf("New power?\n");
    fflush(stdout);
    if (scanf("%lf%*[^\n]", power) == 1) {
      getchar();
      break;
    }
  }
}  /*initpower*/


void initdatasets(long *datasets)
{
  /* handle multi-data set option */
  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("How many data sets?\n");
    fflush(stdout);
    if (scanf("%ld%*[^\n]", datasets) == 1) {
      getchar();
      if (*datasets > 1)
        break;
      else
        printf("Bad data sets number:  it must be greater than 1\n");
    }
    countup(&loopcount, 10);
  }
} /* initdatasets */


void justweights(long *datasets)
{
  /* handle multi-data set option by weights */
  long loopcount;

  loopcount = 0;
  for (;;) {
    printf("How many sets of weights?\n");
    fflush(stdout);
    if (scanf("%ld%*[^\n]", datasets) == 1) {
      getchar();
      if (*datasets >= 1)
        break;
      else 
        printf("BAD NUMBER:  it must be greater than 1\n");
    }
    countup(&loopcount, 10);
  }
} /* justweights */


void initterminal(boolean *ibmpc, boolean *ansi)
{
  /* handle terminal option */
  if (*ibmpc) {
    *ibmpc = false;
    *ansi = true;
  } else if (*ansi)
      *ansi = false;
    else
      *ibmpc = true;
}  /*initterminal*/


void initnumlines(long *screenlines)
{
  long loopcount;

  loopcount = 0;
  do {
    *screenlines = readlong("Number of lines on screen?\n");
    countup(&loopcount, 10);
  } while (*screenlines <= 12);
}  /*initnumlines*/


void initbestrees(bestelm *bestrees, long maxtrees, boolean glob)
{
  /* initializes either global or local field of each array in bestrees */
  long i;

  if (glob)
    for (i = 0; i < maxtrees; i++)
      bestrees[i].gloreange = false;
  else
    for (i = 0; i < maxtrees; i++)
      bestrees[i].locreange = false;
} /* initbestrees */


void newline(FILE *filename, long i, long j, long k)
{
  /* go to new line if i is a multiple of j, indent k spaces */
  long m;

  if ((i - 1) % j != 0 || i <= 1)
    return;
  putc('\n', filename);
  for (m = 1; m <= k; m++)
    putc(' ', filename);
}  /* newline */


void inputnumbersold(long *spp, long *chars, long *nonodes, long n)
{
  /* input the numbers of species and of characters */

  if (fscanf(infile, "%ld%ld", spp, chars) != 2 || *spp <= 0 || *chars <= 0) {
    printf(
    "ERROR: Unable to read the number of species or characters in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  *nonodes = *spp * 2 - n;
}  /* inputnumbersold */


void inputnumbers(long *spp, long *chars, long *nonodes, long n)
{
  /* Read numbers of species and characters from first line of a data set.
   * Return the results in *spp and *chars, respectively. Also returns
   * (*spp * 2 - n)  in *nonodes */

  if (fscanf(infile, "%ld%ld", spp, chars) != 2 || *spp <= 0 || *chars <= 0) {
    printf(
    "ERROR: Unable to read the number of species or characters in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  *nonodes = *spp * 2 - n;
}  /* inputnumbers */


void inputnumbers2(long *spp, long *nonodes, long n)
{
  /* read species number */

  if (fscanf(infile, "%ld", spp) != 1 || *spp <= 0) {
    printf("ERROR: Unable to read the number of species in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  fprintf(outfile, "\n%4ld Populations\n", *spp);
  *nonodes = *spp * 2 - n;
}  /* inputnumbers2 */


void inputnumbers3(long *spp, long *chars)
{
  /* input the numbers of species and of characters */

  if (fscanf(infile, "%ld%ld", spp, chars) != 2 || *spp <= 0 || *chars <= 0) {
    printf(
    "ERROR: Unable to read the number of species or characters in data set\n");
    printf(
       "The input file is incorrect (perhaps it was not saved text only).\n");
    exxit(-1);
  }
}  /* inputnumbers3 */


void samenumsp(long *chars, long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp, curchs;

  if (eoln(infile)) 
    scan_eoln(infile);
  if (fscanf(infile, "%ld%ld", &cursp, &curchs) == 2) {
    if (cursp != spp) {
      printf(
           "\n\nERROR: Inconsistent number of species in data set %ld\n\n", ith);
      exxit(-1);
    }
  }
  else {
    printf(
        "Unable to read number of species and sites from data set %ld\n\n", ith);
    exxit(-1);
  }
  *chars = curchs;
} /* samenumsp */


void samenumsp2(long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp;

  if (eoln(infile)) 
    scan_eoln(infile);
  if (fscanf(infile, "%ld", &cursp) != 1) {
    printf("\n\nERROR: Unable to read number of species in data set %ld\n",
      ith);
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
    exxit(-1);
  }
  if (cursp != spp) {
    printf(
      "\n\nERROR: Inconsistent number of species in data set %ld\n\n", ith);
    exxit(-1);
  }
} /* samenumsp2 */


void readoptions(long *extranum, const char *options)
{ /* read option characters from input file */
  Char ch;

  while (!(eoln(infile))) {
    ch = gettc(infile);
    uppercase(&ch);
    if (strchr(options, ch) != NULL)
     (* extranum)++;
    else if (!(ch == ' ' || ch == '\t')) {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(infile);
}  /* readoptions */


void matchoptions(Char *ch, const char *options)
{  /* match option characters to those in auxiliary options line */

  *ch = gettc(infile);
  uppercase(ch);
  if (strchr(options, *ch) == NULL) {
    printf("ERROR: Incorrect auxiliary options line");
    printf(" which starts with %c\n", *ch);
    exxit(-1);
  }
}  /* matchoptions */


void inputweightsold(long chars, steptr weight, boolean *weights)
{
  Char ch;
  int i;
  for (i = 1; i < nmlngth ; i++)
    getc(infile);
 
  for (i = 0; i < chars; i++) {
    do {
      if (eoln(infile)) 
        scan_eoln(infile);
      ch = gettc(infile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    weight[i] = 1;
    if (isdigit(ch))
      weight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      weight[i] = ch - 'A' + 10;
    } else {
      printf("\n\nERROR: Bad weight character: %c\n\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(infile);
  *weights = true;
} /*inputweightsold*/


void inputweights(long chars, steptr weight, boolean *weights)
{
  /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
  Char ch;
  long i;

  for (i = 0; i < chars; i++) {
    do {
      if (eoln(weightfile)) 
        scan_eoln(weightfile);
      ch = gettc(weightfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    weight[i] = 1;
    if (isdigit(ch))
      weight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      weight[i] = ch - 'A' + 10;
    } else {
      printf("\n\nERROR: Bad weight character: %c\n\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(weightfile);
  *weights = true;
}  /* inputweights */


void inputweights2(long a, long b, long *weightsum,
        steptr weight, boolean *weights, const char *prog)
{
  /* input the character weights, 0 or 1 */
  Char ch;
  long i;

  *weightsum = 0;
  for (i = a; i < b; i++) {
    do {
      if (eoln(weightfile))
        scan_eoln(weightfile);
      ch = gettc(weightfile);
    } while (ch == ' ');
    weight[i] = 1;
    if (ch == '0' || ch == '1')
      weight[i] = ch - '0';
    else {
      printf("\n\nERROR: Bad weight character: %c -- ", ch);
      printf("weights in %s must be 0 or 1\n", prog);
      exxit(-1);
    }
    *weightsum += weight[i];
  }
  *weights = true;
  scan_eoln(weightfile);
}  /* inputweights2 */


void printweights(FILE *filename, long inc, long chars,
        steptr weight, const char *letters)
{
  /* print out the weights of sites */
  long i, j;
  boolean letterweights;

  letterweights = false;
  for (i = 0; i < chars; i++)
    if (weight[i] > 9)
      letterweights = true;
  fprintf(filename, "\n    %s are weighted as follows:", letters);
  if (letterweights)
    fprintf(filename, " (A = 10, B = 11, etc.)\n");
  else
    putc('\n', filename);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    if (weight[i+inc] < 10)
      fprintf(filename, "%ld", weight[i + inc]);
    else
      fprintf(filename, "%c", 'A'-10+(int)weight[i + inc]);
    if ((i+1) % 5 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printweights */


void inputcategs(long a, long b, steptr category, long categs, const char *prog)
{
  /* input the categories, 1-9 */
  Char ch;
  long i;

  for (i = a; i < b; i++) {
    do {
      if (eoln(catfile)) 
        scan_eoln(catfile);
      ch = gettc(catfile);
    } while (ch == ' ');
    if ((ch >= '1') && (ch <= ('0'+categs)))
      category[i] = ch - '0';
    else {
      printf("\n\nERROR: Bad category character: %c", ch);
      printf(" -- categories in %s are currently 1-%ld\n", prog, categs);
      exxit(-1);
    }
  }
  scan_eoln(catfile);
}  /* inputcategs */


void printcategs(FILE *filename, long chars, steptr category,
                 const char *letters)
{
  /* print out the sitewise categories */
  long i, j;

  fprintf(filename, "\n    %s are:\n", letters);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    fprintf(filename, "%ld", category[i]);
    if ((i+1) % 10 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printcategs */


void inputfactors(long chars, Char *factor, boolean *factors)
{
  /* reads the factor symbols */
  long i;

  for (i = 0; i < (chars); i++) {
    if (eoln(factfile)) 
      scan_eoln(factfile);
    factor[i] = gettc(factfile);
    if (factor[i] == '\n')
      factor[i] = ' ';
  }
  scan_eoln(factfile);
  *factors = true;
}  /* inputfactors */


void printfactors(FILE *filename, long chars, Char *factor, const char *letters)
{
  /* print out list of factor symbols */
  long i;

  fprintf(filename, "Factors%s:\n\n", letters);
  for (i = 1; i <= nmlngth - 5; i++)
    putc(' ', filename);
  for (i = 1; i <= (chars); i++) {
    newline(filename, i, 55, nmlngth + 3);
    putc(factor[i - 1], filename);
    if (i % 5 == 0)
      putc(' ', filename);
  }
  putc('\n', filename);
}  /* printfactors */


void headings(long chars, const char *letters1, const char *letters2)
{
  long i, j;

  putc('\n', outfile);
  j = nmlngth + (chars + (chars - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  fprintf(outfile, "Name");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n", letters1);
  fprintf(outfile, "----");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n\n", letters2);
}  /* headings */


void initname(long i)
{
  /* read in species name */
  long j;

  for (j = 0; j < nmlngth; j++) {
    if (eoff(infile) | eoln(infile)){
      printf("\n\nERROR: end-of-line or end-of-file");
      printf(" in the middle of species name for species %ld\n\n", i+1);
      exxit(-1);
    }
    nayme[i][j] = gettc(infile);
    if ((nayme[i][j] == '(') || (nayme[i][j] == ')') || (nayme[i][j] == ':')
        || (nayme[i][j] == ',') || (nayme[i][j] == ';') || (nayme[i][j] == '[')
        || (nayme[i][j] == ']')) {
      printf("\nERROR: Species name may not contain characters ( ) : ; , [ ] \n");
      printf("       In name of species number %ld there is character %c\n\n",
              i+1, nayme[i][j]);
      exxit(-1);
    }
  }
} /* initname */


void findtree(boolean *found, long *pos, long nextree, long *place, bestelm *bestrees)
{
  /* finds tree given by array place in array bestrees by binary search */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i, lower, upper;
  boolean below, done;
  
  below = false;
  lower = 1;
  upper = nextree - 1;
  (*found) = false;
  while (!(*found) && lower <= upper) {
    (*pos) = (lower + upper) / 2;
    i = 3;
    done = false;
    while (!done) {
      done = (i > spp);
      if (!done)
        done = (place[i - 1] != bestrees[(*pos) - 1].btree[i - 1]);
      if (!done)
        i++;
    }
    (*found) = (i > spp);
    if (*found)
      break;
    below = (place[i - 1] <  bestrees[(*pos )- 1].btree[i - 1]);
    if (below)
      upper = (*pos) - 1;
    else
      lower = (*pos) + 1;
  }
  if (!(*found) && !below)
    (*pos)++;
}  /* findtree */


void addtree(long pos, long *nextree, boolean collapse, long *place, bestelm *bestrees)
{
  /* puts tree from array place in its proper position in array bestrees */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i;
  
  for (i = *nextree - 1; i >= pos; i--){
    memcpy(bestrees[i].btree, bestrees[i - 1].btree, spp * sizeof(long));
    bestrees[i].gloreange = bestrees[i - 1].gloreange;
    bestrees[i - 1].gloreange = false;
    bestrees[i].locreange = bestrees[i - 1].locreange;
    bestrees[i - 1].locreange = false;
    bestrees[i].collapse = bestrees[i - 1].collapse;
  }
  for (i = 0; i < spp; i++)
    bestrees[pos - 1].btree[i] = place[i];
  bestrees[pos - 1].collapse = collapse;
  (*nextree)++;
}  /* addtree */


long findunrearranged(bestelm *bestrees, long nextree, boolean glob)
{
  /* finds bestree with either global or local field false */
  long i;

  if (glob) {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].gloreange)
        return i;
  } else {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].locreange)
        return i;
  }
  return -1;
} /* findunrearranged */


boolean torearrange(bestelm *bestrees, long nextree)
{ /* sees if any best tree is yet to be rearranged */
  
  if (findunrearranged(bestrees, nextree, true) >= 0)
    return true;
  else if (findunrearranged(bestrees, nextree, false) >= 0)
    return true;
  else
    return false;
} /* torearrange */


void reducebestrees(bestelm *bestrees, long *nextree)
{
  /* finds best trees with collapsible branches and deletes them */
  long i, j;
  i = 0;
  j = *nextree - 2;
  do {
    while (!bestrees[i].collapse && i < *nextree - 1) i++;
    while (bestrees[j].collapse && j >= 0) j--;
    if (i < j) {
      memcpy(bestrees[i].btree, bestrees[j].btree, spp * sizeof(long));
      bestrees[i].gloreange = bestrees[j].gloreange;
      bestrees[i].locreange = bestrees[j].locreange;
      bestrees[i].collapse = false;
      bestrees[j].collapse = true;
    }
  } while (i < j);
  *nextree = i + 1;
} /* reducebestrees */


void shellsort(double *a, long *b, long n)
{ /* Shell sort keeping a, b in same order */
  /* used by dnapenny, dolpenny, & penny */
  long gap, i, j, itemp;
  double rtemp;

  gap = n / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
     if (a[j - 1] > a[j + gap - 1]) {
       rtemp = a[j - 1];
       a[j - 1] = a[j + gap - 1];
       a[j + gap - 1] = rtemp;
       itemp = b[j - 1];
       b[j - 1] = b[j + gap - 1];
       b[j + gap - 1] = itemp;
     }
     j -= gap;
      }
    }
    gap /= 2;
  }
}  /* shellsort */


void getch(Char *c, long *parens, FILE *treefile)
{ /* get next nonblank character */

  do {
    if (eoln(treefile)) 
      scan_eoln(treefile);
    (*c) = gettc(treefile);

    if ((*c) == '\n' || (*c) == '\t')
      (*c) = ' ';
  } while ( *c == ' ' && !eoff(treefile) );
  if ((*c) == '(')
    (*parens)++;
  if ((*c) == ')')
    (*parens)--;
}  /* getch */


void getch2(Char *c, long *parens)
{ /* get next nonblank character */
  do {
    if (eoln(intree)) 
      scan_eoln(intree);
    *c = gettc(intree);
    if (*c == '\n' || *c == '\t')
      *c = ' ';
  } while (!(*c != ' ' || eoff(intree)));
  if (*c == '(')
   (*parens)++;
  if (*c == ')')
    (*parens)--;
}  /* getch2 */


void findch(Char c, Char *ch, long which)
{ /* scan forward until find character c */
  boolean done;
  long dummy_parens;
  done = false;
  while (!done) {
    if (c == ',') {
      if (*ch == '(' || *ch == ')' || *ch == ';') {
        printf(
      "\n\nERROR in user tree %ld: unmatched parenthesis or missing comma\n\n",
          which);
        exxit(-1);
      } else if (*ch == ',')
        done = true;
    } else if (c == ')') {
      if (*ch == '(' || *ch == ',' || *ch == ';') {
        printf("\n\nERROR in user tree %ld: ", which);
        printf("unmatched parenthesis or non-bifurcated node\n\n");
        exxit(-1);
      } else {
        if (*ch == ')')
          done = true;
      }
    } else if (c == ';') {
      if (*ch != ';') {
        printf("\n\nERROR in user tree %ld: ", which);
        printf("unmatched parenthesis or missing semicolon\n\n");
        exxit(-1);
      } else
        done = true;
    }
    if (*ch != ')' && done)
      continue;
   getch(ch, &dummy_parens, intree);
  }
}  /* findch */


void findch2(Char c, long *lparens, long *rparens, Char *ch)
{ /* skip forward in user tree until find character c */
  boolean done;
  long dummy_parens;
  done = false;

  while (!done) {
    if (c == ',') {
      if (*ch == '(' || *ch == ')' || *ch == ':' || *ch == ';') {
        printf("\n\nERROR in user tree: ");
        printf("unmatched parenthesis, missing comma");
        printf(" or non-trifurcated base\n\n");
     exxit(-1);
      } else if (*ch == ',')
        done = true;
    } else if (c == ')') {
      if (*ch == '(' || *ch == ',' || *ch == ':' || *ch == ';') {
        printf(
   "\n\nERROR in user tree: unmatched parenthesis or non-bifurcated node\n\n");
        exxit(-1);
      } else if (*ch == ')') {
        (*rparens)++;
        if ((*lparens) > 0 && (*lparens) == (*rparens)) {
          if ((*lparens) == spp - 2) {
           getch(ch, &dummy_parens, intree);
            if (*ch != ';') {
              printf( "\n\nERROR in user tree: ");
              printf("unmatched parenthesis or missing semicolon\n\n");
              exxit(-1);
            }
          }
        }
     done = true;
      }
    }
    if (*ch != ')' && done)
      continue;
    if (*ch == ')')
     getch(ch, &dummy_parens, intree);
  }
}  /* findch2 */

void processlength(double *valyew, double *divisor, Char *ch, 
        boolean *lengthIsNegative, FILE *treefile, long *parens)
{ /* read a branch length from a treefile */
  long digit, ordzero, exponent, exponentIsNegative;
  boolean pointread, hasExponent;

  ordzero = '0';
  *lengthIsNegative = false;
  pointread = false;
  hasExponent = false;
  exponentIsNegative = -1; // 3 states:  -1 = unassigned, 1 = true, 0 = false
  exponent = 0;
  *valyew = 0.0;
  *divisor = 1.0;
  getch(ch, parens, treefile);
  if ('+' == *ch)
    getch(ch, parens, treefile); // ignore leading '+', because "+1.2345" == "1.2345"
  else if ('-' == *ch)
    {
      *lengthIsNegative = true;
      getch(ch, parens, treefile);
    }
  digit = (long)(*ch - ordzero);
  while ( ((digit <= 9) && (digit >= 0)) || '.' == *ch || '-' == *ch
	  || '+' == *ch || 'E' == *ch || 'e' == *ch) {
    if ('.' == *ch)
      {
	if (!pointread)
	  pointread = true;
	else
	  {
	    printf("\n\nERROR: Branch length found with more than one \'.\' in it.\n\n");
	    exxit(-1);
	  }
      }
    else if ('+' == *ch)
      {
	if (hasExponent && -1 == exponentIsNegative)
	  exponentIsNegative = 0; // 3 states:  -1 = unassigned, 1 = true, 0 = false
	else
	  {
	    printf("\n\nERROR: Branch length found with \'+\' in an unexpected place.\n\n");
	    exxit(-1);
	  }
      }
    else if ('-' == *ch)
      {
	if (hasExponent && -1 == exponentIsNegative)
	  exponentIsNegative = 1; // 3 states:  -1 = unassigned, 1 = true, 0 = false
	else
	  {
	    printf("\n\nERROR: Branch length found with \'-\' in an unexpected place.\n\n");
	    exxit(-1);
	  }
      }
    else if ('E' == *ch || 'e' == *ch)
      {
	if (!hasExponent)
	  hasExponent = true;
	else
	  {
	    printf("\n\nERROR: Branch length found with more than one \'E\' in it.\n\n");
	    exxit(-1);
	  }
      }
    else {
      if (!hasExponent)
	{
	  *valyew = *valyew * 10.0 + digit;
	  if (pointread)
	    *divisor *= 10.0;
	}
      else
	exponent = 10*exponent + digit;
    }
    getch(ch, parens, treefile);
    digit = (long)(*ch - ordzero);
  }
  if (hasExponent)
    {
      if (exponentIsNegative)
	*divisor *= pow(10.,(double)exponent);
      else
	*divisor /= pow(10.,(double)exponent);
    }
  if (*lengthIsNegative)
    *valyew = -(*valyew);
}  /* processlength */

void writename(long start, long n, long *enterorder)
{ /* write species name and number in entry order */
  long i, j;

  for (i = start; i < start+n; i++) {
    printf(" %3ld. ", i+1);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[enterorder[i] - 1][j]);
    putchar('\n');
    fflush(stdout);
  }
}  /* writename */


void memerror()
{
  printf("Error allocating memory\n");
  exxit(-1);
}  /* memerror */


void odd_malloc(long x)
{ /* error message if attempt to malloc too little or too much memory */
  printf ("ERROR: a function asked for an inappropriate amount of memory:");
  printf ("  %ld bytes\n", x);
  printf ("       This can mean one of two things:\n");
  printf ("       1.  The input file is incorrect");
  printf (" (perhaps it was not saved as Text Only),\n");
  printf ("       2.  There is a bug in the program.\n");
  printf ("       Please check your input file carefully.\n");
  printf ("       If it seems to be a bug, please mail joe (at) gs.washington.edu\n");
  printf ("       with the name of the program, your computer system type,\n");
  printf ("       a full description of the problem, and with the input data file.\n");
  printf ("       (which should be in the body of the message, not as an Attachment).\n");

  /* abort() can be used to crash */
  
  exxit(-1);
}


MALLOCRETURN *mymalloc(long x)
{ /* wrapper for malloc, allowing error message if too little, too much */
  MALLOCRETURN *new_block;

  if ((x <= 0) ||
      (x > TOO_MUCH_MEMORY))
    odd_malloc(x);

  new_block = (MALLOCRETURN *)calloc(1, x);

  if (!new_block) {
    memerror();
    return (MALLOCRETURN *) new_block;
  } else
    return (MALLOCRETURN *) new_block;
} /* mymalloc */


void gnu(node **grbg, node **p)
{ /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
  } else
    *p = (node *)Malloc(sizeof(node));

  (*p)->back       = NULL;
  (*p)->next       = NULL;
  (*p)->tip        = false;
  (*p)->times_in_tree = 0.0;
  (*p)->r          = 0.0;
  (*p)->theta      = 0.0;
  (*p)->x          = NULL;
  (*p)->protx           = NULL;        /* for the sake of proml     */
}  /* gnu */


void chuck(node **grbg, node *p)
{
  /* collect garbage on p -- put it on front of garbage list */

  p->back = NULL;
  p->next = *grbg;
  *grbg = p;
}  /* chuck */


void zeronumnuc(node *p, long endsite)
{
  long i,j;

  for (i = 0; i < endsite; i++)
    for (j = (long)A; j <= (long)O; j++)
      p->numnuc[i][j] = 0;
} /* zeronumnuc */


void zerodiscnumnuc(node *p, long endsite)
{
  long i,j;

  for (i = 0; i < endsite; i++)
    for (j = (long)zero; j <= (long)seven; j++)
      p->discnumnuc[i][j] = 0;
} /* zerodiscnumnuc */


void allocnontip(node *p, long *zeros, long endsite)
{ /* allocate an interior node */
  /* used by dnacomp, dnapars, & dnapenny */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->base = (baseptr)Malloc(endsite*sizeof(long));
  p->oldbase = (baseptr)Malloc(endsite*sizeof(long));
  p->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  memcpy(p->base, zeros, endsite*sizeof(long));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->oldbase, zeros, endsite*sizeof(long));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
  zeronumnuc(p, endsite);
}  /* allocnontip */


void allocdiscnontip(node *p, long *zeros, unsigned char *zeros2, long endsite)
{ /* allocate an interior node */
  /* used by pars */
 
  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->discbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  p->olddiscbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  p->discnumnuc = (discnucarray *)Malloc(endsite*sizeof(discnucarray));
  memcpy(p->discbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->olddiscbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
  zerodiscnumnuc(p, endsite);
}  /* allocdiscnontip */


void allocnode(node **anode, long *zeros, long endsite)
{ /* allocate a node */
  /* used by dnacomp, dnapars, & dnapenny */
  *anode = (node *)Malloc(sizeof(node));
  allocnontip(*anode, zeros, endsite);
}  /* allocnode */


void allocdiscnode(node **anode, long *zeros, unsigned char *zeros2, 
        long endsite)
{ /* allocate a node */
  /* used by pars */
  *anode = (node *)Malloc(sizeof(node));
  allocdiscnontip(*anode, zeros, zeros2, endsite);
}  /* allocdiscnontip */


void gnutreenode(node **grbg, node **p, long i, long endsite, long *zeros)
{ /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
    memcpy((*p)->numsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->base, zeros, endsite*sizeof(long));
    memcpy((*p)->oldbase, zeros, endsite*sizeof(long));
    zeronumnuc(*p, endsite);
  } else
    allocnode(p, zeros, endsite);
  (*p)->back = NULL;
  (*p)->next = NULL;
  (*p)->tip = false;
  (*p)->visited = false;
  (*p)->index = i;
  (*p)->numdesc = 0;
  (*p)->sumsteps = 0.0;
}  /* gnutreenode */


void gnudisctreenode(node **grbg, node **p, long i,
        long endsite, long *zeros, unsigned char *zeros2)
{ /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
    memcpy((*p)->numsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->discbase, zeros2, endsite*sizeof(unsigned char));
    memcpy((*p)->olddiscbase, zeros2, endsite*sizeof(unsigned char));
    zerodiscnumnuc(*p, endsite);
  } else
    allocdiscnode(p, zeros, zeros2, endsite);
  (*p)->back = NULL;
  (*p)->next = NULL;
  (*p)->tip = false;
  (*p)->visited = false;
  (*p)->index = i;
  (*p)->numdesc = 0;
  (*p)->sumsteps = 0.0;
}  /* gnudisctreenode */


void setupnode(node *p, long i)
{ /* initialization of node pointers, variables */

  p->next = NULL;
  p->back = NULL;
  p->times_in_tree = (double) i * 1.0;
  p->index = i;
  p->tip = false;
}  /* setupnode */


node *pnode(tree *t, node *p) {
  /* Get the "parent nodelet" of p's node group */
  return t->nodep[p->index - 1];
}


long count_sibs (node *p)
{ /* Count the number of nodes in a ring, return the total number of */
  /* nodes excluding the one passed into the function (siblings)     */
  node *q;
  long return_int = 0;
  
  if (p->tip) {
   printf ("Error: the function count_sibs called on a tip.  This is a bug.\n");
    exxit (-1);
  }

  q = p->next;
  while (q != p) {
    if (q == NULL) {
      printf ("Error: a loop of nodes was not closed.\n");
      exxit (-1);
    } else {
      return_int++;
      q = q->next;
    }
  }

  return return_int;
}  /* count_sibs */


void inittrav (node *p)
{ /* traverse to set pointers uninitialized on inserting */
  long i, num_sibs;
  node *sib_ptr;
  
  if (p == NULL)
    return;
  if (p->tip)
    return;
  num_sibs = count_sibs (p);
  sib_ptr  = p;
  for (i=0; i < num_sibs; i++) {
    sib_ptr              = sib_ptr->next;
    sib_ptr->initialized = false;
    inittrav(sib_ptr->back);
  }
} /* inittrav */


void commentskipper(FILE ***intree, long *bracket)
{ /* skip over comment bracket contents in reading tree */
  char c;
  
  c = gettc(**intree);
  
  while (c != ']') {
    
    if(feof(**intree)) {
      printf("\n\nERROR: Unmatched comment brackets\n\n");
      exxit(-1);
    }

    if(c == '[') {
      (*bracket)++;
      commentskipper(intree, bracket);
    }
    c = gettc(**intree);
  }
  (*bracket)--;
}  /* commentskipper */


long countcomma(FILE **treefile, long *comma)
{
  /* Modified by Dan F. 11/10/96 */ 

  /* countcomma rewritten so it passes back both lparen+comma to allocate nodep
    and a pointer to the comma variable.  This allows the tree to know how many
    species exist, and the tips to be placed in the front of the nodep array */
  /* The next line inserted so this function leaves the file pointing
     to where it found it, not just re-winding it. */
  /* long orig_position = ftell (*treefile); */

  fpos_t orig_position;
  Char c;
  long  lparen = 0;
  long bracket = 0;

  /* Save the file position */
  if ( fgetpos(*treefile, &orig_position) != 0 ) {
    printf("\n\nERROR: Could not save file position!\n\n");
    exxit(-1);
  }

  (*comma) = 0;

  for (;;){
    c = getc(*treefile);
    if (feof(*treefile))
      break;
    if (c == ';')
      break;
    if (c == ',')
      (*comma)++;
    if (c == '(')
         lparen++;
    if (c == '[') {
      bracket++;
      commentskipper(&treefile, &bracket);
    }
  }

  /* Don't just rewind, */
  /* rewind (*treefile); */
  /* Re-set to where it pointed when the function was called */

  /* fseek (*treefile, orig_position, SEEK_SET); */

  fsetpos(*treefile, &orig_position);

  return lparen + (*comma);
}  /*countcomma*/


long countsemic(FILE **treefile)
{ /* Used to determine the number of user trees.  Return
     either a: the number of semicolons in the file outside comments
     or b: the first integer in the file */
  Char c;
  long return_val, semic = 0;
  long bracket = 0;
  
  /* Eat all whitespace */
  c = gettc(*treefile);
  while ((c == ' ')  ||
      (c == '\t') ||
      (c == '\n')) {
    c = gettc(*treefile);
  }

  /* Then figure out if the first non-white character is a digit; if
     so, return it */
  if (isdigit (c)) {
    ungetc(c, *treefile);
    if (fscanf((*treefile), "%ld", &return_val) != 1) {
      printf("Error reading number of trees in tree file.\n\n");
      exxit(-1);
    }
  } else {

    /* Loop past all characters, count the number of semicolons
       outside of comments */
    for (;;){
      c = fgetc(*treefile);
      if (feof(*treefile))
     break;
      if (c == ';')
     semic++;
      if (c == '[') {
     bracket++;
     commentskipper(&treefile, &bracket);
      }
    }
    return_val = semic;
  }

  rewind (*treefile);
  return return_val;
}  /* countsemic */


void hookup(node *p, node *q)
{ /* hook together two nodes */
  assert(p != NULL);
  assert(q != NULL);
  p->back = q;
  q->back = p;
}  /* hookup */


void unhookup(node *p, node *q)
{
  /* unhook two nodes. Not strictly required, but helps check assumptions */
  assert(p != NULL);
  assert(q != NULL);
  assert(p->back != NULL);
  assert(q->back != NULL);
  assert(p->back == q);
  assert(q->back == p);
  p->back = NULL;
  q->back = NULL;
}


void link_trees(long local_nextnum, long nodenum, long local_nodenum,
        pointarray nodep)
{
  if(local_nextnum == 0)
    hookup(nodep[nodenum], nodep[local_nodenum]);
  else if(local_nextnum == 1)
      hookup(nodep[nodenum], nodep[local_nodenum]->next);
    else if(local_nextnum == 2)
        hookup(nodep[nodenum], nodep[local_nodenum]->next->next);
      else
        printf("Error in Link_Trees()");
} /* link_trees() */


void allocate_nodep(pointarray *nodep, FILE **treefile, long  *precalc_tips)  
{ /* pre-compute space and allocate memory for nodep */

  long numnodes;      /* returns number commas & (    */
  long numcom = 0;        /* returns number commas */
  
  numnodes = countcomma(treefile, &numcom) + 1;
  *nodep      = (pointarray)Malloc(2*numnodes*sizeof(node *));

  (*precalc_tips) = numcom + 1;        /* this will be used in placing the
                                          tip nodes in the front region of
                                          nodep.  Used for species check?  */
} /* allocate_nodep -plc */


void malloc_pheno (node *p, long endsite, long rcategs)
{ /* Allocate the phenotype arrays; used by dnaml */
  long i;

  p->x  = (phenotype)Malloc(endsite*sizeof(ratelike));
  p->underflows = (double *)Malloc(endsite * sizeof(double));
  for (i = 0; i < endsite; i++)
    p->x[i]  = (ratelike)Malloc(rcategs*sizeof(sitelike));
} /* malloc_pheno */


void malloc_ppheno (node *p,long endsite, long rcategs)
{
  /* Allocate the phenotype arrays; used by proml */
  long i;

  p->protx  = (pphenotype)Malloc(endsite*sizeof(pratelike));
  p->underflows  = (double *)Malloc(endsite*sizeof(double));
  
  for (i = 0; i < endsite; i++)
    p->protx[i]  = (pratelike)Malloc(rcategs*sizeof(psitelike));
} /* malloc_ppheno */


long take_name_from_tree (Char *ch, Char *str, FILE *treefile)
{
  /* This loop reads a name from treefile and stores it in *str.
     Returns the length of the name string. str must be at
     least MAXNCH bytes, but no effort is made to null-terminate
     the string. Underscores and newlines are converted to spaces.
     Characters beyond MAXNCH are discarded. */

  long name_length = 0;

  do {
    if ((*ch) == '_')
      (*ch) = ' ';
    if ( name_length < MAXNCH )
      str[name_length++] = (*ch);
    if (eoln(treefile)) 
      scan_eoln(treefile);
    (*ch) = gettc(treefile);
    if (*ch == '\n')
      *ch = ' ';
  } while ( strchr(":,)[;", *ch) == NULL );

  return name_length;
}  /* take_name_from_tree */


void match_names_to_data (Char *str, pointarray treenode, node **p, long spp)
{
  /* This loop matches names taken from treefile to indexed names in
     the data file */

  boolean found;
  long i, n;

  n = 1;  
  do {
    found = true;
    for (i = 0; i < nmlngth; i++) {
      found = (found && ((str[i] == nayme[n - 1][i]) ||
        (((nayme[n - 1][i] == '_') && (str[i] == ' ')) ||
        ((nayme[n - 1][i] == ' ') && (str[i] == '\0')))));
    }
    
    if (found)
      *p = treenode[n - 1];
    else
      n++;

  } while (!(n > spp || found));
  
  if (n > spp) {
    printf("\n\nERROR: Cannot find species: ");
    for (i = 0; (str[i] != '\0') && (i < MAXNCH); i++)
      putchar(str[i]);
    printf(" in data file\n\n");
    exxit(-1);
  }
}  /* match_names_to_data */


void addelement(node **p, node *q, Char *ch, long *parens, FILE *treefile,
        pointarray treenode, boolean *goteof, boolean *first, pointarray nodep,
        long *nextnode, long *ntips, boolean *haslengths, node **grbg,
        initptr initnode, boolean unifok, long maxnodes)
{
  /* Recursive procedure adds nodes to user-defined tree
     This is the main (new) tree-reading procedure */
  
  node *pfirst;
  long i, len = 0, nodei = 0;
  boolean notlast;
  Char str[MAXNCH+1];
  node *r;
  long furs = 0;


  if ((*ch) == '(') {
    (*nextnode)++;          /* get ready to use new interior node */
    nodei = *nextnode;      /* do what needs to be done at bottom */
    if ( maxnodes != -1 && nodei > maxnodes) {
      printf("ERROR in input tree file: Attempting to allocate too\n"); 
      printf("many nodes. This is usually caused by a unifurcation.\n");  
      printf("To use this tree with this program  use Retree to read\n");
      printf("and write this tree.\n");
      exxit(-1);
    }
    
    /* do what needs to be done at bottom */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, bottom, treenode, nodep, str, ch, treefile);
    pfirst      = (*p);
    notlast = true;
    while (notlast) {          /* loop through immediate descendants */
      furs++;
      (*initnode)(&(*p)->next, grbg, q,
                   len, nodei, ntips, parens, nonbottom, treenode,
                   nodep, str, ch, treefile);
                       /* ... doing what is done before each */
      r = (*p)->next;
      getch(ch, parens, treefile);      /* look for next character */
      
       /* handle blank names */
      if((*ch) == ',' || (*ch) == ':'){
        ungetc((*ch), treefile);
        *ch = 0;
      } else if((*ch)==')'){
        ungetc((*ch), treefile);
        (*parens)++;
        *ch = 0;
      }
 
      addelement(&(*p)->next->back, (*p)->next, ch, parens, treefile,
        treenode, goteof, first, nodep, nextnode, ntips,
        haslengths, grbg, initnode, unifok, maxnodes);  

      (*initnode)(&r, grbg, q, len, nodei, ntips,
                    parens, hslength, treenode, nodep, str, ch, treefile);
                           /* do what is done after each about length */
      pfirst->numdesc++;               /* increment number of descendants */
      *p = r;                         /* make r point back to p */

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }
    if ( furs <= 1 && !unifok ) {
      printf("ERROR in input tree file: A Unifurcation was detetected.\n");
      printf("To use this tree with this program use retree to read and");
      printf(" write this tree\n");
      exxit(-1);
    }
    
    (*p)->next = pfirst;
    (*p)       = pfirst;

  } else if ((*ch) != ')') {       /* if it's a species name */
    for (i = 0; i < MAXNCH+1; i++)   /* fill string with nulls */
      str[i] = '\0';

    len = take_name_from_tree (ch, str, treefile);  /* get the name */

    if ((*ch) == ')')
      (*parens)--;         /* decrement count of open parentheses */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, tip, treenode, nodep, str, ch, treefile);
                              /* do what needs to be done at a tip */
  } else
    getch(ch, parens, treefile);
  if (q != NULL)
    hookup(q, (*p));                    /* now hook up */
  (*initnode)(p, grbg, q, len, nodei, ntips, 
                parens, iter, treenode, nodep, str, ch, treefile);
                              /* do what needs to be done to variable iter */
  if ((*ch) == ':')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, length, treenode, nodep, str, ch, treefile);
                                   /* do what needs to be done with length */
  else if ((*ch) != ';' && (*ch) != '[')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, hsnolength, treenode, nodep, str, ch, treefile);
                             /* ... or what needs to be done when no length */
  if ((*ch) == '[')
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, treewt, treenode, nodep, str, ch, treefile);
                             /* ... for processing a tree weight */
  else if ((*ch) == ';')     /* ... and at end of tree */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, unittrwt, treenode, nodep, str, ch, treefile);
}  /* addelement */


void treeread (FILE *treefile, node **root, pointarray treenode,
        boolean *goteof, boolean *first, pointarray nodep, 
        long *nextnode, boolean *haslengths, node **grbg, initptr initnode,
        boolean unifok, long maxnodes)
{
  /* read in user-defined tree and set it up */
  /* Eats blank lines and everything up to the first open paren, then
   * calls the recursive function addelement, which builds the
   * tree and calls back to initnode. */
  char  ch;
  long parens = 0;
  long ntips = 0;

  (*goteof) = false;
  (*nextnode) = spp;

  /* eat blank lines */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }
  if (haslengths != NULL)
    *haslengths = true; 
  addelement(root, NULL, &ch, &parens, treefile,
         treenode, goteof, first, nodep, nextnode, &ntips,
         haslengths, grbg, initnode, unifok, maxnodes);

  /* Eat blank lines and end of current line*/
  do {
    scan_eoln(treefile);
  }
  while (eoln(treefile) && !eoff(treefile));

  if (first)
    *first = false;
  if (parens != 0) {
    printf("\n\nERROR in tree file: unmatched parentheses\n\n");
    exxit(-1);
  }
}  /* treeread */


void addelement2(node *q, Char *ch, long *parens, FILE *treefile,
        pointarray treenode, boolean lngths, double *trweight, boolean *goteof,
        long *nextnode, long *ntips, long no_species, boolean *haslengths,
        boolean unifok, long maxnodes)
{
  /* recursive procedure adds nodes to user-defined tree
     -- old-style bifurcating-only version */

  node *pfirst = NULL, *p;
  long i, len, current_loop_index;
  boolean notlast, minusread;
  Char str[MAXNCH];
  double valyew, divisor;
  long furs = 0; 

  if ((*ch) == '(') {

    current_loop_index = (*nextnode) + spp;
    (*nextnode)++;

    if ( maxnodes != -1 && current_loop_index > maxnodes) {
      printf("ERROR in intree file: Attempting to allocate too many nodes\n");
      printf("This is usually caused by a unifurcation.  To use this\n");
      printf("intree with this program  use retree to read and write\n"); 
      printf("this tree.\n");
      exxit(-1);
    }
    /* This is an assignment of an interior node */
    p = treenode[current_loop_index];
    pfirst = p;
    notlast = true;
    while (notlast) {
      furs++;
      /* This while loop goes through a circle (triad for
      bifurcations) of nodes */
      p = p->next;
      /* added to ensure that non base nodes in loops have indices */
      p->index = current_loop_index + 1;
      
      getch(ch, parens, treefile);
      
      addelement2(p, ch, parens, treefile, treenode, lngths, trweight,
        goteof, nextnode, ntips, no_species, haslengths, unifok, maxnodes);

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }
    if ( furs <= 1 && !unifok ) {
      printf("ERROR in intree file: A Unifurcation was detected.\n");
      printf("To use this intree with this program use retree to read and");
      printf(" write this tree\n");
      exxit(-1);
    }

  } else if ((*ch) != ')') {
    for (i = 0; i < MAXNCH; i++) 
      str[i] = '\0';
    len = take_name_from_tree (ch, str, treefile);
    match_names_to_data (str, treenode, &p, spp);
    pfirst = p;
    if ((*ch) == ')')
      (*parens)--;
    (*ntips)++;
    strncpy (p->nayme, str, len);
  } else
    getch(ch, parens, treefile);
  
  if ((*ch) == '[') {    /* getting tree weight from last comment field */
    if (!eoln(treefile)) {
      if (fscanf(treefile, "%lf", trweight) == 1) {
        getch(ch, parens, treefile);
        if (*ch != ']') {
          printf("\n\nERROR: Missing right square bracket\n\n");
          exxit(-1);
        }
        else {
          getch(ch, parens, treefile);
          if (*ch != ';') {
            printf("\n\nERROR: Missing semicolon after square brackets\n\n");
            exxit(-1);
          }
        }
      }
      else {
        printf("\n\nERROR: Expecting tree weight in last comment field.\n\n");
        exxit(-1);
      }
    }
  }
  else if ((*ch) == ';') {
    (*trweight) = 1.0 ;
    if (!eoln(treefile))
      printf("WARNING: tree weight set to 1.0\n");
  }
  else if (haslengths != NULL)
    (*haslengths) = ((*haslengths) && q == NULL);
  
  if (q != NULL)
    hookup(q, pfirst);
  
  if ((*ch) == ':') {
    processlength(&valyew, &divisor, ch,
       &minusread, treefile, parens);
    if (q != NULL) {
      if (!minusread)
        q->oldlen = valyew / divisor;
      else
        q->oldlen = 0.0;
      if (lngths) {
        q->v = valyew / divisor;
        q->back->v = q->v;
        q->iter = false;
        q->back->iter = false;
      }
    }
  }
  
}  /* addelement2 */


void treeread2 (FILE *treefile, node **root, pointarray treenode,
        boolean lngths, double *trweight, boolean *goteof,
        boolean *haslengths, long *no_species, boolean unifok, long maxnodes)
{
  /* read in user-defined tree and set it up
     -- old-style bifurcating-only version */
  char  ch;
  long parens = 0;
  long ntips = 0;
  long nextnode;
 
  (*goteof) = false;
  nextnode = 0;

  /* Eats all blank lines at start of file */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }

  addelement2(NULL, &ch, &parens, treefile, treenode, lngths, trweight,
      goteof, &nextnode, &ntips, (*no_species), haslengths, unifok, maxnodes);
  (*root) = treenode[*no_species];

  /*eat blank lines */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  (*root)->oldlen = 0.0;

  if (parens != 0) {
    printf("\n\nERROR in tree file:  unmatched parentheses\n\n");
    exxit(-1);
  }
}  /* treeread2 */


void exxit(int exitcode)
{ /* Terminate the program with exit code exitcode.
   * On Windows, supplying a nonzero exitcode will print a message and wait
   * for the user to hit enter. */

  exit (exitcode);
} /* exxit */


char gettc(FILE* file) 
{ /* Return the next character in file.
   * If EOF is reached, print an error and die.
   * DOS ('\r\n') and Mac ('\r') newlines are returned as a single '\n'. */
  int ch;

  ch=getc(file);

  if ( ch == EOF )
    EOF_error();

  if ( ch == '\r' ) {
    ch = getc(file);
    if ( ch != '\n' )
      ungetc(ch, file);
    ch = '\n';
  }
  return ch;
} /* gettc */

void unroot(tree *t, long nonodes) 
{
  /* used by fitch, restml and contml */
  if (t->start->back == NULL) { 
    if (t->start->next->back->tip)
      t->start = t->start->next->next->back;
    else  t->start = t->start->next->back;
  }
  if (t->start->next->back == NULL) {
    if (t->start->back->tip)
      t->start = t->start->next->next->back;
    else t->start = t->start->back;
  }
  if (t->start->next->next->back == NULL)  {
    if (t->start->back->tip)
      t->start = t->start->next->back;
    else t->start = t->start->back;
  }
    

  unroot_r(t->start,t->nodep,nonodes);
  unroot_r(t->start->back, t->nodep, nonodes);
}


void unroot_here(node* root, node** nodep, long nonodes)
{
  node* tmpnode;
  double newl;
  /* used by unroot */
  /* assumes bifurcation this is ok in the programs that use it */
 
  newl = root->next->oldlen + root->next->next->oldlen;
  root->next->back->oldlen = newl;
  root->next->next->back->oldlen = newl;

  newl = root->next->v + root->next->next->v;
  root->next->back->v = newl;
  root->next->next->back->v = newl;

  root->next->back->back=root->next->next->back;
  root->next->next->back->back = root->next->back;
  
  while ( root->index != nonodes ) {
    tmpnode = nodep[ root->index ];
    nodep[root->index] = root;
    root->index++;
    root->next->index++;
    root->next->next->index++;
    nodep[root->index - 2] = tmpnode;
    tmpnode->index--;
    tmpnode->next->index--;
    tmpnode->next->next->index--;
  }
}


void unroot_r(node* p, node** nodep, long nonodes) 
{
  /* used by unroot */
  node *q;
  
  if ( p->tip) return;

  q = p->next;
  while ( q != p ) {
    if (q->back == NULL)
      unroot_here(q, nodep, nonodes);
    else unroot_r(q->back, nodep, nonodes);
    q = q->next;
  }
}

void clear_connections(tree *t, long nonodes) 
{
  long i;
  node *p;
  for ( i = 0 ; i < nonodes ; i++) {
    p = t->nodep[i];
    if (p != NULL) {
      p->back = NULL;
      p->v = 0;
      for (p = p->next; p && p != t->nodep[i]; p = p->next) {
        p->next->back = NULL;
        p->next->v    = 0;
      }
    }
  }
}

/* These functions are temporarily used for translating the fixed-width
 * space-padded nayme array to an array of null-terminated char *. */
char **stringnames_new(void)
{
  /* Copy nayme array to null terminated strings and return array of char *.
   * Spaces are stripped from end of naym's.
   * Returned array size is spp+1; last element is NULL. */

  char **names;
  char *ch;
  long i;

  names = (char **)Malloc((spp+1) * sizeof(char *));

  for ( i = 0; i < spp; i++ ) {
    names[i] = (char *)Malloc((MAXNCH+1) * sizeof(char));
    strncpy(names[i], nayme[i], MAXNCH);
    names[i][MAXNCH] = '\0';
    /* Strip trailing spaces */
    for ( ch = names[i] + MAXNCH - 1; *ch == ' ' || *ch == '\0'; ch-- )
      *ch = '\0';
  }
  names[spp] = NULL;

  return names;
}

void stringnames_delete(char **names)
{
  /* Free a string array returned by stringnames_new() */
  long i;

  assert( names != NULL );
  for ( i = 0; i < spp; i++ ) {
    assert( names[i] != NULL );
    free(names[i]);
  }
  free(names);
}

int fieldwidth_double(double val, unsigned int precision)
{
  /* Printf a double to a temporary buffer with specified precision using %g
   * and return its length. Precision must not be greater than 999,999 */

  char format[10];
  char buf[0x200]; /* TODO: What's the largest possible? */

  if (precision > 999999)
    abort();

  sprintf(format, "%%.%uf", precision); /* %.Nf */
  /* snprintf() would be better, but is it avaliable on all systems? */
  return sprintf(buf, format, val);
}

void output_matrix_d(FILE *fp, double **matrix,
    unsigned long rows, unsigned long cols,
    char **row_head, char **col_head, int flags)
{
  /*
   * Print a matrix of double to file. Headings are given in row_head and
   * col_head, either of which may be NULL to indicate that headings should not
   * be printed. Otherwise, they must be null-terminated arrays of pointers to
   * null-terminalted character arrays.
   *
   * The macro OUTPUT_PRECISION defines the number of significant figures to
   * print, and OUTPUT_TEXTWIDTH defines the maximum length of each line.
   * 
   * Optional formatting is specified by flags argument, using macros MAT_*
   * defined in phylip.h.
   */

  unsigned     *colwidth;               /* [0..spp-1] min width of each column */
  unsigned      headwidth;              /* min width of row header column */
  unsigned long linelen;                /* length of current printed line */
  unsigned      fw;
  unsigned long row, col;
  unsigned long i;
  unsigned long cstart, cend;
  unsigned long textwidth = OUTPUT_TEXTWIDTH;
  const unsigned int     gutter = 1;
  boolean       do_block;
  boolean       lower_triangle;
  boolean       border;
  boolean       output_cols;
  boolean       pad_row_head;

  if ( flags & MAT_NOHEAD )
    col_head = NULL;
  if ( flags & MAT_NOBREAK )
    textwidth = 0;
  do_block = (flags & MAT_BLOCK) && (textwidth > 0);
  lower_triangle = flags & MAT_LOWER;
  border = flags & MAT_BORDER;
  output_cols = flags & MAT_PCOLS;
  pad_row_head = flags & MAT_PADHEAD;

  /* Determine minimal width for row headers, if given */
  headwidth = 0;
  if ( row_head != NULL ) { 
    for (row = 0; row < rows; row++) {
      fw = strlen(row_head[row]);
      if ( headwidth < fw )
        headwidth = fw;
    }
  }

  /* Enforce minimum of 10 ch for machine-readable output */
  if ( (pad_row_head) && (headwidth < 10) )
    headwidth = 10;
  
  /* Determine minimal width for each matrix col */
  colwidth = (unsigned int *)Malloc(spp * sizeof(int));
  for (col = 0; col < cols; col++) {
    if ( col_head != NULL )
      colwidth[col] = strlen(col_head[col]);
    else
      colwidth[col] = 0;
    for (row = 0; row < rows; row++) {
      fw = fieldwidth_double(matrix[row][col], PRECISION);
      if ( colwidth[col] < fw )
        colwidth[col] = fw;
    }
  }
  
  /*** Print the matrix ***/
  /* Number of columns if requested */
  if ( output_cols ) {
    fprintf(fp, "%5lu\n", cols);
  }
  
  /* Omit last column for lower triangle */
  if ( lower_triangle )
    cols--;

  /* Blocks */
  cstart = cend = 0;
  while ( cend != cols ) {
    if ( do_block ) {
      linelen = headwidth;
      for ( col = cstart; col < cols; col++ ) {
        if ( linelen + colwidth[col] + gutter > textwidth ) {
          break;
        }
        linelen += colwidth[col] + gutter;
      }
      cend = col;
      /* Always print at least one, regardless of line len */
      if ( cend == cstart )
        cend++;
    } else {
      cend = cols;
    }


    /* Column headers */
    if ( col_head != NULL ) {
      /* corner space */
      for ( i = 0; i < headwidth; i++ )
        putc(' ', fp);
      if ( border ) {
        for ( i = 0; i < gutter+1; i++ )
          putc(' ', fp);
      }
      /* Names */
      for ( col = cstart; col < cend; col++ ) {
        for ( i = 0; i < gutter; i++ )
          putc(' ', fp);
        /* right justify */
        fw = strlen(col_head[col]);
        for ( i = 0; i < colwidth[col] - fw; i++ )
          putc(' ', fp);
        fputs(col_head[col], fp);
      }
      putc('\n', fp);
    }
    
    /* Top border */
    if ( border ) {
      for ( i = 0; i < headwidth + gutter; i++ )
        putc(' ', fp);
      putc('\\', fp);
      for ( col = cstart; col < cend; col++ ) {
        for ( i = 0; i < colwidth[col] + gutter; i++ )
          putc('-', fp);
      }
      putc('\n', fp);
    }
    
    /* Rows */
    for (row = 0; row < rows; row++) {
      /* Row header, if given */
      if ( row_head != NULL ) {
        /* right-justify for non-machine-readable */
        if ( !pad_row_head ) {
          for ( i = strlen(row_head[row]); i < headwidth; i++ )
            putc(' ', fp);
        }
        fputs(row_head[row], fp);
        /* left-justify for machine-readable */
        if ( pad_row_head ) {
          for ( i = strlen(row_head[row]); i < headwidth; i++ )
            putc(' ', fp);
        }
      }
      linelen = headwidth;

      /* Left border */
      if ( border ) {
        for ( i = 0; i < gutter; i++ )
          putc(' ', fp);
        putc('|', fp);
        linelen += 2;
      }
      
      /* Row data */
      for (col = cstart; col < cend; col++) { /* cols */
        /* Stop after col == row for lower triangle */
        if ( lower_triangle && col >= row )
          break;
        /* Break line if going over max text width */
        if ( !do_block && textwidth > 0 ) {
          if ( linelen + colwidth[col] > textwidth ) {
            putc('\n', fp);
            linelen = 0;
          }
          linelen += colwidth[col] + gutter;
        }
        
        for ( i = 0; i < gutter; i++ )
          putc(' ', fp);

        /* Print the datum */
        fprintf(fp, "%*.6f", colwidth[col], matrix[row][col]);
      }
      putc('\n', fp);
    } /* End of row */
    if (col_head != NULL)
      putc('\n', fp); /* blank line */
    cstart = cend;
  } /* End of block */
  free(colwidth);
} /* output_matrix_d */




////////////////////////  functions from seq.c ///////////////////////




void fix_x(node* p,long site, double maxx, long rcategs)
{ /* dnaml dnamlk */
  long i,j;
  p->underflows[site] += log(maxx);

  for ( i = 0 ; i < rcategs ; i++ ) {
    for ( j = 0 ; j < ((long)T - (long)A + 1) ; j++)
      p->x[site][i][j] /= maxx;
  }
} /* fix_x */


void fix_protx(node* p,long site, double maxx, long rcategs) 
{ /* proml promlk */
  long i,m;

  p->underflows[site] += log(maxx);

  for ( i = 0 ; i < rcategs  ; i++ ) 
    for (m = 0; m <= 19; m++)
      p->protx[site][i][m] /= maxx;
} /* fix_protx */


void alloctemp(node **temp, long *zeros, long endsite)
{
  /*used in dnacomp and dnapenny */
  *temp = (node *)Malloc(sizeof(node));
  (*temp)->numsteps = (steptr)Malloc(endsite*sizeof(long));
  (*temp)->base = (baseptr)Malloc(endsite*sizeof(long));
  (*temp)->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  memcpy((*temp)->base, zeros, endsite*sizeof(long));
  memcpy((*temp)->numsteps, zeros, endsite*sizeof(long));
  zeronumnuc(*temp, endsite);
}  /* alloctemp */


void freetemp(node **temp)
{
  /* used in dnacomp, dnapars, & dnapenny */
  free((*temp)->numsteps);
  free((*temp)->base);
  free((*temp)->numnuc);
  free(*temp);
}  /* freetemp */


void freetree2 (pointarray treenode, long nonodes)
{
  /* The natural complement to alloctree2.  Free all elements of all
  the rings (normally triads) in treenode */
  long i;
  node *p, *q;

  /* The first spp elements are just nodes, not rings */
  for (i = 0; i < spp; i++)
    free (treenode[i]);

  /* The rest are rings */
  for (i = spp; i < nonodes; i++) {
    p = treenode[i]->next;
    while (p != treenode[i]) {
      q = p->next;
      free (p);
      p = q;
    }
    /* p should now point to treenode[i], which has yet to be freed */
    free (p);
  }
  free (treenode);
}  /* freetree2 */


void inputdata(long chars)
{
  /* input the names and sequences for each species */
  /* used by dnacomp, dnadist, dnainvar, dnaml, dnamlk, dnapars, & dnapenny */
  long i, j, k, l, basesread, basesnew=0;
  Char charstate;
  boolean allread, done;

  if (printdata)
    headings(chars, "Sequences", "---------");
  basesread = 0;
  allread = false;
  while (!(allread)) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp) {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i-1);
      j = (interleaved) ? basesread : 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          if ((strchr("ABCDGHKMNRSTUVWXY?O-",charstate)) == NULL){
            printf("ERROR: bad base: %c at site %5ld of species %3ld\n",
                   charstate, j+1, i);
            if (charstate == '.') {
              printf("       Periods (.) may not be used as gap characters.\n");
              printf("       The correct gap character is (-)\n");
            }
            exxit(-1);
          }
          j++;
          y[i - 1][j - 1] = charstate;
        }
        if (interleaved)
          continue;
        if (j < chars) 
          scan_eoln(infile);
        else if (j == chars)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j;

      scan_eoln(infile);
    
      if ((interleaved && j != basesnew) ||
          (!interleaved && j != chars)) {
        printf("\nERROR: sequences out of alignment at position %ld", j+1);
        printf(" of species %ld\n\n", i);
        exxit(-1);
      }
      i++;
    }
    
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == chars);
    } else
      allread = (i > spp);
  }
  if (!printdata)
    return;
  for (i = 1; i <= ((chars - 1) / 60 + 1); i++) {
    for (j = 1; j <= spp; j++) {
      for (k = 0; k < nmlngth; k++)
        putc(nayme[j - 1][k], outfile);
      fprintf(outfile, "   ");
      l = i * 60;
      if (l > chars)
        l = chars;
      for (k = (i - 1) * 60 + 1; k <= l; k++) {
        if (dotdiff && (j > 1 && y[j - 1][k - 1] == y[0][k - 1]))
          charstate = '.';
        else
          charstate = y[j - 1][k - 1];
        putc(charstate, outfile);
        if (k % 10 == 0 && k % 60 != 0)
          putc(' ', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* inputdata */


void alloctree(pointarray *treenode, long nonodes, boolean usertree)
{
  /* allocate treenode dynamically */
  /* used in dnapars, dnacomp, dnapenny & dnamove */
  long i, j;
  node *p, *q;

  *treenode = (pointarray)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < spp; i++) {
    (*treenode)[i] = (node *)Malloc(sizeof(node));
    (*treenode)[i]->tip = true;
    (*treenode)[i]->index = i+1;
    (*treenode)[i]->iter = true;
    (*treenode)[i]->branchnum = 0;
    (*treenode)[i]->initialized = true;
  }
  if (!usertree)
    for (i = spp; i < nonodes; i++) {
      q = NULL;
      for (j = 1; j <= 3; j++) {
        p = (node *)Malloc(sizeof(node));
        p->tip = false;
        p->index = i+1;
        p->iter = true;
        p->branchnum = 0;
        p->initialized = false;
        p->next = q;
        q = p;
      }
      p->next->next->next = p;
      (*treenode)[i] = p;
    }
} /* alloctree */


void allocx(long nonodes, long rcategs, pointarray treenode, boolean usertree)
{
  /* allocate x dynamically */
  /* used in dnaml & dnamlk */
  long i, j, k;
  node *p;

  for (i = 0; i < spp; i++){
    treenode[i]->x = (phenotype)Malloc(endsite*sizeof(ratelike));
    treenode[i]->underflows = (double *)Malloc(endsite * sizeof (double));
    for (j = 0; j < endsite; j++)
      treenode[i]->x[j]  = (ratelike)Malloc(rcategs*sizeof(sitelike));
  }
  if (!usertree) {
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
        p->underflows = (double *)Malloc (endsite * sizeof (double));
        p->x = (phenotype)Malloc(endsite*sizeof(ratelike));
        for (k = 0; k < endsite; k++)
          p->x[k] = (ratelike)Malloc(rcategs*sizeof(sitelike));
        p = p->next;
      }
    }
  }
}  /* allocx */


void prot_allocx(long nonodes, long rcategs, pointarray treenode, 
                        boolean usertree)
{
  /* allocate x dynamically */
  /* used in proml          */
  long i, j, k;
  node *p;

  for (i = 0; i < spp; i++){
    treenode[i]->protx = (pphenotype)Malloc(endsite*sizeof(pratelike));
    treenode[i]->underflows = (double *)Malloc(endsite*sizeof(double));
    for (j = 0; j < endsite; j++)
      treenode[i]->protx[j]  = (pratelike)Malloc(rcategs*sizeof(psitelike));
  }  
  if (!usertree) {
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
        p->protx = (pphenotype)Malloc(endsite*sizeof(pratelike));
        p->underflows = (double *)Malloc(endsite*sizeof(double));
        for (k = 0; k < endsite; k++)
          p->protx[k] = (pratelike)Malloc(rcategs*sizeof(psitelike));
        p = p->next;
      } 
    }  
  } 
}  /* prot_allocx */




void setuptree(pointarray treenode, long nonodes, boolean usertree)
{
  /* initialize treenodes */
  long i;
  node *p;

  for (i = 1; i <= nonodes; i++) {
    if (i <= spp || !usertree) {
      treenode[i-1]->back = NULL;
      treenode[i-1]->tip = (i <= spp);
      treenode[i-1]->index = i;
      treenode[i-1]->numdesc = 0;
      treenode[i-1]->iter = true;
      treenode[i-1]->initialized = true;
      treenode[i-1]->tyme =  0.0;
    }
  }
  if (!usertree) {
    for (i = spp + 1; i <= nonodes; i++) {
      p = treenode[i-1]->next;
      while (p != treenode[i-1]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        p->numdesc = 0;
        p->iter = true;
        p->initialized = false;
        p->tyme = 0.0;
        p = p->next;
      }
    }
  }
} /* setuptree */


void setuptree2(tree *a)
{
  /* initialize a tree */
  /* used in dnaml, dnamlk, & restml */

  a->likelihood = -999999.0;
  a->start = a->nodep[0]->back;
  a->root = NULL;
} /* setuptree2 */


void alloctip(node *p, long *zeros)
{ /* allocate a tip node */
  /* used by dnacomp, dnapars, & dnapenny */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->base = (baseptr)Malloc(endsite*sizeof(long));
  p->oldbase = (baseptr)Malloc(endsite*sizeof(long));
  memcpy(p->base, zeros, endsite*sizeof(long));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->oldbase, zeros, endsite*sizeof(long));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
}  /* alloctip */




void getbasefreqs(double freqa, double freqc, double freqg, double freqt,
            double *freqr, double *freqy, double *freqar, double *freqcy,
            double *freqgr, double *freqty, double *ttratio, double *xi,
            double *xv, double *fracchange, boolean freqsfrom,
            boolean printdata)
{
  /* used by dnadist, dnaml, & dnamlk */
  double aa, bb;

  if (printdata) {
    putc('\n', outfile);
    if (freqsfrom)
      fprintf(outfile, "Empirical ");
    fprintf(outfile, "Base Frequencies:\n\n");
    fprintf(outfile, "   A    %10.5f\n", freqa);
    fprintf(outfile, "   C    %10.5f\n", freqc);
    fprintf(outfile, "   G    %10.5f\n", freqg);
    fprintf(outfile, "  T(U)  %10.5f\n", freqt);
    fprintf(outfile, "\n");
  }
  *freqr = freqa + freqg;
  *freqy = freqc + freqt;
  *freqar = freqa / *freqr;
  *freqcy = freqc / *freqy;
  *freqgr = freqg / *freqr;
  *freqty = freqt / *freqy;
  aa = *ttratio * (*freqr) * (*freqy) - freqa * freqg - freqc * freqt;
  bb = freqa * (*freqgr) + freqc * (*freqty);
  *xi = aa / (aa + bb);
  *xv = 1.0 - *xi;
  if (*xi < 0.0) {
    printf("\n WARNING: This transition/transversion ratio\n");
    printf(" is impossible with these base frequencies!\n");
    *xi = 0.0;
    *xv = 1.0;
    (*ttratio) = (freqa*freqg+freqc*freqt)/((*freqr)*(*freqy));
    printf(" Transition/transversion parameter reset\n");
    printf("  so transition/transversion ratio is %10.6f\n\n", (*ttratio));
  }
  if (freqa <= 0.0)
    freqa = 0.000001;
  if (freqc <= 0.0)
    freqc = 0.000001;
  if (freqg <= 0.0)
    freqg = 0.000001;
  if (freqt <= 0.0)
    freqt = 0.000001;
  *fracchange = (*xi) * (2 * freqa * (*freqgr) + 2 * freqc * (*freqty)) +
      (*xv) * (1.0 - freqa * freqa - freqc * freqc - freqg * freqg
      - freqt * freqt);
}  /* getbasefreqs */


void empiricalfreqs(double *freqa, double *freqc, double *freqg,
                        double *freqt, steptr weight, pointarray treenode)
{
  /* Get empirical base frequencies from the data */
  /* used in dnaml & dnamlk */
  long i, j, k;
  double sum, suma, sumc, sumg, sumt, w;

  *freqa = 0.25;
  *freqc = 0.25;
  *freqg = 0.25;
  *freqt = 0.25;
  for (k = 1; k <= 8; k++) {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;
    for (i = 0; i < spp; i++) {
      for (j = 0; j < endsite; j++) {
        w = weight[j];
        sum = (*freqa) * treenode[i]->x[j][0][0];
        sum += (*freqc) * treenode[i]->x[j][0][(long)C - (long)A];
        sum += (*freqg) * treenode[i]->x[j][0][(long)G - (long)A];
        sum += (*freqt) * treenode[i]->x[j][0][(long)T - (long)A];
        suma += w * (*freqa) * treenode[i]->x[j][0][0] / sum;
        sumc += w * (*freqc) * treenode[i]->x[j][0][(long)C - (long)A] / sum;
        sumg += w * (*freqg) * treenode[i]->x[j][0][(long)G - (long)A] / sum;
        sumt += w * (*freqt) * treenode[i]->x[j][0][(long)T - (long)A] / sum;
      }
    }
    sum = suma + sumc + sumg + sumt;
    *freqa = suma / sum;
    *freqc = sumc / sum;
    *freqg = sumg / sum;
    *freqt = sumt / sum;
  }
  if (*freqa <= 0.0)
    *freqa = 0.000001;
  if (*freqc <= 0.0)
    *freqc = 0.000001;
  if (*freqg <= 0.0)
    *freqg = 0.000001;
  if (*freqt <= 0.0)
    *freqt = 0.000001;
}  /* empiricalfreqs */


void sitesort(long chars, steptr weight)
{
  /* Shell sort keeping sites, weights in same order */
  /* used in dnainvar, dnapars, dnacomp & dnapenny */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = chars / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= chars; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        tied = true;
        k = 1;
        while (k <= spp && tied) {
          flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
          tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j - 1];
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = weight[j - 1];
        weight[j - 1] = weight[j + gap - 1];
        weight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort */


void sitecombine(long chars)
{
  /* combine sites that have identical patterns */
  /* used in dnapars, dnapenny, & dnacomp */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < chars) {
    j = i + 1;
    tied = true;
    while (j <= chars && tied) {
      k = 1;
      while (k <= spp && tied) {
        tied = (tied &&
            y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (tied) {
        weight[i - 1] += weight[j - 1];
        weight[j - 1] = 0;
        ally[alias[j - 1] - 1] = alias[i - 1];
      }
      j++;
    }
    i = j - 1;
  }
}  /* sitecombine */


void sitescrunch(long chars)
{
  /* move so one representative of each pattern of
     sites comes first */
  /* used in dnapars & dnacomp */
  long i, j, itemp;
  boolean done, found;

  done = false;
  i = 1;
  j = 2;
  while (!done) {
    if (ally[alias[i - 1] - 1] != alias[i - 1]) {
      if (j <= i)
        j = i + 1;
      if (j <= chars) {
        do {
          found = (ally[alias[j - 1] - 1] == alias[j - 1]);
          j++;
        } while (!(found || j > chars));
        if (found) {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = weight[i - 1];
          weight[i - 1] = weight[j - 1];
          weight[j - 1] = itemp;
        } else
          done = true;
      } else
        done = true;
    }
    i++;
    done = (done || i >= chars);
  }
}  /* sitescrunch */


void sitesort2(long sites, steptr aliasweight)
{
  /* Shell sort keeping sites, weights in same order */
  /* used in dnaml & dnamnlk */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied, samewt;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        samewt = ((weight[jj - 1] != 0) && (weight[jg - 1] != 0))
                 || ((weight[jj - 1] == 0) && (weight[jg - 1] == 0));
        tied = samewt && (category[jj - 1] == category[jg - 1]);
        flip = ((!samewt) && (weight[jj - 1] == 0))
               || (samewt && (category[jj - 1] > category[jg - 1]));
        k = 1;
        while (k <= spp && tied) {
          flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
          tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j - 1];
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = aliasweight[j - 1];
        aliasweight[j - 1] = aliasweight[j + gap - 1];
        aliasweight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort2 */


void sitecombine2(long sites, steptr aliasweight)
{
  /* combine sites that have identical patterns */
  /* used in dnaml & dnamlk */
  long i, j, k;
  boolean tied, samewt;

  i = 1;
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      samewt = ((aliasweight[i - 1] != 0) && (aliasweight[j - 1] != 0))
               || ((aliasweight[i - 1] == 0) && (aliasweight[j - 1] == 0));
      tied = samewt
             && (category[alias[i - 1] - 1] == category[alias[j - 1] - 1]);
      k = 1;
      while (k <= spp && tied) {
        tied = (tied &&
            y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (!tied)
        break;
      aliasweight[i - 1] += aliasweight[j - 1];
      aliasweight[j - 1] = 0;
      ally[alias[j - 1] - 1] = alias[i - 1];
      j++;
    }
    i = j;
  }
}  /* sitecombine2 */


void sitescrunch2(long sites, long i, long j, steptr aliasweight)
{
  /* move so positively weighted sites come first */
  /* used by dnainvar, dnaml, dnamlk, & restml */
  long itemp;
  boolean done, found;

  done = false;
  while (!done) {
    if (aliasweight[i - 1] > 0)
      i++;
    else {
      if (j <= i)
        j = i + 1;
      if (j <= sites) {
        do {
          found = (aliasweight[j - 1] > 0);
          j++;
        } while (!(found || j > sites));
        if (found) {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = aliasweight[i - 1];
          aliasweight[i - 1] = aliasweight[j - 1];
          aliasweight[j - 1] = itemp;
        } else
          done = true;
      } else
        done = true;
    }
    done = (done || i >= sites);
  }
}  /* sitescrunch2 */


void makevalues(pointarray treenode, long *zeros, boolean usertree)
{
  /* set up fractional likelihoods at tips */
  /* used by dnacomp, dnapars, & dnapenny */
  long i, j;
  char ns = 0;
  node *p;

  setuptree(treenode, nonodes, usertree);
  for (i = 0; i < spp; i++)
    alloctip(treenode[i], zeros);
  if (!usertree) {
    for (i = spp; i < nonodes; i++) {
      p = treenode[i];
      do {
        allocnontip(p, zeros, endsite);
        p = p->next;
      } while (p != treenode[i]);
    }
  }
  for (j = 0; j < endsite; j++) {
    for (i = 0; i < spp; i++) {
      switch (y[i][alias[j] - 1]) {

      case 'A':
        ns = 1 << A;
        break;

      case 'C':
        ns = 1 << C;
        break;

      case 'G':
        ns = 1 << G;
        break;

      case 'U':
        ns = 1 << T;
        break;

      case 'T':
        ns = 1 << T;
        break;

      case 'M':
        ns = (1 << A) | (1 << C);
        break;

      case 'R':
        ns = (1 << A) | (1 << G);
        break;

      case 'W':
        ns = (1 << A) | (1 << T);
        break;

      case 'S':
        ns = (1 << C) | (1 << G);
        break;

      case 'Y':
        ns = (1 << C) | (1 << T);
        break;

      case 'K':
        ns = (1 << G) | (1 << T);
        break;

      case 'B':
        ns = (1 << C) | (1 << G) | (1 << T);
        break;

      case 'D':
        ns = (1 << A) | (1 << G) | (1 << T);
        break;

      case 'H':
        ns = (1 << A) | (1 << C) | (1 << T);
        break;

      case 'V':
        ns = (1 << A) | (1 << C) | (1 << G);
        break;

      case 'N':
        ns = (1 << A) | (1 << C) | (1 << G) | (1 << T);
        break;

      case 'X':
        ns = (1 << A) | (1 << C) | (1 << G) | (1 << T);
        break;

      case '?':
        ns = (1 << A) | (1 << C) | (1 << G) | (1 << T) | (1 << O);
        break;

      case 'O':
        ns = 1 << O;
        break;

      case '-':
        ns = 1 << O;
        break;
      }
      treenode[i]->base[j] = ns;
      treenode[i]->numsteps[j] = 0;
    }
  }
}  /* makevalues */


void makevalues2(long categs, pointarray treenode, long endsite,
                        long spp, sequence y, steptr alias)
{
  /* set up fractional likelihoods at tips */
  /* used by dnaml & dnamlk */
  long i, j, k, l;
  bases b;

  for (k = 0; k < endsite; k++) {
    j = alias[k];
    for (i = 0; i < spp; i++) {
      for (l = 0; l < categs; l++) {
        for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
          treenode[i]->x[k][l][(long)b - (long)A] = 0.0;
        switch (y[i][j - 1]) {

        case 'A':
          treenode[i]->x[k][l][0] = 1.0;
          break;

        case 'C':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          break;

        case 'G':
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'T':
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'U':
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'M':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          break;

        case 'R':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'W':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'S':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'Y':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'K':
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'B':
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'D':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'H':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)T - (long)A] = 1.0;
          break;

        case 'V':
          treenode[i]->x[k][l][0] = 1.0;
          treenode[i]->x[k][l][(long)C - (long)A] = 1.0;
          treenode[i]->x[k][l][(long)G - (long)A] = 1.0;
          break;

        case 'N':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case 'X':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case '?':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case 'O':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;

        case '-':
          for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
            treenode[i]->x[k][l][(long)b - (long)A] = 1.0;
          break;
        }
      }
    }
  }
}  /* makevalues2 */


void fillin(node *p, node *left, node *rt)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes.  */
  long i, j, k, n, purset, pyrset;
  node *q;

  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  if (!left) {
    memcpy(p->base, rt->base, endsite*sizeof(long));
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
    q = rt;
  } else if (!rt) {
    memcpy(p->base, left->base, endsite*sizeof(long));
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
    q = left;
  } else {
    for (i = 0; i < endsite; i++) {
      p->base[i] = left->base[i] & rt->base[i];
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (p->base[i] == 0) {
        p->base[i] = left->base[i] | rt->base[i];
        if (transvp) {
          if (!((p->base[i] == purset) || (p->base[i] == pyrset)))
            p->numsteps[i] += weight[i];
        }
        else p->numsteps[i] += weight[i];
      }
    }
    q = rt;
  }
  if (left && rt) n = 2;
  else n = 1;
  for (i = 0; i < endsite; i++)
    for (j = (long)A; j <= (long)O; j++)
      p->numnuc[i][j] = 0;
  for (k = 1; k <= n; k++) {
    if (k == 2) q = left;
    for (i = 0; i < endsite; i++) {
      for (j = (long)A; j <= (long)O; j++) {
        if (q->base[i] & (1 << j))
          p->numnuc[i][j]++;
      }
    }
  }
}  /* fillin */


long getlargest(long *numnuc)
{
  /* find the largest in array numnuc */
  long i, largest;

  largest = 0;
  for (i = (long)A; i <= (long)O; i++)
    if (numnuc[i] > largest)
      largest = numnuc[i];
  return largest;
} /* getlargest */


void multifillin(node *p, node *q, long dnumdesc)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes according to the
     changes in q's base */
  long i, j, b, largest, descsteps, purset, pyrset;

  memcpy(p->oldbase, p->base, endsite*sizeof(long));
  memcpy(p->oldnumsteps, p->numsteps, endsite*sizeof(long));
  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  for (i = 0; i < endsite; i++) {
    descsteps = 0;
    for (j = (long)A; j <= (long)O; j++) {
      b = 1 << j;
      if ((descsteps == 0) && (p->base[i] & b)) 
        descsteps = p->numsteps[i] 
                      - (p->numdesc - dnumdesc - p->numnuc[i][j]) * weight[i];
    }
    if (dnumdesc == -1)
      descsteps -= q->oldnumsteps[i];
    else if (dnumdesc == 0)
      descsteps += (q->numsteps[i] - q->oldnumsteps[i]);
    else
      descsteps += q->numsteps[i];
    if (q->oldbase[i] != q->base[i]) {
      for (j = (long)A; j <= (long)O; j++) {
        b = 1 << j;
        if (transvp) {
          if (b & purset) b = purset;
          if (b & pyrset) b = pyrset;
        }
        if ((q->oldbase[i] & b) && !(q->base[i] & b))
          p->numnuc[i][j]--;
        else if (!(q->oldbase[i] & b) && (q->base[i] & b))
          p->numnuc[i][j]++;
      }
    }
    largest = getlargest(p->numnuc[i]);
    if (q->oldbase[i] != q->base[i]) {
      p->base[i] = 0;
      for (j = (long)A; j <= (long)O; j++) {
        if (p->numnuc[i][j] == largest)
            p->base[i] |= (1 << j);
      }
    }
    p->numsteps[i] = (p->numdesc - largest) * weight[i] + descsteps;
  }
} /* multifillin */


void sumnsteps(node *p, node *left, node *rt, long a, long b)
{
  /* sets up for each node in the tree the base sequence
     at that point and counts the changes. */
  long i;
  long ns, rs, ls, purset, pyrset;

  if (!left) {
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
    memcpy(p->base, rt->base, endsite*sizeof(long));
  } else if (!rt) {
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
    memcpy(p->base, left->base, endsite*sizeof(long));
  } else  {
    purset = (1 << (long)A) + (1 << (long)G);
    pyrset = (1 << (long)C) + (1 << (long)T);
    for (i = a; i < b; i++) {
      ls = left->base[i];
      rs = rt->base[i];
      ns = ls & rs;
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (ns == 0) {
        ns = ls | rs;
        if (transvp) {
          if (!((ns == purset) || (ns == pyrset)))
            p->numsteps[i] += weight[i];
        }
        else p->numsteps[i] += weight[i];
      }
      p->base[i] = ns;
      }
    }
}  /* sumnsteps */


void sumnsteps2(node *p,node *left,node *rt,long a,long b,long *threshwt)
{
  /* counts the changes at each node.  */
  long i, steps;
  long ns, rs, ls, purset, pyrset;
  long term;

  if (a == 0) p->sumsteps = 0.0;
  if (!left)
    memcpy(p->numsteps, rt->numsteps, endsite*sizeof(long));
  else if (!rt)
    memcpy(p->numsteps, left->numsteps, endsite*sizeof(long));
  else {
    purset = (1 << (long)A) + (1 << (long)G);
    pyrset = (1 << (long)C) + (1 << (long)T);
    for (i = a; i < b; i++) {
      ls = left->base[i];
      rs = rt->base[i];
      ns = ls & rs;
      p->numsteps[i] = left->numsteps[i] + rt->numsteps[i];
      if (ns == 0) {
        ns = ls | rs;
        if (transvp) {
          if (!((ns == purset) || (ns == pyrset)))
            p->numsteps[i] += weight[i];
        }
        else p->numsteps[i] += weight[i];
      }
    }
  }
  for (i = a; i < b; i++) {
    steps = p->numsteps[i];
    if ((long)steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    p->sumsteps += (double)term;
  }
}  /* sumnsteps2 */


void multisumnsteps(node *p, node *q, long a, long b, long *threshwt)
{
  /* computes the number of steps between p and q */
  long i, j, steps, largest, descsteps, purset, pyrset, b1;
  long term;

  if (a == 0) p->sumsteps = 0.0;
  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  for (i = a; i < b; i++) {
    descsteps = 0;
    for (j = (long)A; j <= (long)O; j++) {
      if ((descsteps == 0) && (p->base[i] & (1 << j))) 
        descsteps = p->numsteps[i] -
                        (p->numdesc - 1 - p->numnuc[i][j]) * weight[i];
    }
    descsteps += q->numsteps[i];
    largest = 0;
    for (j = (long)A; j <= (long)O; j++) {
      b1 = (1 << j);
      if (transvp) {
        if (b1 & purset) b1 = purset;
        if (b1 & pyrset) b1 = pyrset;
      }
      if (q->base[i] & b1)
        p->numnuc[i][j]++;
      if (p->numnuc[i][j] > largest)
        largest = p->numnuc[i][j];
    }
    steps = (p->numdesc - largest) * weight[i] + descsteps;
    if ((long)steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    p->sumsteps += (double)term;
  }
} /* multisumnsteps */


void multisumnsteps2(node *p)
{
  /* counts the changes at each multi-way node. Sums up
     steps of all descendants */
  long i, j, largest, purset, pyrset, b1;
  node *q;
  baseptr b;

  purset = (1 << (long)A) + (1 << (long)G);
  pyrset = (1 << (long)C) + (1 << (long)T);
  for (i = 0; i < endsite; i++) {
    p->numsteps[i] = 0;
    q = p->next;
    while (q != p) {
      if (q->back) {
        p->numsteps[i] += q->back->numsteps[i];
        b = q->back->base;
        for (j = (long)A; j <= (long)O; j++) {
          b1 = (1 << j);   
          if (transvp) {
            if (b1 & purset) b1 = purset;
            if (b1 & pyrset) b1 = pyrset;
          }
          if (b[i] & b1) p->numnuc[i][j]++;
        }
      }
      q = q->next;
    }
    largest = getlargest(p->numnuc[i]);
    p->base[i] = 0;
    for (j = (long)A; j <= (long)O; j++) {
      if (p->numnuc[i][j] == largest)
        p->base[i] |= (1 << j);
    }
    p->numsteps[i] += ((p->numdesc - largest) * weight[i]);
  }
}  /* multisumnsteps2 */

boolean alltips(node *forknode, node *p)
{
  /* returns true if all descendants of forknode except p are tips; 
     false otherwise.  */
  node *q, *r;
  boolean tips;

  tips = true;
  r = forknode;
  q = forknode->next;
  do {
    if (q->back && q->back != p && !q->back->tip)
      tips = false;
    q = q->next;
  } while (tips && q != r);
  return tips;
} /* alltips */


void gdispose(node *p, node **grbg, pointarray treenode)
{
  /* go through tree throwing away nodes */
  node *q, *r;

  p->back = NULL;
  if (p->tip)
    return;
  treenode[p->index - 1] = NULL;
  q = p->next;
  while (q != p) {
    gdispose(q->back, grbg, treenode);
    q->back = NULL;
    r = q;
    q = q->next;
    chuck(grbg, r);
  }
  chuck(grbg, q);
}  /* gdispose */


void preorder(node *p, node *r, node *root, node *removing, node *adding,
                        node *changing, long dnumdesc)
{
  /* recompute number of steps in preorder taking both ancestoral and
     descendent steps into account. removing points to a node being 
     removed, if any */
  node *q, *p1, *p2;

  if (p && !p->tip && p != adding) {
    q = p;
    do {
      if (p->back != r) {
        if (p->numdesc > 2) {
          if (changing)
            multifillin (p, r, dnumdesc);
          else
            multifillin (p, r, 0);
        } else {
          p1 = p->next;
          if (!removing)
            while (!p1->back)
              p1 = p1->next;
          else
            while (!p1->back || p1->back == removing)
              p1 = p1->next;
          p2 = p1->next;
          if (!removing)
            while (!p2->back)
              p2 = p2->next;
          else
            while (!p2->back || p2->back == removing)
              p2 = p2->next;
          p1 = p1->back;
          p2 = p2->back;
          if (p->back == p1) p1 = NULL;
          else if (p->back == p2) p2 = NULL;
          memcpy(p->oldbase, p->base, endsite*sizeof(long));
          memcpy(p->oldnumsteps, p->numsteps, endsite*sizeof(long));
          fillin(p, p1, p2);
        }
      }
      p = p->next;
    } while (p != q);
    q = p;
    do {
      preorder(p->next->back, p->next, root, removing, adding, NULL, 0);
      p = p->next;
    } while (p->next != q);
  }
} /* preorder */


void updatenumdesc(node *p, node *root, long n)
{
  /* set p's numdesc to n.  If p is the root, numdesc of p's
  descendants are set to n-1. */
  node *q;

  q = p;
  if (p == root && n > 0) {
    p->numdesc = n;
    n--;
    q = q->next;
  }
  do {
    q->numdesc = n;
    q = q->next;
  } while (q != p);
}  /* updatenumdesc */


void add(node *below,node *newtip,node *newfork,node **root,
        boolean recompute,pointarray treenode,node **grbg,long *zeros)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant.
     if newfork is NULL, newtip is added as below's sibling */
  /* used in dnacomp & dnapars */
  node *p;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (newfork) {
    if (below->back != NULL)
      below->back->back = newfork;
    newfork->back = below->back;
    below->back = newfork->next->next;
    newfork->next->next->back = below;
    newfork->next->back = newtip;
    newtip->back = newfork->next;
    if (*root == below)
      *root = newfork;
    updatenumdesc(newfork, *root, 2);
  } else {
    gnutreenode(grbg, &p, below->index, endsite, zeros);
    p->back = newtip;
    newtip->back = p;
    p->next = below->next;
    below->next = p;
    updatenumdesc(below, *root, below->numdesc + 1);
  }
  if (!newtip->tip)
    updatenumdesc(newtip, *root, newtip->numdesc);
  (*root)->back = NULL;
  if (!recompute)
    return;
  if (!newfork) {
    memcpy(newtip->back->base, below->base, endsite*sizeof(long));
    memcpy(newtip->back->numsteps, below->numsteps, endsite*sizeof(long));
    memcpy(newtip->back->numnuc, below->numnuc, endsite*sizeof(nucarray));
    if (below != *root) {
      memcpy(below->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(below->back->oldnumsteps, zeros, endsite*sizeof(long));
      multifillin(newtip->back, below->back, 1);
    }
    if (!newtip->tip) {
      memcpy(newtip->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(newtip->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(newtip, newtip->back, *root, NULL, NULL, below, 1);
    }
    memcpy(newtip->oldbase, zeros, endsite*sizeof(long));
    memcpy(newtip->oldnumsteps, zeros, endsite*sizeof(long));
    preorder(below, newtip, *root, NULL, newtip, below, 1);
    if (below != *root)
      preorder(below->back, below, *root, NULL, NULL, NULL, 0);
  } else {
    fillin(newtip->back, newtip->back->next->back,
             newtip->back->next->next->back);
    if (!newtip->tip) {
      memcpy(newtip->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(newtip->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(newtip, newtip->back, *root, NULL, NULL, newfork, 1);
    }
    if (newfork != *root) {
      memcpy(below->back->base, newfork->back->base, endsite*sizeof(long));
      memcpy(below->back->numsteps, newfork->back->numsteps, endsite*sizeof(long));
      preorder(newfork, newtip, *root, NULL, newtip, NULL, 0);
    } else {
      fillin(below->back, newtip, NULL);
      fillin(newfork, newtip, below);
      memcpy(below->back->oldbase, zeros, endsite*sizeof(long));
      memcpy(below->back->oldnumsteps, zeros, endsite*sizeof(long));
      preorder(below, below->back, *root, NULL, NULL, newfork, 1);
    }
    if (newfork != *root) {
      memcpy(newfork->oldbase, below->base, endsite*sizeof(long));
      memcpy(newfork->oldnumsteps, below->numsteps, endsite*sizeof(long));
      preorder(newfork->back, newfork, *root, NULL, NULL, NULL, 0);
    }
  }
}  /* add */


void findbelow(node **below, node *item, node *fork)
{
  /* decide which of fork's binary children is below */

  if (fork->next->back == item)
    *below = fork->next->next->back;
  else
    *below = fork->next->back;
} /* findbelow */


void re_move(node *item, node **fork, node **root, boolean recompute,
                        pointarray treenode, node **grbg, long *zeros)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork.
     If item belongs to a node with more than 2 descendants,
     fork will not be deleted */
  /* used in dnacomp & dnapars */
  node *p, *q, *other = NULL, *otherback = NULL;

  if (item->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[item->back->index - 1];
  if ((*fork)->numdesc == 2) {
    updatenumdesc(*fork, *root, 0);
    findbelow(&other, item, *fork);
    otherback = other->back;
    if (*root == *fork) {
      *root = other;
      if (!other->tip)
        updatenumdesc(other, *root, other->numdesc);
    }
    p = item->back->next->back;
    q = item->back->next->next->back;
    if (p != NULL)
      p->back = q;
    if (q != NULL)
      q->back = p;
    (*fork)->back = NULL;
    p = (*fork)->next;
    while (p != *fork) {
      p->back = NULL;
      p = p->next;
    }
  } else {
    updatenumdesc(*fork, *root, (*fork)->numdesc - 1);
    p = *fork;
    while (p->next != item->back)
      p = p->next;
    p->next = item->back->next;
  }
  if (!item->tip) {
    updatenumdesc(item, item, item->numdesc);
    if (recompute) {
      memcpy(item->back->oldbase, item->back->base, endsite*sizeof(long));
      memcpy(item->back->oldnumsteps, item->back->numsteps, endsite*sizeof(long));
      memcpy(item->back->base, zeros, endsite*sizeof(long));
      memcpy(item->back->numsteps, zeros, endsite*sizeof(long));
      preorder(item, item->back, *root, item->back, NULL, item, -1);
    }
  }
  if ((*fork)->numdesc >= 2)
    chuck(grbg, item->back);
  item->back = NULL;
  if (!recompute)
    return;
  if ((*fork)->numdesc == 0) {
    memcpy(otherback->oldbase, otherback->base, endsite*sizeof(long));  
    memcpy(otherback->oldnumsteps, otherback->numsteps, endsite*sizeof(long));
    if (other == *root) {
      memcpy(otherback->base, zeros, endsite*sizeof(long));
      memcpy(otherback->numsteps, zeros, endsite*sizeof(long));
    } else {
      memcpy(otherback->base, other->back->base, endsite*sizeof(long));
      memcpy(otherback->numsteps, other->back->numsteps, endsite*sizeof(long));
    }
    p = other->back;
    other->back = otherback;
    if (other == *root)
      preorder(other, otherback, *root, otherback, NULL, other, -1);
    else
      preorder(other, otherback, *root, NULL, NULL, NULL, 0);
    other->back = p;
    if (other != *root) {
      memcpy(other->oldbase,(*fork)->base, endsite*sizeof(long));
      memcpy(other->oldnumsteps,(*fork)->numsteps, endsite*sizeof(long));
      preorder(other->back, other, *root, NULL, NULL, NULL, 0);
    }
  } else {
    memcpy(item->oldbase, item->base, endsite*sizeof(long));
    memcpy(item->oldnumsteps, item->numsteps, endsite*sizeof(long));
    memcpy(item->base, zeros, endsite*sizeof(long));
    memcpy(item->numsteps, zeros, endsite*sizeof(long));
    preorder(*fork, item, *root, NULL, NULL, *fork, -1);
    if (*fork != *root)
      preorder((*fork)->back, *fork, *root, NULL, NULL, NULL, 0);
    memcpy(item->base, item->oldbase, endsite*sizeof(long));
    memcpy(item->numsteps, item->oldnumsteps, endsite*sizeof(long));
  }
}  /* remove */


void postorder(node *p)
{
  /* traverses an n-ary tree, suming up steps at a node's descendants */
  /* used in dnacomp, dnapars, & dnapenny */
  node *q;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    postorder(q->back);
    q = q->next;
  }
  zeronumnuc(p, endsite);
  if (p->numdesc > 2)
    multisumnsteps2(p);
  else
    fillin(p, p->next->back, p->next->next->back);
}  /* postorder */


void getnufork(node **nufork,node **grbg,pointarray treenode,long *zeros)
{
  /* find a fork not used currently */
  long i;

  i = spp;
  while (treenode[i] && treenode[i]->numdesc > 0) i++;
  if (!treenode[i])
    gnutreenode(grbg, &treenode[i], i, endsite, zeros);
  *nufork = treenode[i];
} /* getnufork */


void reroot(node *outgroup, node *root)
{
  /* reorients tree, putting outgroup in desired position. used if
     the root is binary. */
  /* used in dnacomp & dnapars */
  node *p, *q;

  if (outgroup->back->index == root->index)
    return;
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = q;
  outgroup->back = p;
}  /* reroot */


void reroot2(node *outgroup, node *root)
{
  /* reorients tree, putting outgroup in desired position. */
  /* used in dnacomp & dnapars */
  node *p;

  p = outgroup->back->next;
  while (p->next != outgroup->back)
    p = p->next;
  root->next = outgroup->back;
  p->next = root;
}  /* reroot2 */


void reroot3(node *outgroup, node *root, node *root2, node *lastdesc,
                        node **grbg)
{
  /* reorients tree, putting back outgroup in original position. */
  /* used in dnacomp & dnapars */
  node *p;

  p = root->next;
  while (p->next != root)
    p = p->next;
  chuck(grbg, root);
  p->next = outgroup->back;
  root2->next = lastdesc->next;
  lastdesc->next = root2;
}  /* reroot3 */


void savetraverse(node *p)
{
  /* sets BOOLEANs that indicate which way is down */
  node *q;

  p->bottom = true;
  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    q->bottom = false;
    savetraverse(q->back);
    q = q->next;
  }
}  /* savetraverse */


void newindex(long i, node *p)
{
  /* assigns index i to node p */

  while (p->index != i) {
    p->index = i;
    p = p->next;
  }
} /* newindex */


void flipindexes(long nextnode, pointarray treenode)
{
  /* flips index of nodes between nextnode and last node.  */
  long last;
  node *temp;

  last = nonodes;
  while (treenode[last - 1]->numdesc == 0)
    last--;
  if (last > nextnode) {
    temp = treenode[nextnode - 1];
    treenode[nextnode - 1] = treenode[last - 1];
    treenode[last - 1] = temp;
    newindex(nextnode, treenode[nextnode - 1]);
    newindex(last, treenode[last - 1]);
  }
} /* flipindexes */  


boolean parentinmulti(node *anode)
{
  /* sees if anode's parent has more than 2 children */
  node *p;

  while (!anode->bottom) anode = anode->next;
  p = anode->back;
  while (!p->bottom)
    p = p->next;
  return (p->numdesc > 2);
} /* parentinmulti */


long sibsvisited(node *anode, long *place)
{
  /* computes the number of nodes which are visited earlier than anode among 
  its siblings */
  node *p;
  long nvisited;

  while (!anode->bottom) anode = anode->next;
  p = anode->back->next;
  nvisited = 0;
  do {
    if (!p->bottom && place[p->back->index - 1] != 0)
      nvisited++;
    p = p->next;
  } while (p != anode->back);
  return nvisited;
}  /* sibsvisited */


long smallest(node *anode, long *place)
{
  /* finds the smallest index of sibling of anode */
  node *p;
  long min;

  while (!anode->bottom) anode = anode->next;
  p = anode->back->next;
  if (p->bottom) p = p->next;
  min = nonodes;
  do {
    if (p->back && place[p->back->index - 1] != 0) {
      if (p->back->index <= spp) {
        if (p->back->index < min)
          min = p->back->index;
      } else {
        if (place[p->back->index - 1] < min)
          min = place[p->back->index - 1];
      }
    }
    p = p->next;
    if (p->bottom) p = p->next;
  } while (p != anode->back);
  return min;
}  /* smallest */


void bintomulti(node **root, node **binroot, node **grbg, long *zeros)
{  /* attaches root's left child to its right child and makes
      the right child new root */
  node *left, *right, *newnode, *temp;

  right = (*root)->next->next->back;
  left = (*root)->next->back;
  if (right->tip) {
    (*root)->next = right->back;
    (*root)->next->next = left->back;
    temp = left;
    left = right;
    right = temp;
    right->back->next = *root;
  }
  gnutreenode(grbg, &newnode, right->index, endsite, zeros);
  newnode->next = right->next;
  newnode->back = left;
  left->back = newnode;
  right->next = newnode;
  (*root)->next->back = (*root)->next->next->back = NULL;
  *binroot = *root;
  (*binroot)->numdesc = 0;
  *root = right;
  (*root)->numdesc++;
  (*root)->back = NULL;
} /* bintomulti */


void backtobinary(node **root, node *binroot, node **grbg)
{ /* restores binary root */
  node *p;

  binroot->next->back = (*root)->next->back;
  (*root)->next->back->back = binroot->next;
  p = (*root)->next;
  (*root)->next = p->next;
  binroot->next->next->back = *root;
  (*root)->back = binroot->next->next;
  chuck(grbg, p);
  (*root)->numdesc--;
  *root = binroot;
  (*root)->numdesc = 2;
} /* backtobinary */


boolean outgrin(node *root, node *outgrnode)
{ /* checks if outgroup node is a child of root */
  node *p;

  p = root->next;
  while (p != root) {
    if (p->back == outgrnode)
      return true;
    p = p->next;
  }
  return false;
} /* outgrin */


void flipnodes(node *nodea, node *nodeb)
{ /* flip nodes */
  node *backa, *backb;

  backa = nodea->back;
  backb = nodeb->back;
  backa->back = nodeb;
  backb->back = nodea;
  nodea->back = backb;
  nodeb->back = backa;
} /* flipnodes */


void moveleft(node *root, node *outgrnode, node **flipback)
{ /* makes outgroup node to leftmost child of root */
  node *p;
  boolean done;

  p = root->next;
  done = false;
  while (p != root && !done) {
    if (p->back == outgrnode) {
      *flipback = p;
      flipnodes(root->next->back, p->back);
      done = true;
    }
    p = p->next;
  }
} /* moveleft */


void savetree(node *root,  long *place, pointarray treenode,
                        node **grbg, long *zeros)
{ /* record in place where each species has to be
     added to reconstruct this tree */
  /* used by dnacomp & dnapars */
  long i, j, nextnode, nvisited;
  node *p, *q, *r = NULL, *root2, *lastdesc, 
                *outgrnode, *binroot, *flipback;
  boolean done, newfork;

  binroot = NULL;
  lastdesc = NULL;
  root2 = NULL;
  flipback = NULL;
  outgrnode = treenode[outgrno - 1];
  if (root->numdesc == 2)
    bintomulti(&root, &binroot, grbg, zeros);
  if (outgrin(root, outgrnode)) {
    if (outgrnode != root->next->back)
      moveleft(root, outgrnode, &flipback);
  } else {
    root2 = root;
    lastdesc = root->next;
    while (lastdesc->next != root)
      lastdesc = lastdesc->next;
    lastdesc->next = root->next;
    gnutreenode(grbg, &root, outgrnode->back->index, endsite, zeros);
    root->numdesc = root2->numdesc;
    reroot2(outgrnode, root);
  }
  savetraverse(root);
  nextnode = spp + 1;
  for (i = nextnode; i <= nonodes; i++)
    if (treenode[i - 1]->numdesc == 0)
      flipindexes(i, treenode);
  for (i = 0; i < nonodes; i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= spp; i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      while (!p->bottom)
        p = p->next;
      r = p;
      p = p->back;
    }
    if (i > 1) {
      q = treenode[i - 1]; 
      newfork = true;
      nvisited = sibsvisited(q, place);
      if (nvisited == 0) {
        if (parentinmulti(r)) {
          nvisited = sibsvisited(r, place);
          if (nvisited == 0)
            place[i - 1] = place[p->index - 1];
          else if (nvisited == 1)
            place[i - 1] = smallest(r, place);
          else {
            place[i - 1] = -smallest(r, place);
            newfork = false;
          }
        } else
          place[i - 1] = place[p->index - 1];
      } else if (nvisited == 1) {
        place[i - 1] = place[p->index - 1];
      } else {
        place[i - 1] = -smallest(q, place);
        newfork = false;
      }
      if (newfork) {
        j = place[p->index - 1];
        done = false;
        while (!done) {
          place[p->index - 1] = nextnode;
          while (!p->bottom)
            p = p->next;
          p = p->back;
          done = (p == NULL);
          if (!done)
            done = (place[p->index - 1] != j);
          if (done) {
            nextnode++;
          }
        }
      }
    }
  }
  if (flipback)
    flipnodes(outgrnode, flipback->back);
  else {
    if (root2) {
      reroot3(outgrnode, root, root2, lastdesc, grbg);
      root = root2;
    }
  }
  if (binroot)
    backtobinary(&root, binroot, grbg);
}  /* savetree */ 


void addnsave(node *p, node *item, node *nufork, node **root, node **grbg,
                boolean multf, pointarray treenode, long *place, long *zeros)
{  /* adds item to tree and save it.  Then removes item. */
  node *dummy;

  if (!multf)
    add(p, item, nufork, root, false, treenode, grbg, zeros);
  else
    add(p, item, NULL, root, false, treenode, grbg, zeros);
  savetree(*root, place, treenode, grbg, zeros);
  if (!multf)
    re_move(item, &nufork, root, false, treenode, grbg, zeros);
  else
    re_move(item, &dummy, root, false, treenode, grbg, zeros);
} /* addnsave */


void addbestever(long *pos, long *nextree, long maxtrees, boolean collapse,
                        long *place, bestelm *bestrees)
{ /* adds first best tree */

  *pos = 1;
  *nextree = 1;
  initbestrees(bestrees, maxtrees, true);
  initbestrees(bestrees, maxtrees, false);
  addtree(*pos, nextree, collapse, place, bestrees);
} /* addbestever */


void addtiedtree(long pos, long *nextree, long maxtrees, boolean collapse,
                        long *place, bestelm *bestrees)
{ /* add tied tree */

  if (*nextree <= maxtrees)
    addtree(pos, nextree, collapse, place, bestrees);
} /* addtiedtree */


void clearcollapse(pointarray treenode)
{
  /* clears collapse status at a node */
  long i;
  node *p;

  for (i = 0; i < nonodes; i++) {
    treenode[i]->collapse = undefined;
    if (!treenode[i]->tip) {
      p = treenode[i]->next;
      while (p != treenode[i]) {
        p->collapse = undefined;
        p = p->next;
      }
    }
  }
} /* clearcollapse */


void clearbottom(pointarray treenode)
{
  /* clears boolean bottom at a node */
  long i;
  node *p;

  for (i = 0; i < nonodes; i++) {
    treenode[i]->bottom = false;
    if (!treenode[i]->tip) {
      p = treenode[i]->next;
      while (p != treenode[i]) {
        p->bottom = false;
        p = p->next;
      }
    }
  }
} /* clearbottom */


void collabranch(node *collapfrom, node *tempfrom, node *tempto)
{ /* collapse branch from collapfrom */
  long i, j, b, largest, descsteps;
  boolean done;

  for (i = 0; i < endsite; i++) {
    descsteps = 0;
    for (j = (long)A; j <= (long)O; j++) {
      b = 1 << j;
      if ((descsteps == 0) && (collapfrom->base[i] & b)) 
        descsteps = tempfrom->oldnumsteps[i] 
                     - (collapfrom->numdesc - collapfrom->numnuc[i][j])
                       * weight[i];
    }
    done = false;
    for (j = (long)A; j <= (long)O; j++) {
      b = 1 << j;
      if (!done && (tempto->base[i] & b)) {
        descsteps += (tempto->numsteps[i] 
                      - (tempto->numdesc - collapfrom->numdesc
                        - tempto->numnuc[i][j]) * weight[i]);
        done = true;
      }
    }
    for (j = (long)A; j <= (long)O; j++)
      tempto->numnuc[i][j] += collapfrom->numnuc[i][j];
    largest = getlargest(tempto->numnuc[i]);
    tempto->base[i] = 0;
    for (j = (long)A; j <= (long)O; j++) {
      if (tempto->numnuc[i][j] == largest)
        tempto->base[i] |= (1 << j);
    }
    tempto->numsteps[i] = (tempto->numdesc - largest) * weight[i] + descsteps;
  }
} /* collabranch */


boolean allcommonbases(node *a, node *b, boolean *allsame)
{  /* see if bases are common at all sites for nodes a and b */    
  long i;
  boolean allcommon;

  allcommon = true;
  *allsame = true;
  for (i = 0; i < endsite; i++) {
    if ((a->base[i] & b->base[i]) == 0)
      allcommon = false;
    else if (a->base[i] != b->base[i])
      *allsame = false;
  }
  return allcommon;
} /* allcommonbases */


void findbottom(node *p, node **bottom)
{ /* find a node with field bottom set at node p */
  node *q;

  if (p->bottom)
    *bottom = p;
  else {
    q = p->next;
    while(!q->bottom && q != p)
      q = q->next;
    *bottom = q;
  }
} /* findbottom */


boolean moresteps(node *a, node *b)
{  /* see if numsteps of node a exceeds those of node b */    
  long i;

  for (i = 0; i < endsite; i++)
    if (a->numsteps[i] > b->numsteps[i])
      return true;
  return false;
} /* moresteps */


boolean passdown(node *desc, node *parent, node *start, node *below,
                        node *item, node *added, node *total, node *tempdsc,
            node *tempprt, boolean multf)
{ /* track down to node start to see if an ancestor branch can be collapsed */
  node *temp;
  boolean done, allsame;

  done = (parent == start);
  while (!done) {
    desc = parent;
    findbottom(parent->back, &parent);
    if (multf && start == below && parent == below)
      parent = added;
    memcpy(tempdsc->base, tempprt->base, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, desc->base, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, desc->numsteps, endsite*sizeof(long));
    memcpy(tempprt->base, parent->base, endsite*sizeof(long));
    memcpy(tempprt->numsteps, parent->numsteps, endsite*sizeof(long));
    memcpy(tempprt->numnuc, parent->numnuc, endsite*sizeof(nucarray));
    tempprt->numdesc = parent->numdesc;
    multifillin(tempprt, tempdsc, 0);
    if (!allcommonbases(tempprt, parent, &allsame))
      return false;
    else if (moresteps(tempprt, parent))
      return false;
    else if (allsame)
      return true;
    if (parent == added)
      parent = below;
    done = (parent == start);
    if (done && ((start == item) || (!multf && start == below))) {
      memcpy(tempdsc->base, tempprt->base, endsite*sizeof(long));
      memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
      memcpy(tempdsc->oldbase, start->base, endsite*sizeof(long));
      memcpy(tempdsc->oldnumsteps, start->numsteps, endsite*sizeof(long));
      multifillin(added, tempdsc, 0);
      tempprt = added;
    }
  } 
  temp = tempdsc;
  if (start == below || start == item)
    fillin(temp, tempprt, below->back);
  else
    fillin(temp, tempprt, added);
  return !moresteps(temp, total);
} /* passdown */


boolean trycollapdesc(node *desc, node *parent, node *start,
                        node *below, node *item, node *added, node *total,
            node *tempdsc, node *tempprt, boolean multf, long *zeros)
  { /* see if branch between nodes desc and parent can be collapsed */
  boolean allsame;

  if (desc->numdesc == 1)
    return true;
  if (multf && start == below && parent == below)
    parent = added;
  memcpy(tempdsc->base, zeros, endsite*sizeof(long));
  memcpy(tempdsc->numsteps, zeros, endsite*sizeof(long));
  memcpy(tempdsc->oldbase, desc->base, endsite*sizeof(long));
  memcpy(tempdsc->oldnumsteps, desc->numsteps, endsite*sizeof(long));
  memcpy(tempprt->base, parent->base, endsite*sizeof(long));
  memcpy(tempprt->numsteps, parent->numsteps, endsite*sizeof(long));
  memcpy(tempprt->numnuc, parent->numnuc, endsite*sizeof(nucarray));
  tempprt->numdesc = parent->numdesc - 1;
  multifillin(tempprt, tempdsc, -1);
  tempprt->numdesc += desc->numdesc;
  collabranch(desc, tempdsc, tempprt);
  if (!allcommonbases(tempprt, parent, &allsame) || 
        moresteps(tempprt, parent)) {
    if (parent != added) {
      desc->collapse = nocollap;
      parent->collapse = nocollap;
    }
    return false;
  } else if (allsame) {
    if (parent != added) {
      desc->collapse = tocollap;
      parent->collapse = tocollap;
    }
    return true;
  }
  if (parent == added)
    parent = below;
  if ((start == item && parent == item) ||
        (!multf && start == below && parent == below)) {
    memcpy(tempdsc->base, tempprt->base, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, tempprt->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, start->base, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, start->numsteps, endsite*sizeof(long));
    memcpy(tempprt->base, added->base, endsite*sizeof(long));
    memcpy(tempprt->numsteps, added->numsteps, endsite*sizeof(long));
    memcpy(tempprt->numnuc, added->numnuc, endsite*sizeof(nucarray));
    tempprt->numdesc = added->numdesc;
    multifillin(tempprt, tempdsc, 0);
    if (!allcommonbases(tempprt, added, &allsame))
      return false;
    else if (moresteps(tempprt, added))
      return false;
    else if (allsame)
      return true;
  }
  return passdown(desc, parent, start, below, item, added, total, tempdsc,
                    tempprt, multf);
} /* trycollapdesc */


void setbottom(node *p)
{ /* set field bottom at node p */
  node *q;

  p->bottom = true;
  q = p->next;
  do {
    q->bottom = false;
    q = q->next;
  } while (q != p);
} /* setbottom */

boolean zeroinsubtree(node *subtree, node *start, node *below, node *item,
                        node *added, node *total, node *tempdsc, node *tempprt,
                        boolean multf, node* root, long *zeros)
{ /* sees if subtree contains a zero length branch */
  node *p;

  if (!subtree->tip) {
    setbottom(subtree);
    p = subtree->next;
    do {
      if (p->back && !p->back->tip && 
         !((p->back->collapse == nocollap) && (subtree->collapse == nocollap))
           && (subtree->numdesc != 1)) {
        if ((p->back->collapse == tocollap) && (subtree->collapse == tocollap)
            && multf && (subtree != below))
          return true;
        /* when root->numdesc == 2
         * there is no mandatory step at the root, 
         * instead of checking at the root we check around it 
         * we only need to check p because the first if 
         * statement already gets rid of it for the subtree */
        else if ((p->back->index != root->index || root->numdesc > 2) && 
            trycollapdesc(p->back, subtree, start, below, item, added, total, 
                tempdsc, tempprt, multf, zeros))
          return true;
        else if ((p->back->index == root->index && root->numdesc == 2) && 
            !(root->next->back->tip) && !(root->next->next->back->tip) && 
            trycollapdesc(root->next->back, root->next->next->back, start, 
                below, item,added, total, tempdsc, tempprt, multf, zeros))
          return true;
      }
      p = p->next;
    } while (p != subtree);
    p = subtree->next;
    do {
      if (p->back && !p->back->tip) {
        if (zeroinsubtree(p->back, start, below, item, added, total, 
                            tempdsc, tempprt, multf, root, zeros))
          return true;
      }
      p = p->next;
    } while (p != subtree);
  }
  return false;
} /* zeroinsubtree */


boolean collapsible(node *item, node *below, node *temp, node *temp1,
                        node *tempdsc, node *tempprt, node *added, node *total,
            boolean multf, node *root, long *zeros, pointarray treenode)
{
  /* sees if any branch can be collapsed */
  node *belowbk;
  boolean allsame;

  if (multf) {
    memcpy(tempdsc->base, item->base, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, item->numsteps, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, zeros, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy(added->base, below->base, endsite*sizeof(long));
    memcpy(added->numsteps, below->numsteps, endsite*sizeof(long));
    memcpy(added->numnuc, below->numnuc, endsite*sizeof(nucarray));
    added->numdesc = below->numdesc + 1;
    multifillin(added, tempdsc, 1);
  } else {
    fillin(added, item, below);
    added->numdesc = 2;
  }
  fillin(total, added, below->back);
  clearbottom(treenode);
  if (below->back) {
    if (zeroinsubtree(below->back, below->back, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros))
      return true;
  }
  if (multf) {
    if (zeroinsubtree(below, below, below, item, added, total,
                       tempdsc, tempprt, multf, root, zeros))
      return true;
  } else if (!below->tip) {
    if (zeroinsubtree(below, below, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros))
      return true;
  }
  if (!item->tip) {
    if (zeroinsubtree(item, item, below, item, added, total,
                        tempdsc, tempprt, multf, root, zeros))
      return true;
  }
  if (multf && below->back && !below->back->tip) {
    memcpy(tempdsc->base, zeros, endsite*sizeof(long));
    memcpy(tempdsc->numsteps, zeros, endsite*sizeof(long));
    memcpy(tempdsc->oldbase, added->base, endsite*sizeof(long));
    memcpy(tempdsc->oldnumsteps, added->numsteps, endsite*sizeof(long));
    if (below->back == treenode[below->back->index - 1])
      belowbk = below->back->next;
    else
      belowbk = treenode[below->back->index - 1];
    memcpy(tempprt->base, belowbk->base, endsite*sizeof(long));
    memcpy(tempprt->numsteps, belowbk->numsteps, endsite*sizeof(long));
    memcpy(tempprt->numnuc, belowbk->numnuc, endsite*sizeof(nucarray));
    tempprt->numdesc = belowbk->numdesc - 1;
    multifillin(tempprt, tempdsc, -1);
    tempprt->numdesc += added->numdesc;
    collabranch(added, tempdsc, tempprt);
    if (!allcommonbases(tempprt, belowbk, &allsame))
      return false;
    else if (allsame && !moresteps(tempprt, belowbk))
      return true;
    else if (belowbk->back) {
      fillin(temp, tempprt, belowbk->back);
      fillin(temp1, belowbk, belowbk->back);
      return !moresteps(temp, temp1);
    }
  }
  return false;
} /* collapsible */


void replaceback(node **oldback, node *item, node *forknode,
                        node **grbg, long *zeros)
{ /* replaces back node of item with another */
  node *p;

  p = forknode;
  while (p->next->back != item)
    p = p->next;
  *oldback = p->next;
  gnutreenode(grbg, &p->next, forknode->index, endsite, zeros);
  p->next->next = (*oldback)->next;
  p->next->back = (*oldback)->back;
  p->next->back->back = p->next;
  (*oldback)->next = (*oldback)->back = NULL;
} /* replaceback */


void putback(node *oldback, node *item, node *forknode, node **grbg)
{ /* restores node to back of item */
  node *p, *q;

  p = forknode;
  while (p->next != item->back)
    p = p->next;
  q = p->next;
  oldback->next = p->next->next;
  p->next = oldback;
  oldback->back = item;
  item->back = oldback;
  oldback->index = forknode->index;
  chuck(grbg, q);
} /* putback */


void savelocrearr(node *item, node *forknode, node *below, node *tmp,
        node *tmp1, node *tmp2, node *tmp3, node *tmprm, node *tmpadd,
        node **root, long maxtrees, long *nextree, boolean multf,
        boolean bestever, boolean *saved, long *place,
        bestelm *bestrees, pointarray treenode, node **grbg,
        long *zeros)
{ /* saves tied or better trees during local rearrangements by removing
     item from forknode and adding to below */
  node *other, *otherback = NULL, *oldfork, *nufork, *oldback;
  long pos;
  boolean found, collapse;

  if (forknode->numdesc == 2) {
    findbelow(&other, item, forknode);
    otherback = other->back;
    oldback = NULL;
  } else {
    other = NULL;
    replaceback(&oldback, item, forknode, grbg, zeros);
  }
  re_move(item, &oldfork, root, false, treenode, grbg, zeros);
  if (!multf)
    getnufork(&nufork, grbg, treenode, zeros);
  else
    nufork = NULL;
  addnsave(below, item, nufork, root, grbg, multf, treenode, place, zeros);
  pos = 0;
  findtree(&found, &pos, *nextree, place, bestrees);
  if (other) {
    add(other, item, oldfork, root, false, treenode, grbg, zeros);
    if (otherback->back != other)
      flipnodes(item, other);
  } else
    add(forknode, item, NULL, root, false, treenode, grbg, zeros);
  *saved = false;
  if (found) {
    if (oldback)
      putback(oldback, item, forknode, grbg);
  } else {
    if (oldback)
      chuck(grbg, oldback);
    re_move(item, &oldfork, root, true, treenode, grbg, zeros);
    collapse = collapsible(item, below, tmp, tmp1, tmp2, tmp3, tmprm,
                     tmpadd, multf, *root, zeros, treenode);
    if (!collapse) {
      if (bestever)
        addbestever(&pos, nextree, maxtrees, collapse, place, bestrees);
      else
        addtiedtree(pos, nextree, maxtrees, collapse, place, bestrees);
    }
    if (other)
      add(other, item, oldfork, root, true, treenode, grbg, zeros);
    else
      add(forknode, item, NULL, root, true, treenode, grbg, zeros);
    *saved = !collapse;
  }
} /* savelocrearr */


void clearvisited(pointarray treenode)
{
  /* clears boolean visited at a node */
  long i;
  node *p;

  for (i = 0; i < nonodes; i++) {
    treenode[i]->visited = false;
    if (!treenode[i]->tip) {
      p = treenode[i]->next;
      while (p != treenode[i]) {
        p->visited = false;
        p = p->next;
      }
    }
  }
} /* clearvisited */


void hyprint(long b1, long b2, struct LOC_hyptrav *htrav,
                        pointarray treenode, Char *basechar)
{
  /* print out states in sites b1 through b2 at node */
  long i, j, k, n;
  boolean dot;
  bases b;

  if (htrav->bottom) {
    if (!outgropt)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  } else
    fprintf(outfile, "%4ld   ", htrav->r->back->index - spp);
  if (htrav->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[htrav->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", htrav->r->index - spp);
  if (htrav->bottom)
    fprintf(outfile, "          ");
  else if (htrav->nonzero)
    fprintf(outfile, "   yes    ");
  else if (htrav->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1; i <= b2; i++) {
    j = location[ally[i - 1] - 1];
    htrav->tempset = htrav->r->base[j - 1];
    htrav->anc = htrav->hypset[j - 1];
    if (!htrav->bottom)
      htrav->anc = treenode[htrav->r->back->index - 1]->base[j - 1];
    dot = dotdiff && (htrav->tempset == htrav->anc && !htrav->bottom);
    if (dot)
      putc('.', outfile); 
    else if (htrav->tempset == (1 << A))
      putc('A', outfile);
    else if (htrav->tempset == (1 << C))
      putc('C', outfile);
    else if (htrav->tempset == (1 << G))
      putc('G', outfile);
    else if (htrav->tempset == (1 << T))
      putc('T', outfile);
    else if (htrav->tempset == (1 << O))
      putc('-', outfile);
    else {
      k = 1;
      n = 0;
      for (b = A; b <= O; b = b + 1) {
        if (((1 << b) & htrav->tempset) != 0)
          n += k;
        k += k;
      }
      putc(basechar[n - 1], outfile);
    }
    if (i % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */


void gnubase(gbases **p, gbases **garbage, long endsite)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (*garbage != NULL) {
    *p = *garbage;
    *garbage = (*garbage)->next;
  } else {
    *p = (gbases *)Malloc(sizeof(gbases));
    (*p)->base = (baseptr)Malloc(endsite*sizeof(long));
  }
  (*p)->next = NULL;
}  /* gnubase */


void chuckbase(gbases *p, gbases **garbage)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = *garbage;
  *garbage = p;
}  /* chuckbase */


void hyptrav(node *r_, long *hypset_, long b1, long b2, boolean bottom_,
                        pointarray treenode, gbases **garbage, Char *basechar)
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  long i, j, k;
  long largest;
  gbases *ancset;
  nucarray *tempnuc;
  node *p, *q;

  Vars.bottom = bottom_;
  Vars.r = r_;
  Vars.hypset = hypset_;
  gnubase(&ancset, garbage, endsite);
  tempnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  Vars.maybe = false;
  Vars.nonzero = false;
  if (!Vars.r->tip)
    zeronumnuc(Vars.r, endsite);
  for (i = b1 - 1; i < b2; i++) {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip) {
      p = Vars.r->next;
      for (k = (long)A; k <= (long)O; k++)
        if (Vars.anc & (1 << k))
          Vars.r->numnuc[j - 1][k]++;
      do {
        for (k = (long)A; k <= (long)O; k++)
          if (p->back->base[j - 1] & (1 << k))
            Vars.r->numnuc[j - 1][k]++;
        p = p->next;
      } while (p != Vars.r);
      largest = getlargest(Vars.r->numnuc[j - 1]);
      Vars.tempset = 0;
      for (k = (long)A; k <= (long)O; k++) {
        if (Vars.r->numnuc[j - 1][k] == largest)
          Vars.tempset |= (1 << k);
      }
      Vars.r->base[j - 1] = Vars.tempset;
    }
    if (!Vars.bottom)
      Vars.anc = treenode[Vars.r->back->index - 1]->base[j - 1];
    Vars.nonzero = (Vars.nonzero || (Vars.r->base[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || Vars.r->base[j - 1] != Vars.anc);
  }
  hyprint(b1, b2, &Vars, treenode, basechar);
  Vars.bottom = false;
  if (!Vars.r->tip) {
    memcpy(tempnuc, Vars.r->numnuc, endsite*sizeof(nucarray));
    q = Vars.r->next;
    do {
      memcpy(Vars.r->numnuc, tempnuc, endsite*sizeof(nucarray));
      for (i = b1 - 1; i < b2; i++) {
        j = location[ally[i] - 1];
        for (k = (long)A; k <= (long)O; k++)
          if (q->back->base[j - 1] & (1 << k))
            Vars.r->numnuc[j - 1][k]--;
        largest = getlargest(Vars.r->numnuc[j - 1]);
        ancset->base[j - 1] = 0;
        for (k = (long)A; k <= (long)O; k++)
          if (Vars.r->numnuc[j - 1][k] == largest)
            ancset->base[j - 1] |= (1 << k);
        if (!Vars.bottom)
          Vars.anc = ancset->base[j - 1];
      }
      hyptrav(q->back, ancset->base, b1, b2, Vars.bottom,
                treenode, garbage, basechar);
      q = q->next;
    } while (q != Vars.r);
  }
  chuckbase(ancset, garbage);
}  /* hyptrav */


void hypstates(long chars, node *root, pointarray treenode,
                        gbases **garbage, Char *basechar)
{
  /* fill in and describe states at interior nodes */
  /* used in dnacomp, dnapars, & dnapenny */
  long i, n;
  baseptr nothing;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                            ");
  if (dotdiff)
    fprintf(outfile, " ( . means same as in the node below it on tree)\n");
  nothing = (baseptr)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++)
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++) {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    hyptrav(root, nothing, i * 40 - 39, n, true, treenode, garbage, basechar);
  }
  free(nothing);
}  /* hypstates */


void initbranchlen(node *p)
{
  node *q;

  p->v = 0.0;
  if (p->back)
    p->back->v = 0.0;
  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    initbranchlen(q->back);
    q = q->next;
  }
  q = p->next;
  while (q != p) {
    q->v = 0.0;
    q = q->next;
  }
} /* initbranchlen */


void initmin(node *p, long sitei, boolean internal)
{
  long i;

  if (internal) {
    for (i = (long)A; i <= (long)O; i++) {
      p->cumlengths[i] = 0;
      p->numreconst[i] = 1;
    }
  } else {
    for (i = (long)A; i <= (long)O; i++) {
      if (p->base[sitei - 1] & (1 << i)) {
        p->cumlengths[i] = 0;
        p->numreconst[i] = 1;
      } else {
        p->cumlengths[i] = -1;
        p->numreconst[i] = 0;
      }
    }
  }
} /* initmin */


void initbase(node *p, long sitei)
{
  /* traverse tree to initialize base at internal nodes */
  node *q;
  long i, largest;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    if (q->back) {
      memcpy(q->numnuc, p->numnuc, endsite*sizeof(nucarray));
      for (i = (long)A; i <= (long)O; i++) {
        if (q->back->base[sitei - 1] & (1 << i))
          q->numnuc[sitei - 1][i]--;
      }
      if (p->back) {
        for (i = (long)A; i <= (long)O; i++) {
          if (p->back->base[sitei - 1] & (1 << i))
            q->numnuc[sitei - 1][i]++;
        }
      }
      largest = getlargest(q->numnuc[sitei - 1]);
      q->base[sitei - 1] = 0;
      for (i = (long)A; i <= (long)O; i++) {
        if (q->numnuc[sitei - 1][i] == largest)
          q->base[sitei - 1] |= (1 << i);
      }
    }
    q = q->next;
  }
  q = p->next;
  while (q != p) {
    initbase(q->back, sitei);
    q = q->next;
  }
} /* initbase */


void inittreetrav(node *p, long sitei)
{
  /* traverse tree to clear boolean initialized and set up base */
  node *q;

  if (p->tip) {
    initmin(p, sitei, false);
    p->initialized = true;
    return;
  }
  q = p->next;
  while (q != p) {
    inittreetrav(q->back, sitei);
    q = q->next;
  }
  initmin(p, sitei, true);
  p->initialized = false;
  q = p->next;
  while (q != p) {
    initmin(q, sitei, true);
    q->initialized = false;
    q = q->next;
  }
} /* inittreetrav */


void compmin(node *p, node *desc)
{
  /* computes minimum lengths up to p */
  long i, j, minn, cost, desclen, descrecon=0, maxx;

  maxx = 10 * spp;
  for (i = (long)A; i <= (long)O; i++) {
    minn = maxx;
    for (j = (long)A; j <= (long)O; j++) {
      if (transvp) {
        if (
               (
                ((i == (long)A) || (i == (long)G))
             && ((j == (long)A) || (j == (long)G))
               )
            || (
                ((j == (long)C) || (j == (long)T))
             && ((i == (long)C) || (i == (long)T))
               )
           )
          cost = 0;
        else
          cost = 1;
      } else {
        if (i == j)
          cost = 0;
        else
          cost = 1;
      }
      if (desc->cumlengths[j] == -1) {
        desclen = maxx;
      } else {
        desclen = desc->cumlengths[j];
      }
      if (minn > cost + desclen) {
        minn = cost + desclen;
        descrecon = 0;
      }
      if (minn == cost + desclen) {
        descrecon += desc->numreconst[j];
      }
    }
    p->cumlengths[i] += minn;
    p->numreconst[i] *= descrecon;
  }
  p->initialized = true;
} /* compmin */


void minpostorder(node *p, pointarray treenode)
{
  /* traverses an n-ary tree, computing minimum steps at each node */
  node *q;

  if (p->tip) {
    return;
  }
  q = p->next;
  while (q != p) {
    if (q->back)
      minpostorder(q->back, treenode);
    q = q->next;
  }
  if (!p->initialized) {
    q = p->next;
    while (q != p) {
      if (q->back)
        compmin(p, q->back);
      q = q->next;
    }
  }
}  /* minpostorder */


void branchlength(node *subtr1, node *subtr2, double *brlen,
                        pointarray treenode)
{
  /* computes a branch length between two subtrees for a given site */
  long i, j, minn, cost, nom, denom;
  node *temp;

  if (subtr1->tip) {
    temp = subtr1;
    subtr1 = subtr2;
    subtr2 = temp;
  }
  if (subtr1->index == outgrno) {
    temp = subtr1;
    subtr1 = subtr2;
    subtr2 = temp;
  }
  minpostorder(subtr1, treenode);
  minpostorder(subtr2, treenode);
  minn = 10 * spp;
  nom = 0;
  denom = 0;
  for (i = (long)A; i <= (long)O; i++) {
    for (j = (long)A; j <= (long)O; j++) {
      if (transvp) {
        if (
               (
                ((i == (long)A) || (i == (long)G))
             && ((j == (long)A) || (j == (long)G))
               )
            || (
                ((j == (long)C) || (j == (long)T))
             && ((i == (long)C) || (i == (long)T))
               )
           )
          cost = 0;
        else
          cost = 1;
      } else {
        if (i == j)
          cost = 0;
        else
          cost = 1;
      }
      if (subtr1->cumlengths[i] != -1 && (subtr2->cumlengths[j] != -1)) {
        if (subtr1->cumlengths[i] + cost + subtr2->cumlengths[j] < minn) {
          minn = subtr1->cumlengths[i] + cost + subtr2->cumlengths[j];
          nom = 0;
          denom = 0;
        }
        if (subtr1->cumlengths[i] + cost + subtr2->cumlengths[j] == minn) {
          nom += subtr1->numreconst[i] * subtr2->numreconst[j] * cost;
          denom += subtr1->numreconst[i] * subtr2->numreconst[j];
        }
      }
    }
  }
  *brlen = (double)nom/(double)denom;
} /* branchlength */  


void printbranchlengths(node *p)
{
  node *q;
  long i;

  if (p->tip)
    return;
  q = p->next;
  do {
    fprintf(outfile, "%6ld      ",q->index - spp);
    if (q->back->tip) {
      for (i = 0; i < nmlngth; i++)
        putc(nayme[q->back->index - 1][i], outfile);
    } else
      fprintf(outfile, "%6ld    ", q->back->index - spp);
    fprintf(outfile, "   %f\n",q->v);
    if (q->back)
      printbranchlengths(q->back);
    q = q->next;
  } while (q != p);
} /* printbranchlengths */


void branchlentrav(node *p, node *root, long sitei, long chars,
                        double *brlen, pointarray treenode)
  {
  /*  traverses the tree computing tree length at each branch */
  node *q;

  if (p->tip)
    return;
  if (p->index == outgrno)
    p = p->back;
  q = p->next;
  do {
    if (q->back) {
      branchlength(q, q->back, brlen, treenode);
      q->v += ((weight[sitei - 1] / 10.0) * (*brlen)/chars);
      q->back->v += ((weight[sitei - 1] / 10.0) * (*brlen)/chars);
      if (!q->back->tip)
        branchlentrav(q->back, root, sitei, chars, brlen, treenode);
    }
    q = q->next;
  } while (q != p);
}  /* branchlentrav */


void treelength(node *root, long chars, pointarray treenode)
  {
  /*  calls branchlentrav at each site */
  long sitei;
  double trlen;

  initbranchlen(root);
  for (sitei = 1; sitei <= endsite; sitei++) {
    trlen = 0.0;
    initbase(root, sitei);
    inittreetrav(root, sitei);
    branchlentrav(root, root, sitei, chars, &trlen, treenode);
  }
} /* treelength */


void coordinates(node *p, long *tipy, double f, long *fartemp)
{
  /* establishes coordinates of nodes for display without lengths */
  node *q, *first, *last;
  node *mid1 = NULL, *mid2 = NULL;
  long numbranches, numb2;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  numbranches = 0;
  q = p->next;
  do {
    coordinates(q->back, tipy, f, fartemp);
    numbranches += 1;
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p) 
    q = q->next;
  last = q->back;
  numb2 = 1;
  q = p->next;
  while (q != p) {
    if (numb2 == (long)(numbranches + 1)/2)
      mid1 = q->back;
    if (numb2 == (long)(numbranches/2 + 1))
      mid2 = q->back;
    numb2 += 1;
    q = q->next;
  }
  p->xcoord = (long)((double)(last->ymax - first->ymin) * f);
  p->ycoord = (long)((mid1->ycoord + mid2->ycoord) / 2);
  p->ymin = first->ymin;
  p->ymax = last->ymax;
  if (p->xcoord > *fartemp)
    *fartemp = p->xcoord;
}  /* coordinates */


void drawline(long i, double scale, node *root)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done, noplus;

  p = root;
  q = root;
  extra = false;
  noplus = false;
  if (i == (long)p->ycoord && p == root) {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
    noplus = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)(scale * (p->xcoord - q->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if (noplus) {
        putc('-', outfile);
        noplus = false;
      }
      else
        putc('+', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
        noplus = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i
            && i != (long)p->ycoord) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
      noplus = false;
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
      noplus = false;
    }
    if (p != q)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree(node *root, double f)
{
  /* prints out diagram of the tree */
  /* used in dnacomp, dnapars, & dnapenny */
  long i, tipy, dummy;
  double scale;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  dummy = 0;
  coordinates(root, &tipy, f, &dummy);
  scale = 1.5;
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale, root);
  fprintf(outfile, "\n  remember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
}  /* printree */


void writesteps(long chars, boolean weights, steptr oldweight, node *root)
{
  /* used in dnacomp, dnapars, & dnapenny */
  long i, j, k, l;

  putc('\n', outfile);
  if (weights)
    fprintf(outfile, "weighted ");
  fprintf(outfile, "steps in each site:\n");
  fprintf(outfile, "      ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%4ld", i);
  fprintf(outfile, "\n     *------------------------------------");
  fprintf(outfile, "-----\n");
  for (i = 0; i <= (chars / 10); i++) {
    fprintf(outfile, "%5ld", i * 10);
    putc('|', outfile);
    for (j = 0; j <= 9; j++) {
      k = i * 10 + j;
      if (k == 0 || k > chars)
        fprintf(outfile, "    ");
      else {
        l = location[ally[k - 1] - 1];
        if (oldweight[k - 1] > 0)
          fprintf(outfile, "%4ld",
                  oldweight[k - 1] *
                  (root->numsteps[l - 1] / weight[l - 1]));
        else
          fprintf(outfile, "   0");
      }
    }
    putc('\n', outfile);
  }
} /* writesteps */


void treeout(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in dnacomp, dnamove, dnapars, & dnapenny */
  node *q;
  long i, n;
  Char c;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  } else {
    putc('(', outtree);
    (*col)++;
    q = p->next;
    while (q != p) {
      treeout(q->back, nextree, col, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 60) {
        putc('\n', outtree);
        *col = 0;
      }
    }
    putc(')', outtree);
    (*col)++;
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout */


void treeout3(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in dnapars -- writes branch lengths */
  node *q;
  long i, n, w;
  double x;
  Char c;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  } else {
    putc('(', outtree);
    (*col)++;
    q = p->next;
    while (q != p) {
      treeout3(q->back, nextree, col, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 60) {
        putc('\n', outtree);
        *col = 0;
      }
    }
    putc(')', outtree);
    (*col)++;
  }
  x = p->v;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p != root) {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    *col += w + 8; 
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout3 */


/* FIXME curtree should probably be passed by reference */
void drawline2(long i, double scale, tree curtree)
{
  fdrawline2(outfile, i, scale, &curtree);
}

void fdrawline2(FILE *fp, long i, double scale, tree *curtree)
{
  /* draws one row of the tree diagram by moving up tree */
  /* used in dnaml, proml, & restml */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first =NULL, *last =NULL;
  boolean done;

  p = curtree->start;
  q = curtree->start;
  extra = false;
  if (i == (long)p->ycoord && p == curtree->start) {
    if (p->index - spp >= 10)
      fprintf(fp, " %2ld", p->index - spp);
    else
      fprintf(fp, "  %ld", p->index - spp);
    extra = true;
  } else
    fprintf(fp, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (p != curtree->start && r == p) ||
                 (p == curtree->start && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (p == curtree->start)
        last = p->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if ((long)p->ycoord != (long)q->ycoord)
        putc('+', fp);
      else
        putc('-', fp);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', fp);
        if (q->index - spp >= 10)
          fprintf(fp, "%2ld", q->index - spp);
        else
          fprintf(fp, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', fp);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          (i != (long)p->ycoord || p == curtree->start)) {
        putc('|', fp);
        for (j = 1; j < n; j++)
          putc(' ', fp);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', fp);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', fp);
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index-1][j], fp);
  }
  putc('\n', fp);
}  /* drawline2 */


void drawline3(long i, double scale, node *start)
{
  /* draws one row of the tree diagram by moving up tree */
  /* used in dnapars */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first =NULL, *last =NULL;
  boolean done;

  p = start;
  q = start;
  extra = false;
  if (i == (long)p->ycoord) {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (r == p))); 
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if ((long)p->ycoord != (long)q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          (i != (long)p->ycoord || p == start)) {
        putc('|', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index-1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline3 */


void copynode(node *c, node *d, long categs)
{
  long i, j;

  for (i = 0; i < endsite; i++)
    for (j = 0; j < categs; j++)
      memcpy(d->x[i][j], c->x[i][j], sizeof(sitelike));
  memcpy(d->underflows,c->underflows,sizeof(double) * endsite);
  d->tyme = c->tyme;
  d->v = c->v;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
  d->iter = c->iter;                   /* iter used in dnaml only */
  d->haslength = c->haslength;         /* haslength used in dnamlk only */
  d->initialized = c->initialized;     /* initialized used in dnamlk only */
}  /* copynode */


void prot_copynode(node *c, node *d, long categs)
{
  /* a version of copynode for proml */
  long i, j;

  for (i = 0; i < endsite; i++)
    for (j = 0; j < categs; j++)
      memcpy(d->protx[i][j], c->protx[i][j], sizeof(psitelike));
  memcpy(d->underflows,c->underflows,sizeof(double) * endsite);
  d->tyme = c->tyme;
  d->v = c->v;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
  d->iter = c->iter;                   /* iter used in dnaml only */
  d->haslength = c->haslength;         /* haslength used in dnamlk only */
  d->initialized = c->initialized;     /* initialized used in dnamlk only */
}  /* prot_copynode */


void copy_(tree *a, tree *b, long nonodes, long categs)
{
  /* used in dnamlk */
  long i;
  node *p, *q, *r, *s, *t;

  for (i = 0; i < spp; i++) {
    copynode(a->nodep[i], b->nodep[i], categs);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes; i++) {
    if (a->nodep[i]) {
      p = a->nodep[i];
      q = b->nodep[i];
          r = p;
      do {
        copynode(p, q, categs);
        if (p->back) {
          s = a->nodep[p->back->index - 1];
          t = b->nodep[p->back->index - 1];
          if (s->tip) {
            if(p->back == s)
              q->back = t;
          } else {
            do {
              if (p->back == s)
                q->back = t;
              s = s->next;
              t = t->next;
            } while (s != a->nodep[p->back->index - 1]);
          }
        }
        else
          q->back = NULL;
        p = p->next;
        q = q->next;
      } while (p != r);
    }
  }
  b->likelihood = a->likelihood;
  b->start = a->start;               /* start used in dnaml only */
  b->root = a->root;                 /* root used in dnamlk only */
}  /* copy_ */


void prot_copy_(tree *a, tree *b, long nonodes, long categs)
{
  /* used in promlk */
  /* identical to copy_() except for calls to prot_copynode rather */
  /* than copynode.                                                */
  long i;
  node *p, *q, *r, *s, *t;

  for (i = 0; i < spp; i++) {
    prot_copynode(a->nodep[i], b->nodep[i], categs);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1]->next
) 
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes; i++) {
    if (a->nodep[i]) {
      p = a->nodep[i];
      q = b->nodep[i];
          r = p;
      do {
        prot_copynode(p, q, categs);
        if (p->back) {
          s = a->nodep[p->back->index - 1];
          t = b->nodep[p->back->index - 1];
          if (s->tip)
            {
                if(p->back == s)
                  q->back = t;
          } else {
            do {
              if (p->back == s)
                q->back = t;
              s = s->next;
              t = t->next;
            } while (s != a->nodep[p->back->index - 1]);
          }
        }
        else
          q->back = NULL;
        p = p->next;
        q = q->next;
      } while (p != r);
    }
  }
  b->likelihood = a->likelihood;
  b->start = a->start;               /* start used in dnaml only */
  b->root = a->root;                 /* root used in dnamlk only */
}  /* prot_copy_ */


void standev(long chars, long numtrees, long minwhich, double minsteps,
                        double *nsteps, long **fsteps, longer seed)
{  /* do paired sites test (KHT or SH test) on user-defined trees */
   /* used in dnapars & protpars */
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;
  double **covar, *P, *f, *r;

#define SAMPLES 1000
  if (numtrees == 2) {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees) {
      fprintf(outfile, "%3ld%10.1f", which, nsteps[which - 1] / 10);
      if (minwhich == which)
        fprintf(outfile, "  <------ best\n");
      else {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = 0; i < endsite; i++) {
          if (weight[i] > 0) {
            wt = weight[i] / 10.0;
            sumw += wt;
            temp = (fsteps[which - 1][i] - fsteps[minwhich - 1][i]) / 10.0;
            sum += wt * temp;
            sum2 += wt * temp * temp;
          }
        }
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
        fprintf(outfile, "%10.1f%12.4f",
                (nsteps[which - 1] - minsteps) / 10, sd);
        if ((sum > 0.0) && (sum > 1.95996 * sd))
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  } else {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES){
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n"
              , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    } else {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees*sizeof(double *));  
    sumw = 0.0;
    for (i = 0; i < endsite; i++)
      sumw += weight[i] / 10.0;
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = nsteps[i]/(10.0*sumw);
      for (j = 0; j <=i; j++) {
        sum2 = nsteps[j]/(10.0*sumw);
        temp = 0.0;
        for (k = 0; k < endsite; k++) {
          if (weight[k] > 0) {
            wt = weight[k]/10.0;
            temp = temp + wt*(fsteps[i][k]/10.0-sum)
                            *(fsteps[j][k]/10.0-sum2);
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++) { /* in-place Cholesky decomposition
                                        of trees x trees covariance matrix */
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++) {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resampled sums */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0]/10.0;               /* sum2 will be smallest # of steps */
    for (i = 1; i < numtrees; i++)
      if (sum2 > nsteps[i]/10.0)
        sum2 = nsteps[i]/10.0;
    for (i = 1; i <= SAMPLES; i++) {          /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++) {        /* compute vectors */
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get min of vector */
        if (f[j] < sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (nsteps[j]/10.0-sum2 <= f[j] - sum)
          P[j] += 1.0/SAMPLES;
    }
    fprintf(outfile, "Tree    Steps   Diff Steps   P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++) {
      fprintf(outfile, "%3ld%10.1f", i+1, nsteps[i]/10);
      if ((minwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else {
        fprintf(outfile, " %9.1f  %10.3f", nsteps[i]/10.0-sum2, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
  fprintf(outfile, "\n");
  free(P);             /* free the variables we Malloc'ed */
  free(f);
  free(r);
  for (i = 0; i < numtrees; i++)
    free(covar[i]);
  free(covar);
  }
}  /* standev */


void standev2(long numtrees, long maxwhich, long a, long b, double maxlogl,
              double *l0gl, double **l0gf, steptr aliasweight, longer seed)
{  /* do paired sites test (KHT or SH) for user-defined trees */
  /* used in dnaml, dnamlk, proml, promlk, and restml */
  double **covar, *P, *f, *r;
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;

#define SAMPLES 1000
  if (numtrees == 2) {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    logL    Diff logL    Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees) {
      fprintf(outfile, "%3ld %9.1f", which, l0gl[which - 1]);
      if (maxwhich == which)
        fprintf(outfile, "  <------ best\n");
      else {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = a; i <= b; i++) {
          if (aliasweight[i] > 0) {
            wt = aliasweight[i];
            sumw += wt;
            temp = l0gf[which - 1][i] - l0gf[maxwhich - 1][i];
            sum += temp * wt;
            sum2 += wt * temp * temp;
          }
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw ));
        fprintf(outfile, "%10.1f %11.4f", (l0gl[which - 1])-maxlogl, sd);
        if ((sum < 0.0) && ((-sum) > 1.95996 * sd))
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  } else {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES){
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n"
              , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    } else {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees*sizeof(double *));  
    sumw = 0.0;
    for (i = a; i <= b; i++)
      sumw += aliasweight[i];
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = l0gl[i]/sumw;
      for (j = 0; j <=i; j++) {
        sum2 = l0gl[j]/sumw;
        temp = 0.0;
        for (k = a; k <= b ; k++) {
          if (aliasweight[k] > 0) {
            wt = aliasweight[k];
            temp = temp + wt*(l0gf[i][k]-sum)
                            *(l0gf[j][k]-sum2);
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++) { /* in-place Cholesky decomposition
                                        of trees x trees covariance matrix */
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++) {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resampled likelihoods */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    for (i = 1; i <= SAMPLES; i++) {          /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++) {        /* compute vectors */
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get max of vector */
        if (f[j] > sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (maxlogl-l0gl[j] <= sum-f[j])
          P[j] += 1.0/SAMPLES;
    }
    fprintf(outfile, "Tree    logL    Diff logL    P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++) {
      fprintf(outfile, "%3ld%10.1f", i+1, l0gl[i]);
      if ((maxwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else {
        fprintf(outfile, " %9.1f  %10.3f", l0gl[i]-maxlogl, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
  fprintf(outfile, "\n");
  free(P);             /* free the variables we Malloc'ed */
  free(f);
  free(r);
  for (i = 0; i < numtrees; i++)
    free(covar[i]);
  free(covar);
  }
}  /* standev */


void freetip(node *anode)
{
  /* used in dnacomp, dnapars, & dnapenny */

  free(anode->numsteps);
  free(anode->oldnumsteps);
  free(anode->base);
  free(anode->oldbase);
}  /* freetip */


void freenontip(node *anode)
{
  /* used in dnacomp, dnapars, & dnapenny */

  free(anode->numsteps);
  free(anode->oldnumsteps);
  free(anode->base);
  free(anode->oldbase);
  free(anode->numnuc);
}  /* freenontip */


void freenodes(long nonodes, pointarray treenode)
{
  /* used in dnacomp, dnapars, & dnapenny */
  long i;
  node *p;

  for (i = 0; i < spp; i++)
    freetip(treenode[i]);
  for (i = spp; i < nonodes; i++) {
    if (treenode[i] != NULL) {
      p = treenode[i]->next;
      do {
        freenontip(p);
        p = p->next;
      } while (p != treenode[i]);
      freenontip(p);
    }
  }
}  /* freenodes */


void freenode(node **anode)
{
  /* used in dnacomp, dnapars, & dnapenny */

  freenontip(*anode);
  free(*anode);
}  /* freenode */


void freetree(long nonodes, pointarray treenode)
{
  /* used in dnacomp, dnapars, & dnapenny */
  long i;
  node *p, *q;

  for (i = 0; i < spp; i++)
    free(treenode[i]);
  for (i = spp; i < nonodes; i++) {
    if (treenode[i] != NULL) {
      p = treenode[i]->next;
      do {
        q = p->next;
        free(p);
        p = q;
      } while (p != treenode[i]);
      free(p);
    }
  }
  free(treenode);
}  /* freetree */


void prot_freex_notip(long nonodes, pointarray treenode)
{
  /* used in proml */
  long i, j;
  node *p;

  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    if ( p == NULL ) continue;
    do {
      for (j = 0; j < endsite; j++){
        free(p->protx[j]);
        p->protx[j] = NULL;
      }
      free(p->underflows);
      p->underflows = NULL;
      free(p->protx);
      p->protx = NULL;
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* prot_freex_notip */


void prot_freex(long nonodes, pointarray treenode)
{
  /* used in proml */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++) {
    for (j = 0; j < endsite; j++)
      free(treenode[i]->protx[j]);
    free(treenode[i]->protx);
    free(treenode[i]->underflows);
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    do {
      for (j = 0; j < endsite; j++)
        free(p->protx[j]);
      free(p->protx);
      free(p->underflows);
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* prot_freex */


void freex_notip(long nonodes, pointarray treenode)
{
  /* used in dnaml & dnamlk */
  long i, j;
  node *p;

  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    if ( p == NULL ) continue;
    do {
      for (j = 0; j < endsite; j++)
        free(p->x[j]);
      free(p->underflows);
      free(p->x);
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* freex_notip */


void freex(long nonodes, pointarray treenode)
{
  /* used in dnaml & dnamlk */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++) {
    for (j = 0; j < endsite; j++)
      free(treenode[i]->x[j]);
    free(treenode[i]->x);
    free(treenode[i]->underflows);
  }
  for (i = spp; i < nonodes; i++) {
    if(treenode[i]){
      p = treenode[i];
      do {
        for (j = 0; j < endsite; j++)
          free(p->x[j]);
        free(p->x);
        free(p->underflows);
        p = p->next;
      } while (p != treenode[i]);
    }
  }
}  /* freex */



void freegarbage(gbases **garbage)
{
  /* used in dnacomp, dnapars, & dnapenny */
  gbases *p;

  while (*garbage) {
    p = *garbage;
    *garbage = (*garbage)->next;
    free(p->base);
    free(p);
  }
}  /*freegarbage */


void freegrbg(node **grbg)
{
  /* used in dnacomp, dnapars, & dnapenny */
  node *p;

  while (*grbg) {
    p = *grbg;
    *grbg = (*grbg)->next;
    freenontip(p);
    free(p);
  }
} /*freegrbg */


void collapsetree(node *p, node *root, node **grbg, pointarray treenode, 
                  long *zeros)
{
  /*  Recurse through tree searching for zero length brances between */
  /*  nodes (not to tips).  If one exists, collapse the nodes together, */
  /*  removing the branch. */
  node *q, *x1, *y1, *x2, *y2;
  long i, index, index2, numd;
  if (p->tip)
    return;
  q = p->next;
  do {
    if (!q->back->tip && q->v == 0.000000) {
      /* merge the two nodes. */
      x1 = y2 = q->next;
      x2 = y1 = q->back->next;
      while(x1->next != q)
        x1 = x1-> next;
      while(y1->next != q->back)
        y1 = y1-> next;
      x1->next = x2;
      y1->next = y2;

      index = q->index;
      index2 = q->back->index;
      numd = treenode[index-1]->numdesc + q->back->numdesc -1;
      chuck(grbg, q->back);
      chuck(grbg, q);
      q = x2;

      /* update the indicies around the node circle */
      do{
        if(q->index != index){
          q->index = index;
        }
        q = q-> next;
      }while(x2 != q);
      updatenumdesc(treenode[index-1], root, numd);
       
      /* Alter treenode to point to real nodes, and update indicies */
      /* acordingly. */
      i=0;
      for(i = (index2-1); i < nonodes-1 && treenode[i+1]; i++){ 
        treenode[i]=treenode[i+1];
        treenode[i+1] = NULL;
        x1=x2=treenode[i]; 
        do{ 
          x1->index = i+1; 
          x1 = x1 -> next; 
        } while(x1 != x2); 
      }

      /* Create a new empty fork in the blank spot of treenode */
      x1=NULL;
      for(i=1; i <=3 ; i++){
        gnutreenode(grbg, &x2, index2, endsite, zeros);
        x2->next = x1;
        x1 = x2;
      }
      x2->next->next->next = x2;
      treenode[nonodes-1]=x2;
      if (q->back)
        collapsetree(q->back, root, grbg, treenode, zeros);
    } else {
      if (q->back)
        collapsetree(q->back, root, grbg, treenode, zeros);
      q = q->next;
    }
  } while (q != p);
} /* collapsetree */


void collapsebestrees(node **root, node **grbg, pointarray treenode, 
                      bestelm *bestrees, long *place, long *zeros,
                      long chars, boolean recompute, boolean progress)
     
{
  /* Goes through all best trees, collapsing trees where possible, and  */
  /* deleting trees that are not unique.    */
  long i,j, k, pos, nextnode, oldnextree;
  boolean found;
  node *dummy;

  oldnextree = nextree;
  for(i = 0 ; i < (oldnextree - 1) ; i++){
    bestrees[i].collapse = true;
  }

  if(progress)
    printf("Collapsing best trees\n   ");
  k = 0;
  for(i = 0 ; i < (oldnextree - 1) ; i++){
    if(progress){
      if(i % (((oldnextree-1) / 72) + 1) == 0)
        putchar('.');
      fflush(stdout);
    }
    while(!bestrees[k].collapse)
      k++;
    /* Reconstruct tree. */
    *root = treenode[0];
    add(treenode[0], treenode[1], treenode[spp], root, recompute,
        treenode, grbg, zeros);
    nextnode = spp + 2;
    for (j = 3; j <= spp; j++) {
      if (bestrees[k].btree[j - 1] > 0)
        add(treenode[bestrees[k].btree[j - 1] - 1], treenode[j - 1],
            treenode[nextnode++ - 1], root, recompute, treenode, grbg,
            zeros);
      else
          add(treenode[treenode[-bestrees[k].btree[j - 1]-1]->back->index-1],
              treenode[j - 1], NULL, root, recompute, treenode, grbg, zeros);
    }
    reroot(treenode[outgrno - 1], *root);

    treelength(*root, chars, treenode);
    collapsetree(*root, *root, grbg, treenode, zeros);
    savetree(*root, place, treenode, grbg, zeros);
    /* move everything down in the bestree list */
    for(j = k ; j < (nextree - 2) ; j++){
      memcpy(bestrees[j].btree, bestrees[j + 1].btree, spp * sizeof(long));
      bestrees[j].gloreange = bestrees[j + 1].gloreange;
      bestrees[j + 1].gloreange = false;
      bestrees[j].locreange = bestrees[j + 1].locreange;
      bestrees[j + 1].locreange = false;
      bestrees[j].collapse = bestrees[j + 1].collapse;
    }
    pos=0;
    findtree(&found, &pos, nextree-1, place, bestrees);    

    /* put the new tree at the end of the list if it wasn't found */
    nextree--;
    if(!found)
      addtree(pos, &nextree, false, place, bestrees);

    /* Deconstruct the tree */
    for (j = 1; j < spp; j++){
      re_move(treenode[j], &dummy, root, recompute, treenode,
              grbg, zeros);
    }
  }
  if (progress) {
    putchar('\n');
#ifdef WIN32
#endif
  }
}
