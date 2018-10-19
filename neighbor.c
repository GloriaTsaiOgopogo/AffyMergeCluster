
#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
#endif /* OSX_CARBON */

#include <stdio.h>
#include <signal.h>
#include <float.h>

#include "neighbor.h"
//#include "R.h"

#ifdef WIN32
#include <windows.h>
/* for console code (clear screen, text color settings) */
CONSOLE_SCREEN_BUFFER_INFO      savecsbi;
boolean savecsbi_valid = false;
HANDLE  hConsoleOutput;

#endif /* WIN32 */

#ifndef OLDC
/* function prototypes */
void getoptions(void);
void allocrest(void);
void doinit(void);
void inputoptions(void);
void getinput(void);
void describe(node *, double);
void summarize(void);
void jointree(void);
void maketree(void);
void freerest(void);

void alloctree(pointptr *treenode, long nonodes);
void freetree(pointptr *treenode, long nonodes);
void allocd(long nonodes, pointptr treenode);
void freed(long nonodes, pointptr treenode);
void allocw(long nonodes, pointptr treenode);
void freew(long nonodes, pointptr treenode);
void setuptree(tree *a, long nonodes);
void inputdata(boolean replicates, boolean printdata, boolean lower,
                        boolean upper, vector *x, intvector *reps);
void coordinates(node *p, double lengthsum, long *tipy, double *tipmax,
                        node *start, boolean njoin);
void drawline(long i, double scale, node *start, boolean rooted);
void printree(node *start, boolean treeprint,
                        boolean        njoin, boolean rooted);
void treeoutr(node *p, long *col, tree *curtree);
void treeout(node *p, long *col, double m, boolean njoin, node *start);

static void crash_handler(int signum);
/* function prototypes */
#endif

#if defined(OSX_CARBON) && defined(__MWERKS__)
boolean fixedpath = false;
#endif

Char infilename[FNMLNGTH], outtreename[FNMLNGTH];
long nonodes2, outgrno, col, datasets, ith, inseed, spp, words, bits;
vector *x;
intvector *reps;
tree curtree;
longer seed;
long *enterorder;
Char progname[20];
boolean ibmpc, ansi, tranvsp;
naym *nayme;                     /* names of species */

/* variables for maketree, propagated globally for C version: */
node **cluster;
FILE *outfile, *infile, *intree, *intree2, *outtree, *weightfile, *catfile, *ancfile, *mixfile, *factfile;

void getoptions()
{
  /* interactively set options */
  long loopcount;
    
  outgrno = 1;
  loopcount = 0;
  for(;;) {
    //cleerhome();
    printf("\nNeighbor-Joining method revised from PHYLIP, version %s\n\n",VERSION);
    fflush(stdout);
    break;    
    countup(&loopcount, 100);
  }
}  /* getoptions */


void allocrest()
{
  long i;
  //vector *x; //9.17
  //intvector *reps; //9.17  
  
  x = (vector *)Malloc(spp*sizeof(vector));
  for (i = 0; i < spp; i++)
    x[i] = (vector)Malloc(spp*sizeof(double));
  reps = (intvector *)Malloc(spp*sizeof(intvector));
  for (i = 0; i < spp; i++)
    reps[i] = (intvector)Malloc(spp*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
  cluster = (node **)Malloc(spp*sizeof(node *));
}  /* allocrest */


void freerest()
{
  long i;

  for (i = 0; i < spp; i++)
    free(x[i]);
  free(x);
  for (i = 0; i < spp; i++)
    free(reps[i]);
  free(reps);
  free(nayme);
  free(enterorder);
  free(cluster);
}  /* freerest */


void doinit()
{
  /* initializes variables */
  node *p;

  inputnumbers2(&spp, &nonodes2, 2);
  nonodes2 += 0;
  getoptions();
  alloctree(&curtree.nodep, nonodes2+1);
  p = curtree.nodep[nonodes2]->next;
  curtree.nodep[nonodes2]->next = curtree.nodep[nonodes2];
  free(p->next);
  free(p);
  allocrest();
  
}  /* doinit */


void inputoptions()
{
  /* read options information */

  if (ith != 1)
    samenumsp2(ith);
}  /* inputoptions */


void describe(node *p, double height)
{
  /* print out information for one branch */
  node *q;

  q = p->back;
  if (!p->tip) {
    describe(p->next->back, height+q->v);
    describe(p->next->next->back, height+q->v);
  }
}  /* describe */


void summarize()
{
  /* print out branch lengths etc. */
  describe(curtree.start->next->back, 0.0);
  describe(curtree.start->next->next->back, 0.0);
  describe(curtree.start->back, 0.0);
}  /* summarize */


void jointree()
{
  /* calculate the tree */
  long nc, nextnode, mini=0, minj=0, i, j, ia, ja, ii, jj, nude, iter;
  double fotu2, total, tmin, dio, djo, bi, bj, bk, dmin=0, da;
  long el[3];
  vector av;
  intvector oc;

  double *R;   /* added in revisions by Y. Ina */
  R = (double *)Malloc(spp * sizeof(double));

  for (i = 0; i <= spp - 2; i++) {
    for (j = i + 1; j < spp; j++) {
      da = (x[i][j] + x[j][i]) / 2.0;
      x[i][j] = da;
      x[j][i] = da;
    }
  }
  /* First initialization */
  fotu2 = spp - 2.0;
  nextnode = spp + 1;
  av = (vector)Malloc(spp*sizeof(double));
  oc = (intvector)Malloc(spp*sizeof(long));
  for (i = 0; i < spp; i++) {
    av[i] = 0.0;
    oc[i] = 1;
  }
  /* Enter the main cycle */
  iter = spp - 3;
  for (nc = 1; nc <= iter; nc++) {
    for (j = 2; j <= spp; j++) {
      for (i = 0; i <= j - 2; i++)
        x[j - 1][i] = x[i][j - 1];
    }
    tmin = DBL_MAX;
    /* Compute sij and minimize */
        for (i = 0; i < spp; i++)
      R[i] = 0.0;
    for (ja = 2; ja <= spp; ja++) {
      jj = enterorder[ja - 1];
      if (cluster[jj - 1] != NULL) {
        for (ia = 0; ia <= ja - 2; ia++) {
          ii = enterorder[ia];
          if (cluster[ii - 1] != NULL) {
            R[ii - 1] += x[ii - 1][jj - 1];
            R[jj - 1] += x[ii - 1][jj - 1];
          }
        }
      }
    }
	
    for (ja = 2; ja <= spp; ja++) {
      jj = enterorder[ja - 1];
      if (cluster[jj - 1] != NULL) {
        for (ia = 0; ia <= ja - 2; ia++) {
          ii = enterorder[ia];
          if (cluster[ii - 1] != NULL) {
            total = fotu2 * x[ii - 1][jj - 1] - R[ii - 1] - R[jj - 1];            
            if (total < tmin) {
              tmin = total;
              mini = ii;
              minj = jj;
            }
          }
        }
      }
    }
    /* compute lengths and print */
    
	dio = 0.0;
    djo = 0.0;
    for (i = 0; i < spp; i++) {
      dio += x[i][mini - 1];
      djo += x[i][minj - 1];
    }
    dmin = x[mini - 1][minj - 1];
    dio = (dio - dmin) / fotu2;
    djo = (djo - dmin) / fotu2;
    bi = (dmin + dio - djo) * 0.5;
    bj = dmin - bi;
    bi -= av[mini - 1];
    bj -= av[minj - 1];
	
    hookup(curtree.nodep[nextnode - 1]->next, cluster[mini - 1]);
    hookup(curtree.nodep[nextnode - 1]->next->next, cluster[minj - 1]);
    cluster[mini - 1]->v = bi;
    cluster[minj - 1]->v = bj;
    cluster[mini - 1]->back->v = bi;
    cluster[minj - 1]->back->v = bj;
    cluster[mini - 1] = curtree.nodep[nextnode - 1];
    cluster[minj - 1] = NULL;
    nextnode++;
    av[mini - 1] = dmin * 0.5;
    /* re-initialization */
    fotu2 -= 1.0;
    for (j = 0; j < spp; j++) {
      if (cluster[j] != NULL) {
        da = (x[mini - 1][j] + x[minj - 1][j]) * 0.5;
        if (mini - j - 1 < 0)
          x[mini - 1][j] = da;
        if (mini - j - 1 > 0)
          x[j][mini - 1] = da;
      }
    }
    for (j = 0; j < spp; j++) {
      x[minj - 1][j] = 0.0;
      x[j][minj - 1] = 0.0;
    }
    oc[mini - 1] += oc[minj - 1];
  }
  /* the last cycle */
  nude = 1;
  for (i = 1; i <= spp; i++) {
    if (cluster[i - 1] != NULL) {
      el[nude - 1] = i;
      nude++;
    }
  }

  bi = (x[el[0] - 1][el[1] - 1] + x[el[0] - 1][el[2] - 1] - x[el[1] - 1]
        [el[2] - 1]) * 0.5;
  bj = x[el[0] - 1][el[1] - 1] - bi;
  bk = x[el[0] - 1][el[2] - 1] - bi;
  bi -= av[el[0] - 1];
  bj -= av[el[1] - 1];
  bk -= av[el[2] - 1];

  hookup(curtree.nodep[nextnode - 1], cluster[el[0] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next, cluster[el[1] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next->next, cluster[el[2] - 1]);
  cluster[el[0] - 1]->v = bi;
  cluster[el[1] - 1]->v = bj;
  cluster[el[2] - 1]->v = bk;
  cluster[el[0] - 1]->back->v = bi;
  cluster[el[1] - 1]->back->v = bj;
  cluster[el[2] - 1]->back->v = bk;
  curtree.start = cluster[el[0] - 1]->back;
  free(av);
  free(oc);
  free(R);
}  /* jointree */


void maketree()
{
  /* construct the tree */
  long i ;
  inputdata(false, false, true, false, x, reps);
  if (spp < 3) {
    printf("\nERROR: Neighbor-Joining runs must have at least 3 species\n\n");
    exxit(-1);
  }
  if (ith == 1)
    setuptree(&curtree, nonodes2 + 1);
  for (i = 1; i <= spp; i++)
    enterorder[i - 1] = i;
  randumize(seed, enterorder);
  for (i = 0; i < spp; i++)
    cluster[i] = curtree.nodep[i];
  jointree();
  curtree.start = curtree.nodep[outgrno - 1]->back;
  printree(curtree.start, true, true, false);
  summarize();

  col = 0;
  treeout(curtree.start, &col, 0.43429448222, true, curtree.start);  
}  /* maketree */

///////////////////////  functions from dist.c //////////////

void alloctree(pointptr *treenode, long nonodes)
{
  /* allocate spp tips and (nonodes - spp) forks, each containing three
   * nodes. Fill in treenode where 0..spp-1 are pointers to tip nodes, and
   * spp..nonodes-1 are pointers to one node in each fork. */

  /* used in fitch, kitsch, neighbor */

  long i, j;
  node *p, *q;

  *treenode = (pointptr)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < spp; i++)
    (*treenode)[i] = (node *)Malloc(sizeof(node));
  for (i = spp; i < nonodes; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    (*treenode)[i] = p;
  }
} /* alloctree */


void freetree(pointptr *treenode, long nonodes)
{
  long i;
  node *p, *q;

  for (i = 0; i < spp; i++)
    free((*treenode)[i]);
  for (i = spp; i < nonodes; i++) {
    p = (*treenode)[i];
    q = p->next;
    while(q != p)
    {
        node * r = q;
        q = q->next;
        free(r);
    }
    free(p);
  }
  free(*treenode);
} /* freetree */


void allocd(long nonodes, pointptr treenode)
{
  /* used in fitch & kitsch */
  long i, j;
  node *p;

  for (i = 0; i < (spp); i++) {
    treenode[i]->d = (vector)Malloc(nonodes*sizeof(double));
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    for (j = 1; j <= 3; j++) {
      p->d = (vector)Malloc(nonodes*sizeof(double));
      p = p->next;
    }
  }
}


void freed(long nonodes, pointptr treenode)
{
  /* used in fitch */
  long i, j;
  node *p;

  for (i = 0; i < (spp); i++) {
    free(treenode[i]->d);
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    for (j = 1; j <= 3; j++) {
      free(p->d);
      p = p->next;
    }
  }
}


void allocw(long nonodes, pointptr treenode)
{
  /* used in fitch & kitsch */
  long i, j;
  node *p;

  for (i = 0; i < (spp); i++) {
      treenode[i]->w = (vector)Malloc(nonodes*sizeof(double));
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    for (j = 1; j <= 3; j++) {
      p->w = (vector)Malloc(nonodes*sizeof(double));
      p = p->next;
    }
  }
}


void freew(long nonodes, pointptr treenode)
{
  /* used in fitch */
  long i, j;
  node *p;

  for (i = 0; i < (spp); i++) {
    free(treenode[i]->w);
  }
  for (i = spp; i < nonodes; i++) {
    p = treenode[i];
    for (j = 1; j <= 3; j++) {
      free(p->w);
      p = p->next;
    }
  }
}


void setuptree(tree *a, long nonodes)
{
  /* initialize a tree */
  /* used in fitch, kitsch, & neighbor */
  long i=0;
  node *p;

  for (i = 1; i <= nonodes; i++) {
    a->nodep[i - 1]->back = NULL;
    a->nodep[i - 1]->tip = (i <= spp);
    a->nodep[i - 1]->iter = true;
    a->nodep[i - 1]->index = i;
    a->nodep[i - 1]->t = 0.0;
    a->nodep[i - 1]->sametime = false;
    a->nodep[i - 1]->v = 0.0;
    if (i > spp) {
      p = a->nodep[i-1]->next;
      while (p != a->nodep[i-1]) {
        p->back = NULL;
        p->tip = false;
        p->iter = true;
        p->index = i;
        p->t = 0.0;
        p->sametime = false;
        p = p->next;
      }
    }
  }
  a->likelihood = -1.0;
  a->start = a->nodep[0];
  a->root = NULL;
}  /* setuptree */


void inputdata(boolean replicates, boolean printdata, boolean lower,
                        boolean upper, vector *x, intvector *reps)
{
  /* read in distance matrix */
  /* used in fitch & neighbor */
  long i=0, j=0, k=0, columns=0;
  boolean skipit=false, skipother=false;

  if (replicates)
    columns = 4;
  else
    columns = 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  for (i = 0; i < spp; i++) {
    x[i][i] = 0.0;
    scan_eoln(infile);
    initname(i);
    for (j = 0; j < spp; j++) {
      skipit = ((lower && j + 1 >= i + 1) || (upper && j + 1 <= i + 1));
      skipother = ((lower && i + 1 >= j + 1) || (upper && i + 1 <= j + 1));
      if (!skipit) {
        if (eoln(infile))
          scan_eoln(infile);
        if (fscanf(infile, "%lf", &x[i][j]) != 1)  {
          printf("The infile is of the wrong type\n");
          exxit(-1);
        }
        if (replicates) {
          if (eoln(infile))
            scan_eoln(infile);
          if (fscanf(infile, "%ld", &reps[i][j]) != 1) {
            printf("The infile is of the wrong type\n");
            exxit(-1);
          }
        } else
          reps[i][j] = 1;
      }
      if (!skipit && skipother) {
          x[j][i] = x[i][j];
          reps[j][i] = reps[i][j];
      }
      if ((i == j) && (fabs(x[i][j]) > 0.000000001)) {
       printf("\nERROR: diagonal element of row %ld of distance matrix ", i+1);
        printf("is not zero.\n");
        printf("       Is it a distance matrix?\n\n");
        exxit(-1);        
      }
      if ((j < i) && (fabs(x[i][j]-x[j][i]) > 0.000000001)) {
        printf("ERROR: distance matrix is not symmetric:\n");
        printf("       (%ld,%ld) element and (%ld,%ld) element are unequal.\n",
               i+1, j+1, j+1, i+1);
        printf("       They are %10.6f and %10.6f, respectively.\n",
               x[i][j], x[j][i]);
        printf("       Is it a distance matrix?\n\n");
        exxit(-1);
      }
    }
  }
  scan_eoln(infile);
  if (!printdata)
    return;
  for (i = 0; i < spp; i++) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[i][j], outfile);
    putc(' ', outfile);
    for (j = 1; j <= spp; j++) {
      fprintf(outfile, "%10.5f", x[i][j - 1]);
      if (replicates)
        fprintf(outfile, " (%3ld)", reps[i][j - 1]);
      if (j % columns == 0 && j < spp) {
        putc('\n', outfile);
        for (k = 1; k <= nmlngth + 1; k++)
          putc(' ', outfile);
      }
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* inputdata */


void coordinates(node *p, double lengthsum, long *tipy, double *tipmax,
                        node *start, boolean njoin)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = (long)(60 * lengthsum + 0.5);
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    if (lengthsum > *tipmax)
      *tipmax = lengthsum;
    return;
  }
  q = p->next;
  do {
    if (q->back)
      coordinates(q->back, lengthsum + q->v, tipy,tipmax, start, njoin);
    q = q->next;
  } while ((p == start || p != q) && (p != start || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p && q->next->back)  /* is this right ? */
    q = q->next;
  last = q->back;
  p->xcoord = (long)(60 * lengthsum + 0.5);
  if (p == start && p->back) 
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void drawline(long i, double scale, node *start, boolean rooted)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n=0, j=0;
  boolean extra=false, trif=false;
  node *r, *first =NULL, *last =NULL;
  boolean done=false;

  p = start;
  q = start;
  extra = false;
  trif = false;
  if (i == (long)p->ycoord && p == start) {  /* display the root */
    if (rooted) {
      if (p->index - spp >= 10)
        fprintf(outfile, "-");
      else
        fprintf(outfile, "--");
    }
    else {
      if (p->index - spp >= 10)
        fprintf(outfile, " ");
      else
        fprintf(outfile, "  ");
    }
    if (p->index - spp >= 10)
      fprintf(outfile, "%2ld", p->index - spp);
    else
      fprintf(outfile, "%ld", p->index - spp);
    extra = true;
    trif = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) { /* internal nodes */
      r = p->next;
      /* r->back here is going to the same node. */
      do {
        if (!r->back) {
          r = r->next;
          continue;
        }
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          break;
        }
        r = r->next;
      } while (!((p != start && r == p) || (p == start && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p) 
        r = r->next;
      last = r->back;
      if (!rooted && (p == start))
        last = p->back;
    } /* end internal node case... */
    /* draw the line: */
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (!q->tip) {
      if ((n < 3) && (q->index - spp >= 10))
        n = 3;
      if ((n < 2) && (q->index - spp < 10))
        n = 2;
    }
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      if (trif) {
        n++;
        trif = false;
      } 
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
      if ((long)last->ycoord > i && (long)first->ycoord < i
           && i != (long)p->ycoord) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
        trif = false;
      }
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree(node *start, boolean treeprint,
                        boolean        njoin, boolean rooted)
{
  /* prints out diagram of the tree */
  /* used in fitch & neighbor */
  long i;
  long tipy;
  double scale,tipmax;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(start, 0.0, &tipy, &tipmax, start, njoin);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale, start, rooted);
  putc('\n', outfile);
}  /* printree */


void treeoutr(node *p, long *col, tree *curtree)
{
  /* write out file with representation of final tree. 
   * Rooted case. Used in kitsch and neighbor. */
  long i, n, w;
  Char c;
  double x;

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
    (*col) += n;
  } else {
    putc('(', outtree);
    (*col)++;
    treeoutr(p->next->back,col,curtree);
    putc(',', outtree);
    (*col)++;
    if ((*col) > 55) {
      putc('\n', outtree);
      (*col) = 0;
    }
    treeoutr(p->next->next->back,col,curtree);
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
  if (p == curtree->root)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    (*col) += w + 8;
  }
}  /* treeoutr */


void treeout(node *p, long *col, double m, boolean njoin, node *start)
{
  /* write out file with representation of final tree */
  /* used in fitch & neighbor */
  long i=0, n=0, w=0;
  Char c;
  double x=0.0;

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
    treeout(p->next->back, col, m, njoin, start);
    putc(',', outtree);
    (*col)++;
    if (*col > 55) {
      putc('\n', outtree);
      *col = 0;
    }
    treeout(p->next->next->back, col, m, njoin, start);
    if (p == start && njoin) {
      putc(',', outtree);
      treeout(p->back, col, m, njoin, start);
    }
    putc(')', outtree);
    (*col)++;
  }
  x = p->v;
  if (x > 0.0)
    w = (long)(m * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(m * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == start)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%*.5f", (int) w + 7, x);
    *col += w + 8;
  }
}  /* treeout */

///////////////////////////// functions from phylip.c /////////////////



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
  cleerhome();
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




int main(int argc, Char *argv[])
{  /* main program */
#ifdef MAC
  argc = 1;                /* macsetup("Neighbor","");                */
  argv[0] = "Neighbor";
#endif
  init(argc, argv);
  openfile(&infile, argv[1],"input file", "r",argv[0],infilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  datasets = atoi(argv[3]);
  doinit();
  openfile(&outtree,argv[2],"output tree file", "w",argv[0],outtreename);
  ith = 1;
  printf("Please wait patiently while constructing trees......\n\n");

  long repdiv10;

  repdiv10 = datasets / 10;
  if (repdiv10 < 1)
    repdiv10 = 1;
  putchar('\n');

  while (ith <= datasets) {	
    if (datasets < 10 || ith % repdiv10 == 0) {
      printf("Processing data set:    %ld\n",ith);
  }
    inputoptions();
    maketree();
    if (eoln(infile) && (ith < datasets)) 
      scan_eoln(infile);
    ith++;
  }

  FClose(infile);
  FClose(outtree);
  printf("\n\nTree(s) written to file \"%s\"\n\n", outtreename);
  freerest();
  freetree(&curtree.nodep, nonodes2+1);
#ifdef MAC
  fixmacfile(outtreename);
#endif
  printf("Finish.\n\n");
  return 0;
}





