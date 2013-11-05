#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "iostream.h"

#define invroot2 0.7071067814

FILE *fpin, *fpout;
char *infile;

int yu(int a, int b);
void DoEncoding(int ImageData[], char ByteData[], int nx, int ny);
void DoDecoding(int ByteData[], double ByteData[]);
void readppm();
void writeppm();
void encode();
void decode();
void fct2d(double f[], int nrows, int ncols);
void ifct2d(double f[], int nrows, int ncols);

int *block_red=NULL, *block_green=NULL, *block_blue=NULL;
int Ns=0, Nx, Ny, X, Y, color;
int nx, ny, ix, iy;
char type[2];

double quality=5;


static int N=0;
static int m=0;
static double two_over_N=0;
static double root2_over_rootN=0;
static double *C=NULL;

int main(int argc, char **argv)
{
  if(argc!=2){
    printf("\nUsage: ctd file.ppm\n");
    exit(1);
  } 
  infile = argv[1];
  readppm();
  encode();
  decode();
  writeppm();
	return 1;
}

void *quantization(double f[], int nrows, int ncols)
{
  int row, col, nx, i;
  double *ftemp=NULL;

  ftemp = (double *)calloc(nrows*ncols,sizeof(double));

  for (row=0; row<=nrows-1; row++){
    for (col=0; col<=ncols-1; col++){
      if (row+col<=nrows-1){
        nx=(row+col)*(row+col+1)/2+col;
      }else{
        nx=nrows*ncols-(((nrows+ncols)-(row+col)-1)*((nrows+ncols)-(row+col)-2)/2+(ncols-col));
      }
      ftemp[nx] =rint(f[row*ncols+col]/100*quality);
    }
  }
  for (i=0;i<=ncols*nrows-1;i++) f[i] = ftemp[i];
}

//--------------------------------------------------------------
/* read the matrix back to the original matrix*/
double *unquantization(double f[], int nrows, int ncols)
{
  int row, col, nx, i;
  double *ftemp=NULL;

  ftemp = (double *)calloc(nrows*ncols,sizeof(double));

  for (row=0; row<=nrows-1; row++){
    for (col=0; col<=ncols-1; col++){
      if (row+col<=nrows-1){
        nx=(row+col)*(row+col+1)/2+col;
      }else{
        nx=nrows*ncols-(((nrows+ncols)-(row+col)-1)*((nrows+ncols)-(row+col)-2)/2+(ncols-col));
      }
      ftemp[row*ncols+col]=100/quality*f[nx];
    }
  }
  for (i=0;i<=nrows*ncols-1;i++){
      f[i]=ftemp[i];
  }
}

//---------------------------------------------------------------------------------
int yu(int a, int b)
{ int c;
  if (a%b==0) c=b; 
  else c = a%b; 
  return c;
}

void readppm(void)
{
  char rec[200];

  fpin = fopen(infile,"r");
  fscanf(fpin, "%s\n", &type);

  if(strcmp(type,"P3")!=0){
    printf("shit! it's not a ppm file, can't be encoded and decoded by my project\n");
    exit(1);
  }


  fscanf(fpin, "%d %d %d", &X, &Y, &color);


  Nx=int(floor(X/8)+1);
  Ny=int(floor(Y/8)+1);

  block_red = (int *)calloc(X*Y+8*(X+Y), sizeof(int));
  block_green = (int *)calloc(X*Y+8*(X+Y), sizeof(int));
  block_blue = (int *)calloc(X*Y+8*(X+Y), sizeof(int));

  cout <<type<<"  X= "<<X<<" Y= "<<Y<<" color="<<color<<"\n";
  cout <<"Nx="<<Nx<<"  Ny="<<Ny<<"\n";
    while(!feof(fpin)){
      Ns=Ns+1;
      iy = int(floor((Ns-1)/X))%8;
      ny = int(floor(floor((Ns-1)/X)/8));
      ix = yu(yu(Ns,X),8)-1;
      nx = int(floor((Ns-1)%X/8));
      fscanf(fpin,"%d %d %d ", &block_red[(ny*Nx+nx)*64+iy*8+ix],
                              &block_green[(ny*Nx+nx)*64+iy*8+ix],
                              &block_blue[(ny*Nx+nx)*64+iy*8+ix]);
    }
  fclose(fpin);
}

void writeppm(void)
{
   fpout = fopen("decode.ppm","w");
   fprintf(fpout, "%s\n", "P3");
   fprintf(fpout, "%d %d\n", X,Y);
   fprintf(fpout, "%d\n", color);
      
    Ns=0;
    while(nx<Nx&&ny<Ny){
      Ns++;
      iy = int(floor((Ns-1)/X))%8;
      ny = int(floor(floor((Ns-1)/X)/8));
      ix = yu(yu(Ns,X),8)-1;
      nx = int(floor((Ns-1)%X/8));
      fprintf(fpout,"%5d %5d %5d ", abs(block_red[(ny*Nx+nx)*64+iy*8+ix]),
                                    abs(block_green[(ny*Nx+nx)*64+iy*8+ix]),
                                    abs(block_blue[(ny*Nx+nx)*64+iy*8+ix]));

      if(Ns%8 == 0) fprintf(fpout,"\n");
    }
   
   fclose(fpout);
   free(block_red); free(block_green); free(block_blue);
} 
  
void encode(void)
{
    double *point_r=NULL, *point_g=NULL, *point_b=NULL;
    int nrows=8, ncols=8, i;
    int *fint_r=NULL, *fint_g=NULL, *fint_b=NULL;
    int size;
    char *f_byte=NULL;

   fpout=fopen("encode.ctd","wb");

// do the fast cos transform, quantization, and encoding.
// read into buffer
    for (ny=0;ny<=Ny-1;ny++){
    for (nx=0;nx<=Nx-1;nx++){
         point_r = (double *)calloc(64, sizeof(double));
         point_g = (double *)calloc(64, sizeof(double));
         point_b = (double *)calloc(64, sizeof(double));
         for (iy=0;iy<=7; iy++){
         for (ix=0;ix<=7; ix++){
           point_r[iy*8+ix] = 1.0*block_red[(ny*Nx+nx)*64+iy*8+ix];   
           point_g[iy*8+ix] = 1.0*block_green[(ny*Nx+nx)*64+iy*8+ix];
           point_b[iy*8+ix] = 1.0*block_blue[(ny*Nx+nx)*64+iy*8+ix];
         }
         }
      // fct2d, point_color changed to the frequency matrix 
       fct2d(point_r,nrows,ncols);
       fct2d(point_g,nrows,ncols);
       fct2d(point_b,nrows,ncols);
      //quantization the frequency Matrix
       quantization(point_r, nrows, ncols);
       quantization(point_g, nrows, ncols);
       quantization(point_b, nrows, ncols);
      // do the encoding and write to output binary file.
       fint_r= (int *)calloc(64, sizeof(int)); 
       fint_g= (int *)calloc(64, sizeof(int));
       fint_b= (int *)calloc(64, sizeof(int));
       f_byte= (char *)calloc(64, sizeof(char));
      
       for (i=0;i<=63;i++){
         fint_r[i]=(int)point_r[i];
         fint_g[i]=(int)point_g[i];
         fint_b[i]=(int)point_b[i];
       }
    
       DoEncoding(fint_r, f_byte, nx, ny); size = f_byte[0];
       fwrite(f_byte,size,1,fpout);
       DoEncoding(fint_g, f_byte, nx, ny); size = f_byte[0];
       fwrite(f_byte,size,1,fpout);
       DoEncoding(fint_b, f_byte, nx, ny); size = f_byte[0];
       fwrite(f_byte,size,1,fpout);
      
       free(f_byte); free(fint_r); free(fint_g); free(fint_b);
       free(point_r); free(point_g); free(point_b);
    }
    }
  free(block_red); free(block_green); free(block_blue);
  fclose(fpout);
}

void DoEncoding(int ImageData[], char ByteData[], int nx, int ny)
{
  int Index,ByteIndex,NumZero;

  ByteIndex = 3;
  NumZero = 0;
  Index = 0;

  ByteData[0]=0;
  ByteData[1]=nx;
  ByteData[2]=ny;

  while(Index<=63){
    if(ImageData[Index] != 0){
      if(ImageData[Index]>127){
        printf("shit!!!!!!, ImageData[Index] = %d",ImageData[Index]); 
        exit(1);
      }
      ByteData[ByteIndex] = ImageData[Index];
      Index++;
    }
    else{
      ByteData[ByteIndex] = 0;
      ByteIndex++;
      NumZero=0;
      while (ImageData[Index] == 0){
        Index++;
        NumZero++;
      }
      ByteData[ByteIndex] = NumZero;
    }
      ByteIndex++;
  }
  ByteData[ByteIndex-1]=ByteData[ByteIndex-1]-1; //a bug here, the last guy always 1 bigger if not 0
                                                 //so I just subtracted it by one;
  ByteData[0]=ByteIndex;
}

void decode()
{
    int nrows=8, ncols=8;
    double *point_r=NULL, *point_g=NULL, *point_b=NULL;
    double *f64;
    int  *fint=NULL;
    int  i, count=0;
    char size;
    char *f_byte=NULL;

      block_red = (int *)calloc(X*Y+8*(X+Y), sizeof(int));
      block_green = (int *)calloc(X*Y+8*(X+Y), sizeof(int));
      block_blue = (int *)calloc(X*Y+8*(X+Y), sizeof(int));

   fpin=fopen("encode.ctd","r");
   while(!feof(fpin)){
      point_r = (double *)calloc(64, sizeof(double));
      point_g = (double *)calloc(64, sizeof(double));
      point_b = (double *)calloc(64, sizeof(double));
      fint = (int *)calloc(64, sizeof(int));
      fscanf(fpin,"%c",&size);
      f_byte = (char *)calloc(size, sizeof(char));
      fread(f_byte, 1, size-1, fpin);
      fint[0]=size; 
      for (i=0; i<=size-2; i++){
        fint[i+1] = f_byte[i];
        if (fint[i]>=127 ) { printf("***shit, fint[i]= %d\n", fint[i]); exit(1);}
      } 
      switch (count%3){
        case 0: DoDecoding(fint, point_r);
                unquantization(point_r, nrows, ncols); 
                ifct2d(point_r,nrows,ncols);
                for (iy=0;iy<=7; iy++){ for (ix=0;ix<=7; ix++){
                  block_red[(ny*Nx+nx)*64+iy*8+ix] = (int)point_r[iy*8+ix];}}
                break;
        case 1: DoDecoding(fint, point_g); 
                unquantization(point_g, nrows, ncols);
                ifct2d(point_g,nrows,ncols);
                for (iy=0;iy<=7; iy++){ for (ix=0;ix<=7; ix++){
                  block_green[(ny*Nx+nx)*64+iy*8+ix] = (int)point_g[iy*8+ix];}}
                break;
        case 2: DoDecoding(fint, point_b); 
                unquantization(point_b, nrows, ncols);
                ifct2d(point_b,nrows,ncols);
                for (iy=0;iy<=7; iy++){ for (ix=0;ix<=7; ix++){
                  block_blue[(ny*Nx+nx)*64+iy*8+ix] = (int)point_b[iy*8+ix];}}
                break; 
      }
      count++;
      if(ny==Ny) break;
      free(f_byte); free(fint); 
      free(point_r); free(point_g); free(point_b);
   }//end of while

   fclose(fpin);
}

//This function transform the encoded bytes vector to 8 x 8 matrix.
void DoDecoding(int ByteData[], double ImageData[]) //infact, the ByteData is integer vector
{
  int size;
  int ByteIndex, Index, ZeroIndex;
  int i;

  size = ByteData[0];
  nx = ByteData[1];
  ny = ByteData[2];
  ByteIndex = 3;
  Index = 0;

  while (ByteIndex <= size-1) {
    if (ByteData[ByteIndex] != 0) {
      ImageData[Index] = ByteData[ByteIndex]*1.0;
      ByteIndex = ByteIndex + 1;
      Index = Index + 1;
    } else {
      for (ZeroIndex = 0; ZeroIndex<= ByteData[ByteIndex + 1] - 1; ZeroIndex++) {
        ImageData[Index] = 0;
        Index = Index + 1;
      }
       ByteIndex = ByteIndex + 2;
    }
  }
}

static void bitrev(double *f, int len)
{
  int i,j,m,halflen;
  double temp;

  if (len<=2) return; /* No action necessary if n=1 or n=2 */
  halflen = len>>1;
  j=1;
  for(i=1; i<=len; i++){
    if(i<j){
      temp=f[j-1];
      f[j-1]=f[i-1];
      f[i-1]=temp;
    }
    m = halflen;
    while(j>m){
      j=j-m;
      m=(m+1)>>1;
    }
    j=j+m;
  }
}

static void inv_sums(double *f)
{
  int ii,stepsize,stage,curptr,nthreads,thread,step,nsteps;

  for(stage=1; stage <=m-1; stage++){
    nthreads = 1<<(stage-1);
    stepsize = nthreads<<1;
    nsteps   = (1<<(m-stage)) - 1;
    for(thread=1; thread<=nthreads; thread++){
      curptr=N-thread;
      for(step=1; step<=nsteps; step++){
        f[curptr] += f[curptr-stepsize];
        curptr -= stepsize;
      }
    }
  }
}

static void fwd_sums(double *f)
{
  int ii,stepsize,stage,curptr,nthreads,thread,step,nsteps;

  for(stage=m-1; stage >=1; stage--){
    nthreads = 1<<(stage-1);
    stepsize = nthreads<<1;
    nsteps   = (1<<(m-stage)) - 1;
    for(thread=1; thread<=nthreads; thread++){
      curptr=nthreads +thread-1;
      for(step=1; step<=nsteps; step++){
        f[curptr] += f[curptr+stepsize];
        curptr += stepsize;
      }
    }
  }
}

static void scramble(double *f,int len){
  double temp;
  int i,ii1,ii2,halflen,qtrlen;

  halflen = len >> 1;
  qtrlen = halflen >> 1;
  bitrev(f,len);
  bitrev(&f[0], halflen);
  bitrev(&f[halflen], halflen);
  ii1=len-1;
  ii2=halflen;
  for(i=0; i<=qtrlen-1; i++){
    temp = f[ii1];
    f[ii1] = f[ii2];
    f[ii2] = temp;
    ii1--;
    ii2++;
  }
}

static void unscramble(double *f,int len)
{
  double temp;
  int i,ii1,ii2,halflen,qtrlen;

  halflen = len >> 1;
  qtrlen = halflen >> 1;
  ii1 = len-1;
  ii2 = halflen;
  for(i=0; i<=qtrlen-1; i++){
    temp = f[ii1];
    f[ii1] = f[ii2];
    f[ii2] = temp;
    ii1--;
    ii2++;
  }
  bitrev(&f[0], halflen);
  bitrev(&f[halflen], halflen);
  bitrev(f,len);
}

static void initcosarray(int length)
{
  int i,group,base,item,nitems,halfN;
  double factor;

  m = -1;
  do{
    m++;
    N = 1<<m;
    if (N>length){
      printf("ERROR in FCT-- length %d not a power of 2\n",length);
      exit(1);
    }
  }while(N<length);
  if(C != NULL) free(C);
  C = (double *)calloc(N,sizeof(double));
  if(C == NULL){
    printf("Unable to allocate C array\n");
    exit(1);
  }
  halfN=N/2;
  two_over_N = 2.0/(double)N;
  root2_over_rootN = sqrt(2.0/(double)N);
  for(i=0;i<=halfN-1;i++) C[halfN+i]=4*i+1;
  for(group=1;group<=m-1;group++){
    base= 1<<(group-1);
    nitems=base;
    factor = 1.0*(1<<(m-group));
    for(item=1; item<=nitems;item++) C[base+item-1]=factor*C[halfN+item-1];
  }

  //printf("before taking cos, C array =\n"); rarrwrt(C,N);
  for(i=1;i<=N-1;i++) C[i] = 1.0/(2.0*cos(C[i]*M_PI/(2.0*N)));
  //printf("After taking cos, Carray = \n"); rarrwrt(C,N);
}

static void inv_butterflies(double *f)
{
  int stage,ii1,ii2,butterfly,ngroups,group,wingspan,increment,baseptr;
  double Cfac,T;

  for(stage=1; stage<=m;stage++){
    ngroups=1<<(m-stage);
    wingspan=1<<(stage-1);
    increment=wingspan<<1;
    for(butterfly=1; butterfly<=wingspan; butterfly++){
      Cfac = C[wingspan+butterfly-1];
      baseptr=0;
      for(group=1; group<=ngroups; group++){
        ii1=baseptr+butterfly-1;
        ii2=ii1+wingspan;
        T=Cfac * f[ii2];
        f[ii2]=f[ii1]-T;
        f[ii1]=f[ii1]+T;
        baseptr += increment;
      }
    }
  }
}

static void fwd_butterflies(double *f)
{
  int stage,ii1,ii2,butterfly,ngroups,group,wingspan,increment,baseptr;
  double Cfac,T;

  for(stage=m; stage>=1;stage--){
    ngroups=1<<(m-stage);
    wingspan=1<<(stage-1);
    increment=wingspan<<1;
    for(butterfly=1; butterfly<=wingspan; butterfly++){
      Cfac = C[wingspan+butterfly-1];
      baseptr=0;
      for(group=1; group<=ngroups; group++){
        ii1=baseptr+butterfly-1;
        ii2=ii1+wingspan;
        T= f[ii2];
        f[ii2]=Cfac *(f[ii1]-T);
        f[ii1]=f[ii1]+T;
        baseptr += increment;
      }
    }
  }
}

static void ifct_noscale(double *f, int length)
{
  if (length != N) initcosarray(length);
  f[0] *= invroot2;
  inv_sums(f);
  bitrev(f,N);
  inv_butterflies(f);
  unscramble(f,N);
}

static void fct_noscale(double *f, int length)
{
  if (length != N) initcosarray(length);
  scramble(f,N);
  fwd_butterflies(f);
  bitrev(f,N);
  fwd_sums(f);
  f[0] *= invroot2;
}

static void ifct_defn_scaling(double *f, int length){
  ifct_noscale(f,length);
}

static void fct_defn_scaling(double *f, int length){
  int i;

  fct_noscale(f,length);
  for(i=0;i<=N-1;i++) f[i] *= two_over_N;
}

void ifct(double *f, int length){
/* CALL THIS FOR INVERSE 1D DCT DON-MONRO PREFERRED SCALING */
  int i;

  if (length != N) initcosarray(length);  /* BGS patch June 1997 */
  for(i=0;i<=N-1;i++) f[i] *= root2_over_rootN;
  ifct_noscale(f,length);
}

void fct(double *f, int length){
/* CALL THIS FOR FORWARD 1D DCT DON-MONRO PREFERRED SCALING */
  int i;

  fct_noscale(f,length);
  for(i=0;i<=N-1;i++) f[i] *= root2_over_rootN;
}

/****************************************************************
    2D FAST DCT SECTION
****************************************************************/

#define VERBOSE 0

static double *g = NULL;
static double two_over_sqrtncolsnrows = 0.0;
static int ncolsvalue = 0;
static int nrowsvalue = 0;

static void initfct2d(int nrows, int ncols){
  if(VERBOSE) printf("FCT2D -- Initialising for new nrows=%d\n",nrows);
  if ((nrows<=0)||(ncols<0)){
    printf("FCT2D -- ncols=%d or nrows=%d is <=0\n",nrows,ncols);
    exit(1);
  }
  if(g != NULL) free(g);
  g = (double *)calloc(nrows,sizeof(double));
  if(g == NULL){
    printf("FCT2D -- Unable to allocate g array\n");
    exit(1);
  }
  ncolsvalue = ncols;
  nrowsvalue = nrows;
  two_over_sqrtncolsnrows = 2.0/sqrt(ncols*1.0*nrows);
}

void fct2d(double f[], int nrows, int ncols)
/* CALL THIS FOR FORWARD 2d DCT DON-MONRO PREFERRED SCALING */
{
  int u,v;

  if ((ncols!=ncolsvalue)||(nrows!=nrowsvalue)){
    initfct2d(nrows,ncols);
  }
  for (u=0; u<=nrows-1; u++){
    fct_noscale(&f[u*ncols],ncols);
  }
  for (v=0; v<=ncols-1; v++){
    for (u=0; u<=nrows-1; u++){
       g[u] = f[u*ncols+v];
    }
    fct_noscale(g,nrows);
    for (u=0; u<=nrows-1; u++){
      f[u*ncols+v] = g[u]*two_over_sqrtncolsnrows;
    }
  }
}

void ifct2d(double f[], int nrows, int ncols)
/* CALL THIS FOR INVERSE 2d DCT DON-MONRO PREFERRED SCALING */
{
  int u,v;

  if ((ncols!=ncolsvalue)||(nrows!=nrowsvalue)){
    initfct2d(nrows,ncols);
  }
  for (u=0; u<=nrows-1; u++){
    ifct_noscale(&f[u*ncols],ncols);
  }
  for (v=0; v<=ncols-1; v++){
    for (u=0; u<=nrows-1; u++){
       g[u] = f[u*ncols+v];
    }
    ifct_noscale(g,nrows);
    for (u=0; u<=nrows-1; u++){
       f[u*ncols+v] = g[u]*two_over_sqrtncolsnrows;
    }
  }
}


static void rarrwrt2d(double f[],int nrows,int ncols)
{
  int row,col;

  for (row=0;row<=nrows-1; row++){
    for (col=0;col<=ncols-1; col++){
      printf("%5.0f  ",f[row*ncols+col]);
    }
    printf("\n");
  }printf("\n");
}

static void rarrwrti2d(int f[],int nrows,int ncols)
{
  int row,col;

  for (row=0;row<=nrows-1; row++){
    for (col=0;col<=ncols-1; col++){
      printf("%5i  ",f[row*ncols+col]);
    }
    printf("\n");
  }
}

