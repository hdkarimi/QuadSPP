/*------------------------------------------------------------------------------
* gnsscmn.c : program common functions
*
*
* references :
*	  (1) IS-GPS-200H, Navstar GPS Space Segment/Navigation User Segment 
*		  Interfaces, 2015.
*     (2) European GNSS (Galileo) open service, Signal In Space Interface 
*         Control Document, 2016.
*     (3) W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 3.03, 2015.
*	  (4) A.E.Niell, Global mapping functions for the atmosphere delay at radio
*         wavelengths, Jounal of geophysical research, 1996.
*	  (5) China Satellite Navigation Office, BeiDou navigation satellite system
*         signal in space interface control document, open service signal B1I
*         (version 1.0), Dec 2012
*	  (6) O.Montenbruck and E.Gill, Satellite Orbit, 2001.
*
* version : yyyy/mm/dd		ver		description
*			2018/03/25		1.0		create  
*			2018/04/27		1.1		fix bug on NMF 
*
* Copyright (C) 2018 by H.KARIMI, All rights reserved.
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include <ctype.h>

#include "gnsslib.h"

/* constants -----------------------------------------------------------------*/
static const double gpst0[]={1980,1, 6,0,0,0}; /* gps time origin     */
static const double gst0 []={1999,8,22,0,0,0}; /* galileo time origin */
static const double bdt0 []={2006,1, 1,0,0,0}; /* beidou time origin  */

static double leaps[MAXLEAPS+1][7] = {/* leap seconds (y,m,d,h,m,s,utc-gpst)*/
	{2017,1,1,0,0,0,-18},
	{2015,1,1,0,0,0,-17},
	{2012,7,1,0,0,0,-16},
    {2009,1,1,0,0,0,-15},
    {2006,1,1,0,0,0,-14},
    {1999,1,1,0,0,0,-13},
    {1997,7,1,0,0,0,-12},
    {1996,1,1,0,0,0,-11},
    {1994,7,1,0,0,0,-10},
    {1993,7,1,0,0,0, -9},
    {1992,7,1,0,0,0, -8},
    {1991,1,1,0,0,0, -7},
    {1990,1,1,0,0,0, -6},
    {1988,1,1,0,0,0, -5},
    {1985,7,1,0,0,0, -4},
    {1983,7,1,0,0,0, -3},
    {1982,7,1,0,0,0, -2},
    {1981,7,1,0,0,0, -1},
    {0}
};

static char *obscodes[]={		/*observation code string */
	""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */ 
	"1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
	"2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
	"6A","6B","6C","6X","6Z", "6S","6L","8I","8Q","8X", /* 30-39 */
	"2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q",""    /* 40-49 */
};
static unsigned char obsfreqs[] = {  /* 1:L1, 2:L2, 3:L5, 4:L6, 5:L7, 6:L8 */

	0, 1, 1, 1, 1,  1, 1, 1, 1, 1, /*  0- 9 */
	1, 1, 1, 1, 2,  2, 2, 2, 2, 2, /* 10-19 */
	2, 2, 2, 2, 3,  3, 3, 5, 5, 5, /* 20-29 */
	4, 4, 4, 4, 4,  4, 4, 6, 6, 6, /* 30-39 */
	2, 2, 4, 4, 3,  3, 3, 1, 1, 0  /* 40-49 */
};
static char codepris[4][MAXFREQ][10]={  /* code priority table */
   
   /* L1/E1      L2/B1        L5/E5a/L3 L6/B3	  E5b/B2    E5(a+b) */
    {"CPYWMNSL","PYWCMNDSLX","IQX"     ,""       ,""       ,""      }, /* GPS */
    {"PC"      ,"PC"        ,"IQX"     ,""       ,""       ,""      }, /* GLO */
    {"CABXZ"   ,""          ,"IQX"     ,"ABCXZ"  ,"IQX"    ,"IQX"   }, /* GAL */
    {"IQX"     ,"IQX"       ,"IQX"     ,"IQX"    ,"IQX"    ,""      }  /* BDS */
};

/* 1_SATELLITE NUMBER/SYSTEM FUNCTIONS ---------------------------------------*/

/* satellite system+prn/slot number to satellite number ------------------------
* convert satelite system+prn/slot number to satellite  number
* args	 : int		sys		I	satellite system (SYS_GPS,SYS_GLO,...)
*		   int		prn		I	satellite prn/slot number
* return : satellite number (0:error)
*-----------------------------------------------------------------------------*/
extern int satno(int sys, int prn)
{
	if(prn<=0) return 0;
	switch(sys) {
		case SYS_GPS:
			if (prn<MINPRNGPS || MAXPRNGPS<prn) return 0;
			return prn-MINPRNGPS+1;
		case SYS_GLO:
			if (prn<MINPRNGLO || MAXPRNGLO<prn) return 0;
			return NSATGPS + prn-MINPRNGLO+1;
		case SYS_GAL:
			if (prn<MINPRNGAL || MAXPRNGAL<prn) return 0;
			return NSATGPS+NSATGLO + prn-MINPRNGAL+1;
		case SYS_BDS:
			if (prn<MINPRNBDS || MAXPRNBDS<prn) return 0;
			return NSATGPS+NSATGLO+NSATGAL + prn-MINPRNBDS+1;
	}
	return 0;
}
/* satellite number to satellite system ----------------------------------------
* convert satellite number to satellite system                              
* args	 : int		sat		I	satellite number (1-MAXSAT)
*          int		*prn	IO	satellite prn/slot number (NULL: no output)
* return : satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/
extern int satsys(int sat, int *prn)
{
	int sys=SYS_NON;

	if (sat<=0 || MAXSAT<sat) sat=0;
	else if (sat<=NSATGPS) {
		sys=SYS_GPS; sat+=MINPRNGPS-1;
	}
	else if ((sat-=NSATGPS)<=NSATGLO) {
		sys=SYS_GLO; sat+=MINPRNGLO-1;
	}
	else if ((sat-=NSATGLO)<=NSATGAL) {
		sys=SYS_GAL; sat+=MINPRNGAL-1;
	}
	else if ((sat-=NSATGAL)<=NSATBDS) {
		sys=SYS_BDS; sat+=MINPRNBDS-1;
	}
	else sat=0;

	if (prn) *prn=sat;
	return sys;
}
/* satellite id to satellite number --------------------------------------------
* convert satellite id to satellite number 
* args	 : char		*id		I	satellite id (Gnn,Rnn,Enn,Cnn)
*return  : satellite number (0:error)
*-----------------------------------------------------------------------------*/
extern int satid2no(const char *id)
{
	int sys, prn;
	char code;

	if (sscanf(id,"%c%d",&code,&prn)<2) return 0;

	switch (code) {
		case 'G': sys=SYS_GPS; prn+=MINPRNGPS-1; break;
		case 'R': sys=SYS_GLO; prn+=MINPRNGLO-1; break;
		case 'E': sys=SYS_GAL; prn+=MINPRNGAL-1; break;
		case 'C': sys=SYS_BDS; prn+=MINPRNBDS-1; break;
		default: return 0;
	}
	return satno(sys, prn);
}
/* satellite number to satellite id --------------------------------------------
* convert satellite number to satellite id
* args	 : int		sat		I	satellite number
*		   char		*id		O	satellite id (Gnn,Rnn,Enn,Cnn)
*-----------------------------------------------------------------------------*/
extern void satno2id(int sat, char *id)
{
	int prn;

	switch (satsys(sat, &prn)) {
		case SYS_GPS: sprintf(id,"G%02d",prn-MINPRNGPS+1); return;
		case SYS_GLO: sprintf(id,"R%02d",prn-MINPRNGLO+1); return;
		case SYS_GAL: sprintf(id,"E%02d",prn-MINPRNGAL+1); return;
		case SYS_BDS: sprintf(id,"C%02d",prn-MINPRNBDS+1); return;
	}
	strcpy(id, "");
}
/* obs type string to obs code -------------------------------------------------
* convert obs code type string to obs code
* arg	 :	char	*str	I	obs code string ("1C","1P","1Y",...)
*			int		*freq	IO	frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8, 0:err)
*								(NULL: no output)
* return : obs code (CODE_xxx)
*-----------------------------------------------------------------------------*/
extern unsigned char obs2code(const char *obs, int *freq)
{
	int i;
	if (freq) *freq=0;
	for (i=1;*obscodes[i];i++) {
		if (strcmp(obscodes[i],obs)) continue;
		if (freq) *freq=obsfreqs[i];
		return (unsigned char)i;
	}
	return CODE_NON;
}
/* obs code to obs code string -------------------------------------------------
* convert obs code to obs code sting 
* args	 : unsigned char code	I	obs code (CODE_xxx)
*		   int		*freq		IO	frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,0:err)
*									(NULL: no output)
* return : obs code string ("1C","1P","1Y",...)
*-----------------------------------------------------------------------------*/
extern char *code2obs(unsigned char code, int *freq)
{
	if (freq) *freq=0;
	if (code<=CODE_NON||MAXCODE<code) return "";
	if (freq) *freq=obsfreqs[code];
	return obscodes[code]; 
}
/* get code priority -----------------------------------------------------------
* get code priority for multiple codes in a frequency
* args   : int    sys     I     system (SYS_???)
*          unsigned char code I obs code (CODE_???)
*          char   *opt    I     code options (NULL:no option)
* return : priority (15:highest-1:lowest,0:error)
*-----------------------------------------------------------------------------*/
extern int getcodepri(int sys, unsigned char code, const char *opt)
{
    const char *p,*optstr;
    char *obs,str[8]="";
    int i,j;
    
    switch (sys) {
        case SYS_GPS: i=0; optstr="-GL%2s"; break;
        case SYS_GLO: i=1; optstr="-RL%2s"; break;
        case SYS_GAL: i=2; optstr="-EL%2s"; break;
        case SYS_BDS: i=3; optstr="-CL%2s"; break;
        default: return 0;
    }
    obs=code2obs(code,&j);
    
    /* parse code options */
    for (p=opt;p&&(p=strchr(p,'-'));p++) {
        if (sscanf(p,optstr,str)<1||str[0]!=obs[0]) continue;
        return str[1]==obs[1]?15:0;
    }
    /* search code priority */
    return (p=strchr(codepris[i][j-1],obs[1]))?14-(int)(p-codepris[i][j-1]):0;
}

/* 2_MATRIX AND VECTOR FUNCTIONS ---------------------------------------------*/

/* new matrix ------------------------------------------------------------------ 
* allocate meomory of matrix
* args	 : int		n,m,		I	number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0 retuen NULL)
*-----------------------------------------------------------------------------*/
extern double *mat(int n, int m)
{
	double *p;

	if (n<=0 || m<=0) return NULL;
	p = (double *)malloc(sizeof(double)*n*m);
	
	return p;
}
/* new integer matrix ----------------------------------------------------------
* allocate memory of integer matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern int *imat(int n, int m)
{
	int *p;

	if (n <= 0 || m <= 0) return NULL;
	p = (int *)malloc(sizeof(int)*n*m);

	return p;
}
/* zero matrix -----------------------------------------------------------------
* create new zero matrix
* args	 : int		n,m		I	number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *zeros(int n, int m)
{
	double *p;

	if (n<=0 || m<=0) return NULL;
	p=(double *)calloc(sizeof(double), n*m); 

	return p;
}
/* inner product ---------------------------------------------------------------
* inner product of vectors
* args	 : double *a,*b		I	vector a,b (nx1)
*		   int	  n			I 	size of vector a,b
* return a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
	double c=0.0;

	while (--n>=0)
		c+=a[n]*b[n];
	return c;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args	 : double	*a		I	vector a (nx1) 
*		   int		n		I	size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a, int n)
{
	return sqrt(dot(a,a,n));
}
/* copy matrix -----------------------------------------------------------------
* copy matrix 
* args	 : double	*A		O	destination matrix A (n x m)
*		   double	*B		I	source matrix B (n x m)
*		   int		n,m		I	number of rows and columns of matrix
* return : none
*-----------------------------------------------------------------------------*/
extern void matcpy(double *A, const double *B, int n, int m)
{
	memcpy(A,B,sizeof(double)*n*m);
}

/* multiply matrix -------------------------------------------------------------
* multiply matrix by matrix (C=alpha*A*B+beta*C);
* args	 : char		*tr		I	transpose flags ("N": Normal, "T":Transpose)
*		   int		n,k,m	I	size of (transposed) matrix A, B
*		   double	alpha	I	alpha
*		   double	*A,*B	I	(transposed) matrix A (n x m), B (m x k)
*		   double	beta	I	beta
*		   double	*C		IO	matrix	C (n x k)
* return : none
*-----------------------------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
	               const double *A, const double *B, double beta, double *C)
{
	double d;
	int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

	for (i=0;i<n;i++)
		for (j=0;j<k;j++) {
			d=0.0;
			switch (f) {
				case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
				case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
				case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
				case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
			}
			if  (beta==0.0) 
				C[i+j*n]=alpha*d;
			else 
				C[i+j*n]=alpha*d+beta*C[i+j*n];
		}
}
/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
{
	double big,s,temp,*vv=mat(n,1);
	int i,j,k,imax=0;

	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0; for (j=0;j<n;j++) if ((temp=fabs(A[i+j*n]))>big) big=temp;
		if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
			if ((temp=vv[i]*fabs(s))>=big) {big=temp; imax=i;}
		}
		if (j!=imax) {
			for (k=0;k<n;k++) {
				temp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=temp;
			}
			*d=-(*d); vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (A[j+j*n]==0.0) {free(vv); return -1;}
		if (j!=n-1) {
			temp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=temp;
		}
	}
	free(vv);
	return 0;
} 
/* LU back substitution ------------------------------------------------------*/
static void lubksb(const double *A, int n, const int *indx, double *b)
{
	double s;
	int i,ii=-1,ip,j;

	for (i=0;i<n;i++) {
		ip=indx[i]; s=b[ip]; b[ip]=b[i];
		if(ii>=00) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
		b[i]=s;
	}
	for (i=n-1;i>=0;i--) {
		s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
	}
}
/* inverse of matrix -----------------------------------------------------------
* inverse of matrix (A=A^-1)
* args	 : double	*A		IO	matrix (n x n)
*		   int		n		I	size of matrix A
* return : status (0:ok, 0>:err)
*-----------------------------------------------------------------------------*/
extern int matinv(double *A, int n)
{
	double d, *B;
	int i,j, *indx;

	indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
	if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
	for (j=0;j<n;j++) {
		for (i=0;i<n;i++) 
			A[i+j*n]=0.0;

		A[j+j*n]=1.0;
		lubksb(B,n,indx,A+j*n);
	}
	free(indx); free(B);
	return 0;
}
/* least square estimation -----------------------------------------------------
* least square estimation by solving normal equation (x=(A'*A)^-1*A'y)
* args	 : double	*A		I	transpose of (weighted) design matrix (n x m)
*		   double	*y		I	(weighted) measurements (m x 1)
*		   int		n,m		I	number of parameters and measurements (n<=m)
*		   double	*x		O	estimated parameters	(n x 1)
*		   double	*Q		O	estimated parameters covariance matrix (n x n)
* return : status (0:ok,0>:err)
* notes  : for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
*		   matrix stored by column-major order
*-----------------------------------------------------------------------------*/
extern int lsq(const double *A, const double *y, int n, int m, double *x,
			   double *Q)
{
	double *Ay;
	int info;

	if (m<n) return -1;
	Ay=mat(n,1);
	matmul("NN",n,1,m,1.0,A,y,0.0,Ay); /* Ay=A*y */
	matmul("NT",n,n,m,1.0,A,A,0.0,Q);  /* Q=A*A' */
	if(!(info=matinv(Q,n))) matmul("NN",n,1,n,1.0,Q,Ay,0.0,x); /* x= Q^-1*Ay */
	free(Ay);
	return info;
}

/* 3_TIME AND STRING FUNCTIONS -----------------------------------------------*/

/* string to number ------------------------------------------------------------
* convert substring in string to number 
* args	 : char		*s		I	string ("... nn.nn ...")
*		   int		i,n		I	substring position and width
* return : converted number (0.0:err)
*-----------------------------------------------------------------------------*/
extern double str2num(const char *s, int i, int n)
{
	double val;
	char str[256],*p=str;

	if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<n) return 0.0;
	for (s+=i;*s&&--n>=0;s++) *p++=*s=='d'||*s=='D'?'E':*s;
	*p='\0';
	return sscanf(str,"%lf",&val)==1?val:0.0;
}
/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct 
* args	 : char		*s		I	string ("... yyyy mm dd hh mm ss ...")
*		   int		i,n		I	substring position and width
*		   gtime_t	*t		O	gtime_t struct
* retirn status (0:ok,0>:err)
*-----------------------------------------------------------------------------*/
extern int str2time(const char *s, int i, int n, gtime_t *t)
{
	double ep[6];
	char str[256],*p=str;

	if(i<0||(int)strlen(s)<i||(int)sizeof(str)-1<i) return -1;
	for(s+=i;*s&&--n>=0;) *p++=*s++;
	*p='\0';
	if (sscanf(str,"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
		return -1;
	if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
	*t=epoch2time(ep);
	return 0;
}/* calendar day/time to time --------------------------------------------------
* convert calendar day/time to gtime_t struct
* args	 : double	*ep		I	day/time {year,month,day,hour,min,sec}
* return : gtime_t struct
* notes  : proper in 1970-2037 (32bit time_t) or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern gtime_t epoch2time(const double *ep)
{
	const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
	gtime_t time={0};
	int days,sec, year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];

	if (year<1970||2099<year||mon<1||12<mon) return time;

	/* leap year if year%4==0 in 1901-2099*/
	days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
	sec=(int)floor(ep[5]);
	time.time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
	time.sec=ep[5]-sec;
	return time;
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time 
* args	 : gtime_t t		I	gtime_t struct
*		   double *ep		O	day/time {year,month,day,hour,min,sec}
* return : none 
* notes  : proper in 1970-2037 (32bit time_t) or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern void time2epoch(gtime_t t, double *ep)
{
	const int mday[]={ /* # of days in a month */ 
		31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
	};

	int days,sec,mon,day;

	/* leap year if year%4==0 in 1901-2099 */
	days=(int)(t.time/86400);
	sec=(int)(t.time-(time_t)days*86400);
	for (day=days%1461,mon=0;mon<48;mon++) {
		if (day>=mday[mon]) day-=mday[mon]; else break;
	}

	ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
	ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t.sec;
}
/* gps time to time ------------------------------------------------------------
* convert week and tow in gps time to gtime_t struct
* args	 : int		week		I	week number in gps time
*          double	sec			I	time of week in gps time (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2time(int week, double sec)
{
	gtime_t t=epoch2time(gpst0);

	if (sec<-1E9||1E9<sec) sec=0.0;
	t.time+=(time_t)86400*7*week+(int)sec;
	t.sec=sec-(int)sec;

	return t;
}
/* time to gps time ------------------------------------------------------------
* convert gtime_t struct to week and tow in gps time
* args	 : gtime_t	t		I	gtime_t struct
*          int		*week	IO	week number in gps time (NULL: no output)
* return : time of week in gps time (s)
*-----------------------------------------------------------------------------*/
extern double time2gpst(gtime_t t, int *week)
{
	gtime_t t0=epoch2time(gpst0);
	time_t sec=t.time-t0.time;
	int w=(int)(sec/(86400*7));

	if(week) *week=w;
	return (double)(sec-(double)w*86400*7)+t.sec;
}
/* add time --------------------------------------------------------------------
* add time to gtime_t struct
* args	 : gtime_t	t		I	gtime_t struct
*          double sec		I	time to add (s)
* return : gtime_t struct (t+sec)
*-----------------------------------------------------------------------------*/
extern gtime_t  timeadd(gtime_t t, double sec)
{
	double tt;

	t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
	return t;
}
/* time difference -------------------------------------------------------------
* difference between gtime_t struct
* args	 : gtime_t	t1,t2		I	gtime_t struct
* return : time difference (t1-t2) (s)
*-----------------------------------------------------------------------------*/
extern double timediff(gtime_t t1, gtime_t t2)
{
	return difftime(t1.time,t2.time) + t1.sec-t2.sec;
}
/* gpstime to utc --------------------------------------------------------------
* convert gpstime to utc considering leap seconds
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in utc
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2utc(gtime_t t)
{
    gtime_t tu;
    int i;
    
    for (i=0;leaps[i][0]>0;i++) {
        tu=timeadd(t,leaps[i][6]);
        if (timediff(tu,epoch2time(leaps[i]))>=0.0) return tu;
    }
    return t;
}
/* utc to gpstime --------------------------------------------------------------
* convert utc to gpstime considering leap seconds
* args   : gtime_t t        I   time expressed in utc
* return : time expressed in gpstime
*-----------------------------------------------------------------------------*/
extern gtime_t utc2gpst(gtime_t t)
{
    int i;
    
    for (i=0;leaps[i][0]>0;i++) {
        if (timediff(t,epoch2time(leaps[i]))>=0.0) return timeadd(t,-leaps[i][6]);
    }
    return t;
}
/* time to string --------------------------------------------------------------
* convert gtime_t struct to string
* args	 : gtime_t	t		I	gtime_t struct
*		   char		*s		O	string ("yyyy/mm/dd hh:mm:ss.sssss")
*          int		n		I	number of decimals
*-----------------------------------------------------------------------------*/
extern void time2str(gtime_t t, char *s, int n)
{
	double ep[6];

	if(n<0) n=0; else if (n>12) n=12;
	if (1.0-t.sec<0.5/pow(10.0,n)) {t.time++; t.sec=0.0;};
	time2epoch(t,ep);
	sprintf(s,"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f",ep[0],ep[1],ep[2],
		    ep[3],ep[4],n<=0?2:n+3,n<=0?0:n,ep[5]);
}
/* time to day of year ---------------------------------------------------------
* convert gtime_t struct to day of year
* args   : gtime_t t        I   gtime_t struct
* return : day of year (days)
*-----------------------------------------------------------------------------*/
extern double time2doy(gtime_t t)
{
    double ep[6];
    
    time2epoch(t,ep);
    ep[1]=ep[2]=1.0; ep[3]=ep[4]=ep[5]=0.0;
    return timediff(t,epoch2time(ep))/86400.0+1.0;
}
/* screen by time --------------------------------------------------------------
* screening by time start, time end, and time interval
* args   : gtime_t time  I      time
*          gtime_t ts    I      time start (ts.time==0:no screening by ts)
*          gtime_t te    I      time end   (te.time==0:no screening by te)
*          double  tint  I      time interval (s) (0.0:no screen by tint)
* return : 1:on condition, 0:not on condition
*-----------------------------------------------------------------------------*/
extern int screent(gtime_t time, gtime_t ts, gtime_t te, double tint)
{
    return (tint<=0.0||fmod(time2gpst(time,NULL)+DTTOL,tint)<=DTTOL*2.0)&&
           (ts.time==0||timediff(time,ts)>=-DTTOL)&&
           (te.time==0||timediff(time,te)<  DTTOL);
}
/* get time string -------------------------------------------------------------
* get time string
* args   : gtime_t t        I   gtime_t struct
*          int    n         I   number of decimals
* return : time string
* notes  : not reentrant, do not use multiple in a function
*-----------------------------------------------------------------------------*/
extern char *time_str(gtime_t t, int n)
{
    static char buff[64];
    time2str(t,buff,n);
    return buff;
}

/*4_COORSINATES FUNCTIONS ----------------------------------------------------*/

/* convert degree to deg-min-sec -----------------------------------------------
* convert degree to degree-minute-second
* args   : double deg       I   degree
*          double *dms      O   degree-minute-second {deg,min,sec}
*          int    ndec      I   number of decimals of second
* return : none
*-----------------------------------------------------------------------------*/
extern void deg2dms(double deg, double *dms, int ndec)
{
    double sign=deg<0.0?-1.0:1.0, a=fabs(deg);
    double unit=pow(0.1,ndec);

    dms[0]=floor(a); a=(a-dms[0])*60.0;
    dms[1]=floor(a); a=(a-dms[1])*60.0;
    dms[2]=floor(a/unit+0.5)*unit;
    if (dms[2]>=60.0) {
        dms[2]=0.0;
        dms[1]+=1.0;
        if (dms[1]>=60.0) {
            dms[1]=0.0;
            dms[0]+=1.0;
        }
    }
    dms[0]*=sign;
}
/* convert deg-min-sec to degree -----------------------------------------------
* convert degree-minute-second to degree
* args   : double *dms      I   degree-minute-second {deg,min,sec}
* return : degree
*-----------------------------------------------------------------------------*/
extern double dms2deg(const double *dms)
{
    double sign=dms[0]<0.0?-1.0:1.0;
    return sign*(fabs(dms[0])+dms[1]/60.0+dms[2]/3600.0);
}
/* transform ecef to geodetic position ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height, see ref (6)
*-----------------------------------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos)
{
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}
/* transform geodetic to ecef position -----------------------------------------
* transform geodetic position to ecef position
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *r        O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void pos2ecef(const double *pos, double *r)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    double e2=FE_WGS84*(2.0-FE_WGS84),v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
    
    r[0]=(v+pos[2])*cosp*cosl;
    r[1]=(v+pos[2])*cosp*sinl;
    r[2]=(v*(1.0-e2)+pos[2])*sinp;
}
/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order 
*-----------------------------------------------------------------------------*/
extern void xyz2enu(const double *pos, double *E)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    
    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}
/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
extern void ecef2enu(const double *pos, const double *r, double *e)
{
    double E[9];
    
    xyz2enu(pos,E);
    matmul("NN",3,1,3,1.0,E,r,0.0,e);
}
/* transform local vector to ecef coordinate -----------------------------------
* transform local tangental coordinate vector to ecef
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *e        I   vector in local tangental coordinate {e,n,u}
*          double *r        O   vector in ecef coordinate {x,y,z}
* return : none
*-----------------------------------------------------------------------------*/
extern void enu2ecef(const double *pos, const double *e, double *r)
{
    double E[9];
    
    xyz2enu(pos,E);
    matmul("TN",3,1,3,1.0,E,e,0.0,r);
}
/* transform covariance to local tangental coordinate --------------------------
* transform ecef covariance to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *P        I   covariance in ecef coordinate
*          double *Q        O   covariance in local tangental coordinate
* return : none
*-----------------------------------------------------------------------------*/
extern void covenu(const double *pos, const double *P, double *Q)
{
    double E[9],EP[9];
    
    xyz2enu(pos,E);
    matmul("NN",3,3,3,1.0,E,P,0.0,EP);
    matmul("NT",3,3,3,1.0,EP,E,0.0,Q);
}

/* 5_DEBUG FUNCTIONS ---------------------------------------------------------*/
static FILE *fp_dbg=NULL;		/* file pointer of debug */

extern void debugopen(const char *file)
{
	if (!(fp_dbg=fopen(file,"w")))
		fp_dbg=stderr;
}
extern void debugclose(void)
{
	if (fp_dbg&&fp_dbg!=stderr)
		fclose(fp_dbg);
	fp_dbg=NULL;
}
extern void writelog(const char *format, ...)
{
	va_list argp;

	if (!fp_dbg) return;
	va_start(argp,format); vfprintf(fp_dbg,format,argp); va_end(argp);
}

/* 6_POSITIONING FUNCTIONS ---------------------------------------------------*/

/* satellite frequency ---------------------------------------------------------
* get satellite frequency
* args   : int    sat       I   satellite number
*          int    frq       I   frequency index (0:L1,1:L2,2:L5/3,...)
* return : frequency value (0.0: error)
*-----------------------------------------------------------------------------*/
static double getfrq(int sat, int frq)
{
    int prn,sys;
    
	sys=satsys(sat,&prn);

    if (sys==SYS_GLO) {
        if      (frq==0) return 9.0;
		else if (frq==1) return 7.0;
	}
    else if (sys==SYS_BDS) {
        if      (frq==0) return FREQ1_BDS; /* B1 */
        else if (frq==1) return FREQ2_BDS; /* B2 */
        else if (frq==2) return FREQ3_BDS; /* B3 */
    }
    else {
        if      (frq==0) return FREQ1; /* L1/E1 */
        else if (frq==1) return FREQ2; /* L2 */
        else if (frq==2) return FREQ5; /* L5/E5a */
        else if (frq==3) return FREQ6; /* E6 */
        else if (frq==4) return FREQ7; /* E5b */
        else if (frq==5) return FREQ8; /* E5a+b */
    }
    return 0.0;
}

extern void fillfrq(nav_t *nav)
{
    int i,j;

    /* update frequency */
    for (i=0;i<MAXSAT;i++) 
		for (j=0;j<NFREQ;j++)
			nav->frq[i][j]=getfrq(i+1,j);
}
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args	 : double	*rs		I	satellite position (ecef at transmission) (m)
*		   double	*rr		I	receiver position  (ecef at reception) (m)
*		   double	*e		O	line of sight vector (ecef)
* return : geometric distance (m) (0>:err/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
	double r;
	int i;

	if (norm(rs,3)<RE_WGS84) return -1.0;
	for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
	r=norm(e,3);
	for (i=0;i<3;i++) e[i]/=r;
	return r+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}
/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle 
* args	 : double *pos		I	geodetic position {lat,lon,h} (rad,m)
*		   double *e		I	receiver-to-satellite unit vector (ecef)
*		   double *azel		IO	azimuth/elevation {az,el} (rad) (NULL: no output)
*								(0.0<=azel[0]<2*pi, -pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
extern  double satazel(const double *pos, const double *e, double *azel)
{
	double az=0.0,el=PI/2.0,enu[3];

	if (pos[2]>-RE_WGS84) {
		ecef2enu(pos,e,enu);
		az=dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
		if (az<0.0) az+=2*PI;
		el=asin(enu[2]);
	}
	if (azel) {azel[0]=az; azel[1]=el;};
	return el;
}
/* compute dops ----------------------------------------------------------------
* compute DOP (dilution of precision)
* args	 : int	  ns		I	number of satellites
*		   double *azel		I	satellite azimuth/elevation angle (rad)
*		   double elmin		I	elevation cutoff angle (rad)
*		   double *dop		O	DOPs  {GDOP,PDOP,HDOP,VDOP}
* return : none
* notes  : dop[0]-[3] return 0 in case of error in computation
*-----------------------------------------------------------------------------*/
#define SQRT(x)		((x)<0.0||(x)!=(x)?0.0:sqrt(x))

extern void dops(int ns, const double *azel, double elmin, double *dop)
{
	double H[4*MAXSAT],Q[16],cosel,sinel;
	int i,n;

	for(i=0;i<4;i++) dop[i]=0.0;
	for (i=n=0;i<ns&&i<MAXSAT;i++) {
		if (azel[1+i*2]<elmin || azel[1+i*2]<=0.0) continue;
		cosel=cos(azel[1+i*2]);
		sinel=sin(azel[1+i*2]);
		H[  4*n]=cosel*sin(azel[i*2]);
		H[1+4*n]=cosel*cos(azel[i*2]);
		H[2+4*n]=sinel;
		H[3+4*n++]=1.0;
	}
	if (n<4) return;

	matmul("NT",4,4,n,1.0,H,H,0.0,Q);
	if(!matinv(Q,4)) {
		dop[0]=SQRT(Q[0]+Q[5]+Q[10]+Q[15]); /* GDOP */
		dop[1]=SQRT(Q[0]+Q[5]+Q[10]);		/* PDOP */
		dop[2]=SQRT(Q[0]+Q[5]);				/* HDOP */
		dop[3]=SQRT(Q[10]);					/* VDOP */
	}
}

/* 7_ATMOSPHERE FUNCTIONS ----------------------------------------------------*/

/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args	 : gtime_t time		I	time 
*		   double *pos		I	receiver position {lat,lan,h} (rad,m)
*		   double *azel		I	azimuth/elevation angle {az,el} (rad)
*		   double humi		I	relative humidity
*		   double trpw		O	wet tropospheric delay (m)
* return : dry tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern double tropmodel(gtime_t time, const double *pos, double humi, double *wet)
{
	const double temp0=15.0; /* temparature at sea level (degC)*/
	double hgt,pres,temp,e,trph,trpw;

	if (pos[2]<-100.0||1E4<pos[2]) return 0.0;

	/* standard atmosphere */
	hgt=pos[2]<0.0?0.0:pos[2];

	pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
	temp=temp0-6.5E-3*hgt+273.15;
	e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

	/* saastamoinen model */
	trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3);
	trpw=0.0022768*(1255.0/temp+0.05)*e;

	if (wet) *wet=trpw;
	return trph;
}

static double interpc(const double coef[], double lat)
{
    int i=(int)(lat/15.0);
    if (i<1) return coef[0]; else if (i>4) return coef[4];
    return coef[i-1]*(1.0-lat/15.0+i)+coef[i]*(lat/15.0-i);
}
static double mapf(double el, double a, double b, double c)
{
    double sinel=sin(el);
    return (1.0+a/(1.0+b/(1.0+c)))/(sinel+(a/(sinel+b/(sinel+c))));
}
static double nmf(gtime_t time, const double pos[], const double azel[], double *mapfw)
{
    /* ref [4] table 3 */
    /* hydro-ave-a,b,c, hydro-amp-a,b,c; wet-a,b,c at latitude 15,30,45,60,75 */
    const double coef[][5]={
        { 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3}, /*a*/
        { 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3}, /*b*/
        { 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3}, /*c*/
        
        { 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5}, /*a*/
        { 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5}, /*b*/
        { 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5}, /*c*/
        
        { 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4}, /*a*/
        { 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3}, /*b*/
        { 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}  /*c*/
    };
    const double aht[]={ 2.53E-5, 5.49E-3, 1.14E-3}; /* height correction */
    
    double y,cosy,ah[3],aw[3],dm,el=azel[1],lat=pos[0]*R2D,hgt=pos[2];
    int i;
    
    if (el<=0.0) {
        if (mapfw) *mapfw=0.0;
        return 0.0;
    }
    /* year from doy 28, added half a year for southern latitudes */
    y=(time2doy(time)-28.0)/365.25+(lat<0.0?0.5:0.0);
    
    cosy=cos(2.0*PI*y);
    lat=fabs(lat);
    
    for (i=0;i<3;i++) {
        ah[i]=interpc(coef[i  ],lat)-interpc(coef[i+3],lat)*cosy;
        aw[i]=interpc(coef[i+6],lat);
    }
    /* ellipsoidal height is used instead of height above sea level */
    dm=(1.0/sin(el)-mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;
    
    if (mapfw) *mapfw=mapf(el,aw[0],aw[1],aw[2]);
    
    return mapf(el,ah[0],ah[1],ah[2])+dm;
}
/* troposphere mapping function ------------------------------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t t        I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *mapfw    IO  wet mapping function (NULL: not output)
* return : dry mapping function
* note   : see ref (4) (NMF)
*-----------------------------------------------------------------------------*/
extern double tropmapf(gtime_t time, const double pos[], const double azel[], double *mapfw)
{
    if (pos[2]<-1000.0||pos[2]>20000.0) {
        if (mapfw) *mapfw=0.0;
        return 0.0;
    }
    return nmf(time,pos,azel,mapfw); /* NMF */
}