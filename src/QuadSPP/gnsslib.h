/*------------------------------------------------------------------------------
* gnsslib.h : program constants, types and function prototypes
* 
* version : yyyy/mm/dd		ver		description
*			2018/03/22		1.0		create
*			2018/04/26		1.1		modify type definition
*			2018/11/26		1.2		modify constants
*
* Copyright (C) 2018 by H.KARIMI, All rights reserved.
*-----------------------------------------------------------------------------*/
#ifndef GNSSLIB_H
#define GNSSLIB_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>


#ifdef __cplusplus
extern "C" {
#endif

/* constants -----------------------------------------------------------------*/
#define PI			3.1415926535897932	/* pi */
#define D2R			(PI/180.0)			/* deg to rad */
#define R2D			(180.0/PI)			/* rad to deg */
#define CLIGHT		299792458.0			/* speed of light (m/s) */

#define OMGE		7.2921151467E-5		/* earth angular velocity (rad/s) */
#define RE_WGS84	6378137.0			/* earth semimajor axis (WGS84) (m) */
#define FE_WGS84	(1.0/298.257223563) /* earth flattening (WGS84) */

#define MAXFREQ		6					/* max NFREQ */
#define FREQ1		1.57542E9			/* L1/E1  frequency (Hz) */
#define FREQ2		1.22760E9			/* L2     frequency (Hz) */
#define FREQ5		1.17645E9			/* L5/E5a frequency (Hz) */
#define FREQ6		1.27875E9			/* E6     frequency (Hz) */
#define FREQ7		1.20714E9			/* E5b    frequency (Hz) */
#define FREQ8		1.191795E9			/* E5a+b  frequency (Hz) */
#define FREQ1_GLO   1.60200E9			/* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO	0.56250E6			/* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO	1.24600E9			/* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO	0.43750E6			/* GLONASS G2 bias frequency (Hz/n) */
#define FREQ1_BDS	1.561098E9			/* BeiDou B1 frequency (Hz) */
#define FREQ2_BDS	1.20714E9			/* BeiDou B2 frequency (Hz) */
#define FREQ3_BDS	1.26852E9			/* BeiDou B3 frequency (Hz) */

#define SYS_NON		0x00				/* system: none */
#define SYS_GPS		0x01				/* system: GPS */
#define SYS_GLO		0x02				/* system: GLONASS */
#define SYS_GAL		0x04				/* system: Galileo */
#define SYS_BDS		0x08				/* system: BeiDou */
#define SYS_ALL		0xFF				/* system: all */

#define TSYS_GPS	0					/* time system: GPS time */
#define TSYS_UTC	1					/* time system: UTC */
#define TSYS_GLO	2					/* time system: GLONASS time */
#define TSYS_GAL	3					/* time system: Galileo time */
#define TSYS_BDS	4					/* time system: BeiDou time */

#define NFREQ       3                   /* number of carrier frequencies */
#define NFREQGLO    2                   /* number of carrier frequencies of GLONASS */

#define MINPRNGPS	1					/* min satellite PRN number of GPS */
#define MAXPRNGPS	32					/* max satellite PRN number of GPS */
#define NSATGPS		(MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS		1

#define MINPRNGLO	1					/* min satellite slot number of GLONASS */
#define MAXPRNGLO	27					/* max satellite slot number of GLONASS */
#define NSATGLO		(MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO		1

#define MINPRNGAL	1					/* min satellite PRN number of Galileo */
#define MAXPRNGAL	36					/* max satellite PRN number of Galileo */
#define NSATGAL		(MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL		1

#define MINPRNBDS	1					/* min satellite PRN number of BeiDou */
#define MAXPRNBDS	35					/* max satellite PRN number of BeiDou */
#define NSATBDS		(MAXPRNBDS-MINPRNBDS+1) /* number of BeiDou satellites */
#define NSYSBDS		1

#define NSYS		(NSYSGPS+NSYSGLO+NSYSGAL+NSYSBDS)/* number of systems */
#define MAXSAT		(NSATGPS+NSATGLO+NSATGAL+NSATBDS)/* max satellite number (1 to MAXSAT) */
#define MAXOBS		64					/* max number of obs in an epoch */
#define MAXOBSTYPE	64					/* max number of obs type in RINEX */
#define MAXLEAPS    64                  /* max number of leap seconds table */
#define MAXSOLMSG   8191                /* max length of solution message */
#define DTTOL       0.005               /* tolerance of time difference (s) */

#define SOLF_LLH    0                   /* solution format: lat/lon/height */
#define SOLF_XYZ    1                   /* solution format: x/y/z-ecef */

#define CODE_NON	0					/* obs code: none or unknown */
#define CODE_L1C	1					/* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL) */
#define CODE_L1P	2					/* obs code: L1P,G1P	(GPS,GLO) */
#define CODE_L1W	3					/* obs code: L1 Z-track (GPS) */
#define CODE_L1Y	4					/* obs code: L1Y		(GPS) */
#define CODE_L1M	5					/* obs code: L1M		(GPS) */
#define CODE_L1N	6					/* obs code: L1codeless (GPS) */
#define CODE_L1S	7					/* obs code: L1C(D)     (GPS) */
#define CODE_L1L	8					/* obs code: L1C(P)     (GPS) */
#define CODE_L1E	9					/* (not used) */
#define CODE_L1A	10					/* obs code: E1A		(GAL) */
#define CODE_L1B	11					/* obs code: E1B		(GAL) */
#define CODE_L1X	12					/* obs code: L1C(D+P),E1B+C,B1I+Q (GPS,GAL,BDS) */
#define CODE_L1Z	13					/* obs code: E1A+B+C    (GAL) */
#define CODE_L2C	14					/* obs code: L2C/A,G2C/A   (GPS,GLO) */
#define CODE_L2D	15					/* obs code: Semi-codeless (GPS) */
#define CODE_L2S	16					/* obs code: L2C(M)		(GPS) */
#define CODE_L2L	17					/* obs code: L2C(L)		(GPS) */
#define CODE_L2X	18					/* obs code: L2C(M+L),B1I+Q (GPS,BDS) */
#define CODE_L2P	19					/* obs code: L2P,G2P	(GPS,GLO) */
#define CODE_L2W	20					/* obs code: L2 Z-track (GPS) */
#define CODE_L2Y	21					/* obs code: L2Y		(GPS) */
#define CODE_L2M	22					/* obs code: L2M		(GPS) */
#define CODE_L2N	23					/* obs code: L2codelese (GPS) */
#define CODE_L5I	24					/* obs code: L5/E5aI	(GPS,GAL) */
#define CODE_L5Q	25					/* obs code: L5/E5aQ	(GPS,GAL) */
#define CODE_L5X	26					/* obs code: L5/E5aI+Q	(GPS,GAL) */
#define CODE_L7I	27					/* obs code: E5bI,B2I	(GAL,BDS) */
#define CODE_L7Q	28					/* obs code: E5bQ,B2Q	(GAL,BDS) */
#define CODE_L7X	29					/* obs code: E5bI+Q,B2I+Q (GAL,BDS) */
#define CODE_L6A    30                  /* obs code: E6A        (GAL) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,B3I+Q(GAL,BDS) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C    (GAL) */
#define CODE_L8I	35					/* obs code: E5(a+b)I	(GAL) */	
#define CODE_L8Q	36					/* obs code: E5(a+b)Q	(GAL) */	
#define CODE_L8X	37					/* obs code: E5(a+b)I+Q	(GAL) */
#define CODE_L2I	38					/* obs code: B1I		(BDS) */
#define CODE_L2Q	39					/* obs code: B1Q		(BDS) */
#define CODE_L6I	40					/* obs code: B3I		(BDS) */
#define CODE_L6Q	41					/* obs code: B3Q		(BDS) */
#define CODE_L3I    42                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    43                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    44                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I	45					/* obs code: B1I		(BDS) */
#define CODE_L1Q	46					/* obs code: B1Q		(BDS) */
#define MAXCODE		46					/* max number of obs codes */

#define COMMENTH    "#"                 /* comment line indicator for solution */

/* type definitions ----------------------------------------------------------*/

typedef struct {				/* time struct */
	time_t time;				/* time (s) expressed by standard time_t */
	double sec;					/* fraction of second under 1s */
} gtime_t;

typedef struct {				/* observation data record */
	gtime_t time;				/* receiver sampling time (GPST) */
	unsigned char sat;			/* satellite number */
	unsigned char code[NFREQ];	/* code indicator (CODE_xxx) */
	double P[NFREQ];			/* observation data pseudorange (m) */
} obsd_t;

typedef struct {				/* observation data */
	int n, nmax;				/* number of observation data/allocated */
	obsd_t *data;				/* observation data record */
} obs_t;

typedef struct {				/* precise ephemeris type */
	gtime_t time;				/* time (GPST) */
	double pos[MAXSAT][4];		/* satellite position/clock (ecef) (m|s) */
	float  std[MAXSAT][4];		/* satellite position/clock std (m|s) */
} peph_t;

typedef struct {				 /* navigation data type*/
	int n, nmax;				 /* number of precise ephemeris */
	peph_t *peph;				 /* precise ephemeris */
	int glo_fcn[MAXPRNGLO+1];    /* glonass frequency channel number */
	double frq[MAXSAT][NFREQ];   /* frequency value */
} nav_t;

typedef struct {				/* solution type */
    gtime_t time;				/* time (GPST) */
    double rr[3];				/* position {x,y,x} (m) */
	double dtr[4];				/* receiver clock bias to time systems (s) */
    float  qr[6];				/* position variance/covariance (m^2) */
								/* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx}    */
	unsigned char ns;			/* number of all   satellites */
    unsigned char nv;			/* number of valid satellites */
} sol_t;

typedef struct {		/* processing option type */
	int navsys;         /* navigation system */
	int posf;           /* solution format (SOLF_xxx) */
	double elmin;		/* elevation mask angle (rad) */
	double maxgdop;     /* reject threshold of gdop */
	char rnxopt[256];	/* rinex option */
	char sep[64];       /* field separator */
} opt_t;

/* functions prototype -------------------------------------------------------*/
/* satellite, systems, code functions --------------------------------------- */
extern int  satno (int sys, int prn);
extern int  satsys(int sat, int *prn);
extern int  satid2no(const char *id);
extern void satno2id(int sat, char *id);
extern unsigned char obs2code(const char *obs, int *freq);
extern char *code2obs(unsigned char code, int *freq);

/* matrix and vector functions -----------------------------------------------*/
extern double *mat  (int n, int m);
extern int    *imat (int n, int m);
extern double *zeros(int n, int m);

extern double dot (const double *a,const double *b, int n);
extern double norm(const double *a, int n);

extern void matcpy(double *A, const double *B, int n, int m);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
				   const double *A, const double *B, double beta, double *C);
extern int  matinv(double *A, int n);
extern int  lsq   (const double *A, const double *y, int n, int m, double *x, double *Q);

/* time and string functoins -------------------------------------------------*/
extern double  str2num (const char *s, int i, int n);
extern int     str2time(const char *s, int i, int n, gtime_t *t);
extern void    time2str(gtime_t t, char *str,int n);

extern gtime_t epoch2time(const double *ep);
extern void    time2epoch(gtime_t, double *ep);

extern gtime_t gpst2time(int week, double sec);
extern double  time2gpst(gtime_t t, int *week);

extern gtime_t  timeadd (gtime_t t, double sec);
extern double	timediff(gtime_t t1, gtime_t t2);

extern gtime_t gpst2utc(gtime_t t);
extern gtime_t utc2gpst(gtime_t t);

extern double time2doy(gtime_t t);

extern int screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
extern char *time_str(gtime_t t, int n);

/* coordinates functions -----------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos);
extern void pos2ecef(const double *pos, double *r);

extern void xyz2enu (const double *pos, double *E);

extern void ecef2enu(const double *pos, const double *r, double *e);
extern void enu2ecef(const double *pos, const double *e, double *r);

extern void covenu(const double *pos, const double *P, double *Q);
/* debug functions -----------------------------------------------------------*/
extern void debugopen (const char *file);
extern void debugclose(void);
extern void writelog  (const char *format, ...);

/* positioning functions -----------------------------------------------------*/
extern void   fillfrq(nav_t *nav);
extern double satazel(const double *pos, const double *e, double *azel);
extern double geodist(const double *rs, const double *rr, double *e);
extern void   dops   (int ns, const double *azel, double elmin, double *dop);
 
/* atmosphere models functions -----------------------------------------------*/
extern double tropmodel(gtime_t time, const double *pos, double humi, double *wet);
extern double tropmapf (gtime_t time, const double *pos, const double *azel, double *mapfw);

/* rinex functions -----------------------------------------------------------*/
extern int readrnx (const char *file, const char *opt, obs_t *obs, nav_t *nav);
extern int readrnxt(const char *file, gtime_t ts, gtime_t te,double tint,
					const char *opt, obs_t *obs, nav_t *nav);

/* precise ephemeris functions -----------------------------------------------*/
extern int readsp3 (const char *file, nav_t *nav);
extern void satposs (const obsd_t *obs, int n, const nav_t *nav,
				     double *rs, double *dts, double *var);

/* standard positioning ------------------------------------------------------*/
extern int pntpos (const obsd_t *obs, int n, const nav_t *nav, const opt_t *opt, sol_t *sol);
extern int postpos(gtime_t ts, gtime_t te, double ti, const opt_t *popt, char **infile);
extern int showmsg(char *format, ...);

#ifdef __cplusplus
}
#endif
#endif /*GNSSLIB_H*/