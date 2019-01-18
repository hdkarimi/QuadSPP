/*------------------------------------------------------------------------------
* postpro.c : post-processing posditioning
*
* version : yyyy/mm/dd		ver		description
*			2018/05/10		1.0		create
*
* Copyright (C) 2018 by H.KARIMI, All rights reserved.
*-----------------------------------------------------------------------------*/
#include "gnsslib.h"

/* constants/global variables ------------------------------------------------*/
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define MAXINFILE   2			/* max number of input files */

static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static int   iobs=0;			/* current reference observation data index */

/* option for seperation -----------------------------------------------------*/
static const char *opt2sep(const opt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}
/* output solution as the form of x/y/z-ecef ---------------------------------*/
static int outecef(unsigned char *buff, const char *s, const sol_t *sol,opt_t *opt)
{
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;
    
    p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f",
               s,sep,sol->rr[0],sep,sol->rr[1],sep,sol->rr[2],sep,sol->ns,sep,sol->nv,sep,
			   SQRT(sol->qr[0]),sep,SQRT(sol->qr[1]),sep,SQRT(sol->qr[2]));
    
    p+=sprintf(p,"\n");
    return p-(char *)buff;
}
/* solution to covariance ----------------------------------------------------*/
static void soltocov(const sol_t *sol, double *P)
{
    P[0]     =sol->qr[0]; /* xx or ee */
    P[4]     =sol->qr[1]; /* yy or nn */
    P[8]     =sol->qr[2]; /* zz or uu */
    P[1]=P[3]=sol->qr[3]; /* xy or en */
    P[5]=P[7]=sol->qr[4]; /* yz or nu */
    P[2]=P[6]=sol->qr[5]; /* zx or ue */
}
/* output solution as the form of lat/lon/height -----------------------------*/
static int outpos(unsigned char *buff, const char *s, const sol_t *sol, opt_t *opt)
{
    double pos[3],P[9],Q[9];
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;

    ecef2pos(sol->rr,pos);
	soltocov(sol,P);
    covenu(pos,P,Q);

	p+=sprintf(p,"%s%s%14.9f%s%14.9f",s,sep,pos[0]*R2D,sep,pos[1]*R2D);

    p+=sprintf(p,"%s%10.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f",
               sep,pos[2],sep,sol->ns,sep,sol->nv,sep,
			   SQRT(Q[4]),sep,SQRT(Q[0]),sep,SQRT(Q[8]));
    
    p+=sprintf(p,"\n");
    return p-(char *)buff;
}
/* output processing options -------------------------------------------------*/
static int outprcopts(unsigned char *buff, const opt_t *opt)
{
    const int sys[]={SYS_GPS,SYS_GLO,SYS_GAL,SYS_BDS,0};
    const char s1[]={"iono-free"};
    const char s2[]={"saastamoinen"};
	const char s3[]={"niell"};
    const char s4[]={"precise"};
    const char *s5[]={"gps","glonass","galileo","beidou",""};
    int i;
    char *p=(char *)buff;
    
    p+=sprintf(p,"%s elev mask : %.1f deg\n",COMMENTH,opt->elmin*R2D);
	p+=sprintf(p,"%s ionos corr: %s\n",COMMENTH,s1);
    p+=sprintf(p,"%s tropo corr: %s\n",COMMENTH,s2);
	p+=sprintf(p,"%s tropo mapf: %s\n",COMMENTH,s3);
    p+=sprintf(p,"%s ephemeris : %s\n",COMMENTH,s4);
    if (opt->navsys!=SYS_NON) {
        p+=sprintf(p,"%s navi sys  :",COMMENTH);
        for (i=0;sys[i];i++) {
            if (opt->navsys&sys[i]) p+=sprintf(p," %s",s5[i]);
        }
        p+=sprintf(p,"\n");
    }
	p+=sprintf(p,"%s \n",COMMENTH);
    return p-(char *)buff;
}
/* output solution header ----------------------------------------------------*/
static int outsolheads(unsigned char *buff, const opt_t *opt)
{
    const char *s[]={"GPST","UTC "},*sep=opt2sep(opt);
    char *p=(char *)buff;
    
    p+=sprintf(p,"%s  %-*s%s",COMMENTH,(15),s[0],sep);

    if      (opt->posf==SOLF_LLH) { /* lat/lon/hgt */
		p+=sprintf(p,"%14s%s%14s%s%10s%s%3s%s%3s%s%8s%s%8s%s%8s",
			       "latitude(deg)",sep,"longitude(deg)",sep,"height(m)",sep,
			       "ns",sep,"nv",sep,"sdn(m)",sep,"sde(m)",sep,"sdu(m)");
	}
    else if (opt->posf==SOLF_XYZ) { /* x/y/z-ecef */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s",
			       "x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,
			       "ns",sep,"nv",sep,"sdx(m)",sep,"sdy(m)",sep,"sdz(m)");
    }
    p+=sprintf(p,"\n");
    return p-(char *)buff;
}
/* output solution body ------------------------------------------------------*/
static int outsols(unsigned char *buff,const sol_t *sol,const opt_t *opt)
{
    gtime_t time;
    char s[64];
    unsigned char *p=buff;

    time=sol->time;
	time2str(time,s,0);
    switch (opt->posf) {
        case SOLF_LLH:  p+=outpos (p,s,sol,opt);   break;
        case SOLF_XYZ:  p+=outecef(p,s,sol,opt);   break;
    }
    return p-buff;
}
/* output processing option --------------------------------------------------*/
static void outprcopt(FILE *fp, const opt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;

    if ((n=outprcopts(buff,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output solution header ----------------------------------------------------*/
static void outsolhead(FILE *fp, const opt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;

    if ((n=outsolheads(buff,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output solution body ------------------------------------------------------*/
static void outsol(FILE *fp, const sol_t *sol, opt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;

    if ((n=outsols(buff,sol,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, char **file, int n, const opt_t *popt)
{
    const char s1[]="GPST";
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];

	fprintf(fp,"%s program   : Single Point Positioning\n",COMMENTH);
	for (i=0;i<n;i++){
		fprintf(fp,"%s inp file  : %s\n",COMMENTH,file[i]);
	}
	ts=obss.data[0].time;
	te=obss.data[obss.n-1].time;
	t1=time2gpst(ts, &w1);
	t2=time2gpst(te, &w2);
	time2str(ts,s2,1);
	time2str(te,s3,1);
	fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n", COMMENTH,s2,s1,w1,t1);
	fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n", COMMENTH,s3,s1,w2,t2);

	outprcopt (fp, popt);
	outsolhead(fp, popt);

}
/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, const opt_t *popt, char **infile)
{
    FILE *fp=stdout;
    
    if (*outfile) {
        if (!(fp=fopen(outfile,"w"))) {
            showmsg("error : open output file %s",outfile);
            return 0;
        }
    }
    /* output header */
    outheader(fp,infile,2,popt);
    
    if (*outfile) fclose(fp); 
    return 1;
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    return !*outfile?stdout:fopen(outfile,"a");
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i)
{
    double tt;
    int n;
    
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (tt>DTTOL) break;
    }
    return n;
}
/* input obs data ------------------------------------------------------------*/
static int inputobs(obsd_t *obs)
{
    gtime_t time={0};
    int i,n;
    
    if (0<=iobs&&iobs<obss.n) {
        time=obss.data[iobs].time;
        showmsg("processing : %s ",time_str(time,0));
		
		n=nextobsf(&obss,&iobs);

		for (i=0;i<n&&n<MAXOBS;i++) obs[i]=obss.data[iobs+i];
		iobs+=n;
        
		return n;
    }
	return 0;
}
/* process positioning -------------------------------------------------------*/
static void procpos(FILE *fp, const opt_t *popt)
{
    sol_t  sol={{0}};
    obsd_t obs[MAXOBS]; /* for each epoch */
    int nobs;

    while ((nobs=inputobs(obs))>0) {

		if (!pntpos(obs,nobs,&navs,popt,&sol)) {
			showmsg("rejected !");
			continue;
		}
		outsol(fp,&sol,popt);
	}
}
/* execute processing session ------------------------------------------------*/
static int execut(gtime_t ts,gtime_t te,double ti,const opt_t *opt,char **infile,char *oufile)
{
	FILE *fp;
	opt_t popt=*opt;

	/* write header to output file */
	if (!outhead(oufile,&popt,infile)) return -1;

	if ((fp=openfile(oufile))) {
		procpos(fp,&popt); 
		fclose(fp);
		return 0;
	}
	return -1; /* err */
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          opt_t *opt	    I   processing options
*          char   **infile  I   input files (see below)
* return : status (0:ok,0>:error)
* notes  : input files should contain observation data, precise ephemeris

*          the type of an input file is recognized by the file extention as 
*          follows:
*			   .rnx,.RNX,.yyo,yyO,yy.obs: rinex observation file
*              .sp3,.SP3,.eph*,.EPH*    : precise ephemeris (sp3)
*-----------------------------------------------------------------------------*/
extern int postpos(gtime_t ts,gtime_t te,double ti, const opt_t *popt,char **inpfile)
{
    int p,stat=0;
    char out1[]="solution.xyz",out2[]="solution.xyz.debug",*ext;

	/* open debug trace file */
	debugopen(out2);

    if (ts.time!=0&&te.time!=0) {
        if (timediff(te,ts)<0.0) {
            showmsg("error : no period");
            return 0;
		}
	}
	showmsg("reading...");
	/* read rinex observation file */
	if (!readrnxt(inpfile[0],ts,te,ti,popt->rnxopt,&obss,&navs)){
		showmsg("rinex file reading not complete...\n");
		return -1;
	}

	/* get frequency value */
	fillfrq(&navs); 

	/* read precise ephemeris file */
	if (!readsp3(inpfile[1],&navs)) {
		showmsg("sp3 file open error...\n");
		return -1;
	}

	/* execute processing session */
	stat=execut(ts,te,ti,popt,inpfile,out1);
	
	/* close debug trace file */
	debugclose();

	return stat;
}