/*------------------------------------------------------------------------------
* rinexobs.c : rinex observation file functions
*
*
* reference :
*     (1) W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 3.01, June 22, 2009
*     (2) RINEX The Receiver Independent Exchange Format Version 3.02,
*         International GNSS Service (IGS), RINEX Working Group and Radio
*         Technical Commission for Maritime Services Special Committee 104
*         (RTCM-SC104), December 10, 2012
*     (2) RINEX The Receiver Independent Exchange Format Version 3.03,
*         International GNSS Service (IGS), RINEX Working Group and Radio
*         Technical Commission for Maritime Services Special Committee 104
*         (RTCM-SC104), July 14, 2015
*		
* notes   :
*		1) These functions support RINEX 3.01, 3.02, 3.03 
*
* version : yyyy/mm/dd		ver		description
*			2018/04/04		1.0		create
*
* Copyright (C) 2018 by H.KARIMI, All rights reserved.
*-----------------------------------------------------------------------------*/
#include "gnsslib.h"

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)		((x)*(x))

#define NUMSYS		4					/* number of systems*/
#define MAXRNXLEN	(16*MAXOBSTYPE+4)	/* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */
#define NINCOBS     262144              /* incrimental number of obs data */

static const int navsys[]={				/* satellite systems */
	SYS_GPS,SYS_GLO,SYS_GAL,SYS_BDS,0
};
static const char syscodes[]="GREC";	/* satellite systems code */
static const char obscodes[]="CLDS";	/* obs type codes */
static const char frqcodes[]="125678";	/* frequency codes */

/* type definition -----------------------------------------------------------*/
typedef struct {						/* signal index type */
	int n;								/* number of index */
	int frq[MAXOBSTYPE];				/* signal frequency (1:L1, 2:L2) */
	int pos[MAXOBSTYPE];				/* signal index in obs data (-1:no) */
	unsigned char pri[MAXOBSTYPE];		/* signal priority (15-0) */
	unsigned char type[MAXOBSTYPE];		/* type (0:C,1:L,2:D,3:S) */
	unsigned char code[MAXOBSTYPE];		/* obs code (CODE_xxx)*/
} sigind_t;  

/* set string without tail space ---------------------------------------------*/
static void setstr(const char *dst, const char *src, int n)
{
	char *p=dst;
	const char *q=src;
	while (*q&&q<src+n) *p++=*q++;
	*p--='\0';
	while( p>=dst&&*p==' ') *p--='\0';
}
/* satellite to satellite code -----------------------------------------------*/
static int sat2code(int sat, char *code)
{
	int prn;
	switch (satsys(sat,&prn)) {
		case SYS_GPS: sprintf(code,"G%2d",prn-MINPRNGPS+1); break;
		case SYS_GLO: sprintf(code,"R%2d",prn-MINPRNGLO+1); break;
		case SYS_GAL: sprintf(code,"E%2d",prn-MINPRNGAL+1); break;
		case SYS_BDS: sprintf(code,"C%2d",prn-MINPRNBDS+1); break;
		default: return 0;
	}
	return 1;
}

/*------------------------------------------------------------------------------
* rinex observation file functions
*-----------------------------------------------------------------------------*/

/* decode obs header ---------------------------------------------------------*/
static void decode_obsh(FILE *fp, char *buff, double ver, int *tsys,
						char tobs[][MAXOBSTYPE][4], nav_t *nav)
{
	/* default codes for unknown code */
    const char *defcodes[]={
        "CWX    ",  /* GPS: L125___ */
        "CC     ",  /* GLO: L12____ */
        "X XXXX ",  /* GAL: L1_5678 */
        "X  XX  ",  /* BDS: L1__67_ */
    };
	int i,j,k,n,nt,prn,fcn;
	const char *p;
	char *label=buff+60;


    if     (strstr(label,"SYS / # / OBS TYPES" )) { /* ver.3 */
        if (!(p=strchr(syscodes,buff[0]))) {
            writelog("invalid system code: sys=%c\n",buff[0]);
            return;
        }
        i=(int)(p-syscodes);
        n=(int)str2num(buff,3,3);
        for (j=nt=0,k=7;j<n;j++,k+=4) {
            if (k>58) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                k=7;
            }
            if (nt<MAXOBSTYPE-1) setstr(tobs[i][nt++],buff+k,3);
        }
        *tobs[i][nt]='\0';
        
        /* change beidou B1 code: 3.02 draft -> 3.02/3.03 */
        if (i==3) {
            for (j=0;j<nt;j++) if (tobs[i][j][1]=='2') tobs[i][j][1]='1';
        }
        /* if unknown code in ver.3, set default code */
        for (j=0;j<nt;j++) {
            if (tobs[i][j][2]) continue;
            if (!(p=strchr(frqcodes,tobs[i][j][1]))) continue;
            tobs[i][j][2]=defcodes[i][(int)(p-frqcodes)];
            writelog("set default for unknown code: sys=%c code=%s\n",buff[0],tobs[i][j]);
        }
    }
    else if (strstr(label,"TIME OF FIRST OBS"   )) {
        if      (!strncmp(buff+48,"GPS",3)) *tsys=TSYS_GPS;
        else if (!strncmp(buff+48,"GLO",3)) *tsys=TSYS_UTC;
        else if (!strncmp(buff+48,"GAL",3)) *tsys=TSYS_GAL;
        else if (!strncmp(buff+48,"BDT",3)) *tsys=TSYS_BDS; /* ver.3.02 */
    }
    else if (strstr(label,"GLONASS SLOT / FRQ #")) {		/* ver.3.02 */
        if (nav) {
            for (i=0,p=buff+4;i<8;i++,p+=7) {
                if (sscanf(p,"R%2d %2d",&prn,&fcn)<2) continue;
                if (1<=prn&&prn<=MAXPRNGLO) nav->glo_fcn[prn-1]=fcn;
            }
        }
    }
}
/* read rinex header ---------------------------------------------------------*/
static int readrnxh(FILE *fp, double *ver, char *type, int *sys, int *tsys,
                    char tobs[][MAXOBSTYPE][4], nav_t *nav)
{
    char buff[MAXRNXLEN],*label=buff+60;
    int i=0;
    
    writelog("readrnxobs header\n");
    
    *ver=3.03; *type=' '; *sys=SYS_GPS;
    
    while (fgets(buff,MAXRNXLEN,fp)) {
        
        if (strlen(buff)<=60) continue;
        
        else if (strstr(label,"RINEX VERSION / TYPE")) {
            *ver=str2num(buff,0,9);
            *type=*(buff+20);
            
            /* satellite system */
            switch (*(buff+40)) {
                case ' ':
                case 'G': *sys=SYS_GPS;  *tsys=TSYS_GPS; break;
                case 'R': *sys=SYS_GLO;  *tsys=TSYS_UTC; break;
                case 'E': *sys=SYS_GAL;  *tsys=TSYS_GAL; break; 
                case 'C': *sys=SYS_BDS;  *tsys=TSYS_BDS; break; 
                case 'M': *sys=SYS_NON;  *tsys=TSYS_GPS; break; /* mixed */
                default :
                    writelog("not supported satellite system: %c\n",*(buff+40));
                    break;
            }
            continue;
        }

        /* file type */
        switch (*type) {
            case 'O': decode_obsh(fp,buff,*ver,tsys,tobs,nav); break;
        }
        if (strstr(label,"END OF HEADER")) return 1;
        
        if (++i>=MAXPOSHEAD) break; /* no rinex obs file */
    }
    return 0;
}
/* decode obs epoch ----------------------------------------------------------*/
static int decode_obsepoch(FILE *fp, char *buff, double ver, gtime_t *time,
                           int *flag, int *sats)
{
    int n;

    if (ver>2.99) { /* ver.3 */
        if ((n=(int)str2num(buff,32,3))<=0) return 0;
        
        *flag=(int)str2num(buff,31,1);
        
        if (buff[0]!='>'||str2time(buff,1,28,time)) {
            writelog("rinex obs invalid epoch: epoch=%29.29s\n",buff);
            return 0;
        }
    }
    return n;
}
/* decode obs data -----------------------------------------------------------*/
static int decode_obsdata(FILE *fp, char *buff, double ver, int mask,
                          sigind_t *index, obsd_t *obs)
{
    sigind_t *ind;
    double val[MAXOBSTYPE]={0};
    char satid[8]="";
    int i,j,n,m,stat=1,p[MAXOBSTYPE],k[16],l[16];

    if (ver>2.99) { /* ver.3 */
        strncpy(satid,buff,3);
        obs->sat=(unsigned char)satid2no(satid);
    }

    if (!obs->sat) {
        writelog("decode_obsdata: unsupported sat sat=%s\n",satid);
        stat=0;
    }
    else if (!(satsys(obs->sat,NULL)&mask)) {
        stat=0;
    }

    /* read obs data fields */
    switch (satsys(obs->sat,NULL)) {
        case SYS_GLO: ind=index+1; break;
        case SYS_GAL: ind=index+2; break;
        case SYS_BDS: ind=index+3; break;
        default:      ind=index  ; break;
    }
    for (i=0,j=3;i<ind->n;i++,j+=16) {
        if (stat) val[i]=str2num(buff,j,14);
    }

    if (!stat) return 0;
    
    for (i=0;i<NFREQ;i++) {
        obs->P[i]=0.0;
		obs->code[i]=0;
    }
    /* assign position in obs data */
    for (i=n=m=0;i<ind->n;i++) {
        
        p[i]=ind->pos[i];
    }
    
    /* save obs data */
    for (i=0;i<ind->n;i++) {
        if (p[i]<0||val[i]==0.0) continue;
        switch (ind->type[i]) {
            case 0: obs->P[p[i]]=val[i]; obs->code[p[i]]=ind->code[i]; break;
        }
    }
    return 1;
}
/* add obs data --------------------------------------------------------------*/
static int addobsdata(obs_t *obs, const obsd_t *data)
{
    obsd_t *obs_data;
    
    if (obs->nmax<=obs->n) {
        if (obs->nmax<=0) obs->nmax=NINCOBS; else obs->nmax*=2;
        if (!(obs_data=(obsd_t *)realloc(obs->data,sizeof(obsd_t)*obs->nmax))) {
            writelog("addobsdata: memalloc error n=%dx%d\n",sizeof(obsd_t),obs->nmax);
            free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
            return -1;
        }
        obs->data=obs_data;
    }
    obs->data[obs->n++]=*data;
    return 1;
}
/* set system mask -----------------------------------------------------------*/
static int set_sysmask(const char *opt)
{
    const char *p;
    int mask=SYS_NON;
    
    if (!(p=strstr(opt,"-sys="))) return SYS_ALL;
    
    for (p+=5;*p&&*p!=' ';p++) {
        switch (*p) {
            case 'G': mask|=SYS_GPS; break;
            case 'R': mask|=SYS_GLO; break;
            case 'E': mask|=SYS_GAL; break;
            case 'C': mask|=SYS_BDS; break;
        }
    }
    return mask;
}
/* set signal index ----------------------------------------------------------*/
static void set_index(double ver, int sys, char tobs[MAXOBSTYPE][4], sigind_t *ind)
{
    const char *p;
    char str[8],*optstr="";
    int i,j,k,n;
    
    for (i=n=0;*tobs[i];i++,n++) {
        ind->code[i]=obs2code(tobs[i]+1,ind->frq+i);
        ind->type[i]=(p=strchr(obscodes,tobs[i][0]))?(int)(p-obscodes):0;
        ind->pri[i] =getcodepri(sys,ind->code[i],NULL);
        ind->pos[i]=-1;
        
        /* frequency index for beidou */
        if (sys==SYS_BDS) {
            if      (ind->frq[i]==5) ind->frq[i]=2; /* B2 */
            else if (ind->frq[i]==4) ind->frq[i]=3; /* B3 */
        }
    }

    /* assign index for highest priority code */
    for (i=0;i<NFREQ;i++) {
        for (j=0,k=-1;j<n;j++) {
            if (ind->frq[j]==i+1&&ind->pri[j]&&(k<0||ind->pri[j]>ind->pri[k])) {
                k=j;
            }
        }
        if (k<0) continue;
        
        for (j=0;j<n;j++) {
            if (ind->code[j]==ind->code[k]) ind->pos[j]=i;
        }
    }
    ind->n=n;
}
/* read rinex obs data body --------------------------------------------------*/
static int readrnxobsb(FILE *fp, const char *opt, double ver, int *tsys,
                       char tobs[][MAXOBSTYPE][4], int *flag, obsd_t *data)
{
    gtime_t time={0};
    sigind_t index[4]={{0}};
    char buff[MAXRNXLEN];
    int i=0,n=0,nsat=0,sats[MAXOBS]={0},mask;
    
    /* set system mask */
    mask=set_sysmask(opt);
    
    /* set signal index */
    set_index(ver,SYS_GPS,tobs[0],index  );
    set_index(ver,SYS_GLO,tobs[1],index+1);
    set_index(ver,SYS_GAL,tobs[2],index+2);
    set_index(ver,SYS_BDS,tobs[3],index+3);
    
    /* read record */
    while (fgets(buff,MAXRNXLEN,fp)) {
        
        /* decode obs epoch */
        if (i==0) {
            if ((nsat=decode_obsepoch(fp,buff,ver,&time,flag,sats))<=0) {
                continue;
            }
        }
        else if (*flag<=1) {
            
            data[n].time=time;
            data[n].sat=(unsigned char)sats[i-1];
            
            /* decode obs data */
            if (decode_obsdata(fp,buff,ver,mask,index,data+n)&&n<MAXOBS) n++;
        }

        if (++i>nsat) return n;
    }
    return -1;
}
/* read rinex obs ------------------------------------------------------------*/
static int readrnxobs(FILE *fp, gtime_t ts, gtime_t te, double tint,
                      const char *opt, double ver, int *tsys,
                      char tobs[][MAXOBSTYPE][4], obs_t *obs)
{
    obsd_t *data;
    int i,n,flag=0,stat=0; /* epoch flag */
    
    writelog("readrnxobs body: ver=%.2f tsys=%d\n",ver,*tsys);
    
    if (!obs) return 0;
    
    if (!(data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))) return 0;
    
    /* read rinex obs data body (for an epoch each time) */
    while ((n=readrnxobsb(fp,opt,ver,tsys,tobs,&flag,data))>=0&&stat>=0) {
        
        for (i=0;i<n;i++) {
            
            /* utc -> gpst */
            if (*tsys==TSYS_UTC) data[i].time=utc2gpst(data[i].time);
        }
        /* screen data by time */
        if (n>0&&!screent(data[0].time,ts,te,tint)) continue;
        
        for (i=0;i<n;i++) {
            
            /* save obs data */
            if ((stat=addobsdata(obs,data+i))<0) break;
        }
    }
    writelog("readrnxobs data body: nobs=%d stat=%d\n",obs->n,stat);
    free(data);
    return stat;
}
/* read rinex file -----------------------------------------------------------*/
static int readrnxfp(FILE *fp, gtime_t ts, gtime_t te, double tint,
                const char *opt, char *type, obs_t *obs, nav_t *nav)
{
    double ver;
    int sys,tsys=TSYS_GPS;
    char tobs[NUMSYS][MAXOBSTYPE][4]={{""}};
    
    /* read rinex header */
    if (!readrnxh(fp,&ver,type,&sys,&tsys,tobs,nav)) return 0;
    
    /* read rinex body */
    switch (*type) {
        case 'O': return readrnxobs(fp,ts,te,tint,opt,ver,&tsys,tobs,obs);
    }
    writelog("unsupported rinex type ver=%.2f type=%c\n",ver,*type);
    return 0;
}
/* read rinex file -----------------------------------------------------------*/
static int readrnxfile(const char *file, gtime_t ts, gtime_t te, double tint,
                const char *opt, char *type, obs_t *obs, nav_t *nav)
{
    FILE *fp;
    int stat;

    if (!(fp=fopen(file,"r"))) {
        showmsg("rinex file open error...");
        return 0;
    }
    /* read rinex file */
    stat=readrnxfp(fp,ts,te,tint,opt,type,obs,nav);
    fclose(fp);
    return stat;
}
/* read rinex obs  files -------------------------------------------------------
* read rinex obs and nav files
* args   : char *file    I      file ("": stdin)
*         (gtime_t ts)   I      observation time start (ts.time==0: no limit)
*         (gtime_t te)   I      observation time end   (te.time==0: no limit)
*         (double tint)  I      observation time interval (s) (0:all)
*          char  *opt    I      rinex options (see below,"": no option)
*          obs_t *obs    IO     observation data   (NULL: no input)
*          nav_t *nav    IO     navigation data    (NULL: no input)
*          sta_t *sta    IO     station parameters (NULL: no input)
* return : status (1:ok,0:no data,-1:error)
* notes  : read data are appended to obs and nav struct
*          before calling the function, obs and nav should be initialized.
*          observation data and navigation data are not sorted.
*
*          read rinex options :
*            
*            -SYS=sys[,sys...]: select navigation systems
*                               (sys=G:GPS,R:GLO,E:GAL,C:BDS)
*
*-----------------------------------------------------------------------------*/
extern int readrnxt(const char *file, gtime_t ts, gtime_t te, double tint,
	                const char *opt, obs_t *obs, nav_t *nav)
{
    int stat=0;
    char type=' ';

	/* read from stin */
    if (!*file) {
        return readrnxfp(stdin,ts,te,tint,opt,&type,obs,nav);
    }

	stat=readrnxfile(file,ts,te,tint,opt,&type,obs,nav);

    return stat;
}
/* read whole rinex file without time span and interval ----------------------*/
extern int readrnx(const char *file, const char *opt, obs_t *obs, nav_t *nav)
{
    gtime_t t={0};
    
    return readrnxt(file,t,t,0.0,opt,obs,nav);
}