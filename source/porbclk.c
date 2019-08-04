/*------------------------------------------------------------------------------
* porbclk.c : precise orbit and clock functions
*
* references :
*     (1) S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c), 2010
*     (1) S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-d), 2016
*
* version : yyyy/mm/dd		ver		description
*			2018/04/16		1.0		create
*
* Copyright (C) 2018 by H.KARIMI, All rights reserved.
*-----------------------------------------------------------------------------*/
#include "gnsslib.h"

#define SQR(x)		((x)*(x))

#define NMAX		10				/* order of polynomial interpolation */
#define MAXDTE		900				/* max time difference to ephem time (s) */

/* satellite code to satellite system ----------------------------------------*/
static int code2sys(char code)
{
	if (code=='G') return SYS_GPS;
	if (code=='R') return SYS_GLO;
	if (code=='E') return SYS_GAL;
	if (code=='C') return SYS_BDS;
	return SYS_NON;
}
/* read sp3-c header ---------------------------------------------------------*/
static int readsp3ch(FILE *fp, gtime_t *time, double *bfact, double *eint, 
					 char *ftype ,char *tsys)
{
    int i,ns=0;
    char buff[1024];
    
    writelog("read sp3-c header: \n");
    
    for (i=0;i<22;i++) {
        if (!fgets(buff,sizeof(buff),fp)) break;
        
        if (i==0) {
			*ftype = buff[0];
            if (str2time(buff,1,28,time)) return 0;
        }
		else if (i==1) {
			*eint=str2num(buff,24,14);
		}
        else if (i==2) {
			ns=(int)str2num(buff,4,2);
        }
        else if (i==12) {
            strncpy(tsys,buff+9,3); tsys[3]='\0';
        }
        else if (i==14) {
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
        }
    } /* reads until begining of comments */
    return ns;
}
/* read sp3-d header ---------------------------------------------------------*/
static int readsp3dh(FILE *fp, gtime_t *time, double *bfact, double *eint,
					 char *ftype, char *tsys)
{
    int i,ns=0,nl;
    char buff[1024];
    
    writelog("read sp3-d header: \n");

	for(i=1;i<=3;i++) {
		if (!fgets(buff,sizeof(buff),fp)) break;

		if (i==1) {
			*ftype=buff[0];
			if (str2time(buff,1,28,time)) return 0;
		}
		else if (i==2) {
			*eint=str2num(buff,24,14);
		}
		else if (i==3) { 
			ns=(int)str2num(buff,3,3);
		}
	}

	nl=(ns/17+1); /* number of line start with +/++ */
	for (i=4;i<=2+2*nl+6;i++) {
		if (!fgets(buff,sizeof(buff),fp)) break;

		if (i==2+2*nl+1) {
			strncpy(tsys,buff+9,3); tsys[3]='\0';
		}
		else if (i==2+2*nl+3) {
			bfact[0]=str2num(buff, 3,10);
			bfact[1]=str2num(buff,14,12);
		}
	} /* reads until begining of comments */
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;
    
    if (nav->n>=nav->nmax) {
        nav->nmax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nmax))) {
            writelog("readsp3b malloc error n=%d\n",nav->nmax);
            free(nav->peph); nav->peph=NULL; nav->n=nav->nmax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->n++]=*peph;
    return 1;
}
/* read sp3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, int ns, double *bfact, char ftype, char *tsys, nav_t *nav)
{
    peph_t peph;
    gtime_t time;
    double val,std,base;
    int i,j,sat,sys,prn,n=ns*(ftype=='P'?1:2),v;
    char buff[1024];
    
    writelog("read sp3 body:\n");
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (!strncmp(buff,"EOF",3)) break;
        
        if (buff[0]!='*'||str2time(buff,3,28,&time)) {
            writelog("sp3 invalid epoch %31.31s\n",buff);
            continue;
        }
        if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
        peph.time =time;
        
        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<4;j++) {
                peph.pos[i][j]=0.0;
                peph.std[i][j]=0.0f;
            }
        }
        for (i=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {
            
            if (strlen(buff)<4||(buff[0]!='P')) continue;
            
            sys=code2sys(buff[1]);
            prn=(int)str2num(buff,2,2);
            
            if (!(sat=satno(sys,prn))) continue;

            for (j=0;j<4;j++) {

                val=str2num(buff, 4+j*14,14);
                std=str2num(buff,61+j* 3,j<3?2:3);
                
                if (buff[0]=='P') { /* position */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                        v=1; /* valid epoch */
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                    }
                }
            }
        } /* end of for */

        if (v) {
            if (!addpeph(nav,&peph)) return;
        }
    } /* end of while */
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c|d precise ephemeris file
*          nav_t  *nav        IO  navigation data structure type
* return : status 0:err 0<:ok
*-----------------------------------------------------------------------------*/
extern int readsp3(const char *file, nav_t *nav)
{
    FILE *fp;
    gtime_t time={0};
    double bfact[2]={0},eint;
    int ns;
    char ftype=' ',htype[3],tsys[4]="";

	if(!(fp=fopen(file,"r"))) {
		return 0;
	}

	fgets(htype,3,fp);

	/* read sp3 header */
	switch (htype[1]) {
		case 'c': ns=readsp3ch(fp,&time,bfact,&eint,&ftype,tsys); break;
		case 'd': ns=readsp3dh(fp,&time,bfact,&eint,&ftype,tsys); break;
	}
       
	/* read sp3 body */
	readsp3b(fp,ns,bfact,ftype,tsys,nav);

	fclose(fp);

	return 1;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n)
{
    int i,j;
    
    for (j=1;j<n;j++) {
        for (i=0;i<n-j;i++) {
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav,
				   double *rs, double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3];
    int i,j,k,index;
    
    rs[0]=rs[1]=rs[2]=dts[0]=0.0;
    
    if (nav->n<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->n-1].time)>MAXDTE) {
        writelog("no precise ephemeris %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search (find previous given epoch)*/
    for (i=0,j=nav->n-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;
    
    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX>=nav->n) i=nav->n-NMAX-1; /*just first OR last*/
    
    for (j=0;j<=NMAX;j++) {
        t[j]=timediff(nav->peph[i+j].time,time);
        if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
            writelog("precise ephemeris outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
    }
    for (j=0;j<=NMAX;j++) {
        pos=nav->peph[i+j].pos[sat-1];

        p[0][j]=pos[0];
        p[1][j]=pos[1];
        p[2][j]=pos[2];
    }
    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);
    }

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];
    
	if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[index+i].std[sat-1][3];
    }
    else {
        dts[0]=0.0;
		return 0;
    }
    if (varc) *varc=SQR(std);

    return 1;
}
/* satellite clock by precise ephemeris-clock ----------------------------------
* compute satellite clock with precise ephemeris-clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dts         O   sat clock {bias} (s)
* return : status (1:ok,0:error or data outage)
*-----------------------------------------------------------------------------*/
static int satclk(gtime_t time, int sat, const nav_t *nav, double *dt)
{
	gtime_t time_tt;
	double rs[6],rss[3],rst[3],dtss[1],dtst[1],tt=1E-3;
	int i;

	if (sat<=0||MAXSAT<sat) return 0;

	/* satellite position and clock bias at t1*/
	if (!pephpos(time,   sat,nav,rss,dtss,NULL,NULL)) return 0;

	/* satellite position and clock bias at t2*/
	time_tt=timeadd(time,tt);
	if (!pephpos(time_tt,sat,nav,rst,dtst,NULL,NULL)) return 0;

	for (i=0;i<3;i++) {
		rs[i]=rss[i];
		rs[i+3]=(rst[i]-rss[i])/tt;
	}
	/* relativistic effect correction */
	if (dtss[0]!=0.0) {
		dt[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
	}

	return 1;
}
/* satellite position-clock by precise ephemeris-clock -------------------------
* compute satellite position-clock with precise ephemeris-clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts         O   sat clock {bias} (s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph and nav->n must be initialized
*		   by callling readsp3()
*-----------------------------------------------------------------------------*/
static int peph2pos(gtime_t time, int sat, const nav_t *nav,
                    double *rs, double *dts, double *var)
{
    gtime_t time_tt;
    double rss[3],rst[3],dtss[1],dtst[1],vare=0.0,varc=0.0,tt=1E-3;
    int i;
    
    if (sat<=0||MAXSAT<sat) return 0;
    
    /* satellite position and clock bias */
    if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)) return 0;
    
    time_tt=timeadd(time,tt);
    if (!pephpos(time_tt,sat,nav,rst,dtst,NULL,NULL)) return 0;
    
    for (i=0;i<3;i++) {
        rs[i  ]=rss[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }
    /* relativistic effect correction */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
        dts[1]=(dtst[0]-dtss[0])/tt;
    }
    else { /* no precise clock */
        dts[0]=dts[1]=0.0;
    }
    if (var) *var=vare+varc;
    
    return 1;
}
/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          double *rs       O   satellite positions and velocities (ecef)
*          double *dts      O   satellite clocks
*          double *var      O   sat position and clock error variances (m^2)
* return : none
* notes  : rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]        = obs[i] sat position and clock error variance (m^2)
*          if no navigation data, set 0 to rs[], dts[], var[]
*          satellite position is referenced to CoM of satellite
*-----------------------------------------------------------------------------*/
extern void satposs(const obsd_t *obs, int n, const nav_t *nav,
                    double *rs, double *dts, double *var)
{
	gtime_t time[MAXOBS]={{0}};
    double pr,dt;
    int i,j;
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<6;j++) rs [j+i*6]=0.0;
        for (j=0;j<2;j++) dts[j+i*2]=0.0;
        var[i]=0.0;
        
        /* search any psuedorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;
        
        if (j>=NFREQ) {
            writelog("no pseudorange %s sat=%2d\n",time_str(obs[i].time,3),obs[i].sat);
            continue;
        }
        /* transmission time by satellite clock */
        time[i]=timeadd(obs[i].time,-pr/CLIGHT);

		/* satellite clock bias by precise ephemeris */
		if (!satclk(time[i],obs[i].sat,nav,&dt)) {
			writelog("no clock bias %s sat=%2d\n",time_str(obs[i].time,3),obs[i].sat);
			continue;
		}
		time[i]=timeadd(time[i],-dt);

        /* satellite position and clock at transmission time */
		if (!peph2pos(time[i],obs[i].sat,nav,rs+i*6,dts+i*2,var+i)) {
            writelog("no ephemeris %s sat=%2d\n",time_str(obs[i].time,3),obs[i].sat);
            continue;
        }
     }
}