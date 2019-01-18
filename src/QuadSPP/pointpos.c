/*------------------------------------------------------------------------------
* pointpos.c : point positioning functions
*
* references : 
*	 (1) J.Kouba, A Guide to using International GNSS Service prouducts, 2009
*    (2) IGS Multi-GNSS Experiments (MGEX) (http://igs.org/mgex)
*
* version : yyyy/mm/dd		ver		description
*			2018/04/20		1.0		create
*
* Copyright (C) 2018 by H.KARIMI, All rights reserved.
*-----------------------------------------------------------------------------*/
#include "gnsslib.h"

/* constants/variable --------------------------------------------------------*/
#define SQR(x)      ((x)*(x))
#define NX          4	        /* # of estimated parameters */
#define MAXITR      10          /* max number of iteration for point pos */
#define REL_HUMI    0.5         /* relative humidity for saastamoinen model */

const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
    10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
    31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
    46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
    61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
    74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
    88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
    101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
    113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
    126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
    138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};

/* pseudorange measurement error variance ------------------------------------*/
static double comerr(double el, int sys)
{
    double fact,varr;
	//fact=sys==SYS_GLO?1.5:1.0;
    //varr=SQR(100.0)*(SQR(0.003)+SQR(0.003)/sin(el));
    //varr*=SQR(3.0); /* iono-free */

	varr=1/sin(el);
	return varr;
}
/* psendorange with iono-free LC ---------------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, double *var)
{
    const double *frq=nav->frq[obs->sat-1];
    double PC,P1,P2,fi,fj;
    int i=0,j=1,sys;
    
    *var=0.0; /* iono-free LC variance */
    
    if (!(sys=satsys(obs->sat,NULL))) return 0.0;
    
    /* L1-L2 for GPS/GLO/BDS, L1-L5 for GAL */
    if (NFREQ>=3&&(sys&(SYS_GAL))) j=2;
    
    if (NFREQ<2||frq[i]==0.0||frq[j]==0.0) return 0.0;
    
    fi=SQR(frq[i]);	/* fi^2 */
	fj=SQR(frq[j]);	/* fj^2 */
    P1=obs->P[i];
    P2=obs->P[j];

	if (P1==0.0||P2==0.0) {
		writelog("iono-free missed at %s sat=%3d\n",
			     time_str(obs->time,0),obs->sat); 
		return 0.0;
	}

	/* iono-free combination */
	PC = (fi*P1-fj*P2)/(fi-fj);
    
	*var=SQR(0.3);

    return PC;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *trp      O   tropospheric delay (m)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
static double tropcorr(gtime_t time, const double *pos, const double *azel,
					   opt_t *opt, double *dtrp, double *vart)
{
    double zhd,zwd,m_h,m_w;
    
    /* zenith hydrostatic delay */
    zhd=tropmodel(time,pos,REL_HUMI,&zwd);
    
    /* mapping function */
    m_h=tropmapf (time,pos,azel,&m_w);

	*dtrp=m_h*zhd + m_w*zwd;			/* tropspheric delay */
	*vart=SQR(0.3/(sin(azel[1])+0.1));	/* saastamoinen model variance */

    return 1.0;
}
/* pseudorange residuals -----------------------------------------------------*/
static int rescode(int iter, const obsd_t *obs, int n, 
			       const double *rs, const double *dts, const double *vare,
                   const nav_t *nav, const double *x, const opt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
				   double *resp, int *ns)
{
    double r,dtrp,vmeas,vtrp,rr[3],pos[3],dtr,e[3],P;
    int i,j,nv=0,sys,mask[4]={0};

    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        
        if (!(sys=satsys(obs[i].sat,NULL))) continue;

        /* geometric distance/azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr,e))<=0.0||
            satazel(pos,e,azel+i*2)<opt->elmin) continue;
        
        /* psudorange with iono-free LC */
        if ((P=prange(obs+i,nav,&vmeas))==0.0) continue;

        /* tropospheric corrections */
        if (!tropcorr(obs[i].time,pos,azel+i*2,opt,&dtrp,&vtrp)) {
            continue;
        }

		/* x, y, z or dts of sat is zero */
        if (dts[i*2]==0.0) continue;

        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dtrp);	/* dl=l-l0 */
        
        /* design matrix */						/* H(x-x0) */
        for (j=0;j<NX;j++) H[j+nv*NX]=j<3?-e[j]:1.0;
        
		/* time system and receiver bias offset correction */
        //if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        //else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        //else if (sys==SYS_BDS) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        //else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* error variance */
		//var[nv++]=comerr(azel[1+i*2],sys)+vare[i]+vmeas+vtrp;
	    var[nv++]=comerr(azel[1+i*2],sys);
    }
	/* constraint to avoid rank-deficient */
    //for (i=0;i<4;i++) {
    //    if (mask[i]) continue;
    //    v[nv]=0.0;
    //    for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
    //   var[nv++]=0.01;
    //}
    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const opt_t *opt, const double *v, int nv, int nx)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;

	/* chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        writelog("chi-square error nv=%d vv=%.1f cs=%.1f\n",nv,vv,chisqr[nv-nx-1]);
        return 0;
    }
    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        writelog("gdop error nv=%d gdop=%.1f\n",nv,dop[0]);
        return 0;
    }
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const nav_t *nav, const opt_t *opt, 
				  sol_t *sol, double *azel, int *vsat, double *resp)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;

    v=mat(n,1); H=mat(n,NX); var=mat(n,1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    
    for (i=0;i<MAXITR;i++) {
        
        /* pseudorange residuals */
        nv=rescode(i,obs,n,rs,dts,vare,nav,x,opt,v,H,var,azel,vsat,resp,&ns);
        
        if (nv<NX) {
            writelog("lack of valid sats ns=%d\n",nv);
            break;
        }
        /* weight by variance */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            writelog("lsq error info=%d\n",info);
            break;
        }
        for (j=0;j<NX;j++) x[j]+=dx[j];
        
        if (norm(dx,NX)<1E-4) {
            sol->time=obs[0].time;
			sol->dtr[0]=x[3]/CLIGHT;					    /* dtr (s) */
            for (j=0;j<3;j++) sol->rr[j]=x[j];				/* x,y,z (m) */
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];	/* var x,y,z */
            sol->qr[3]=(float)Q[1];							/* cov xy */
            sol->qr[4]=(float)Q[2+NX];						/* cov yz */
            sol->qr[5]=(float)Q[2];							/* cov zx */
			sol->ns=(unsigned char)n;						/* all satellite */
            sol->nv=(unsigned char)ns;						/* used satellite */
            
			/* validate solution */
            stat=valsol(azel,vsat,n,opt,v,nv,NX);

            free(v); free(H); free(var);
            
            return stat;
        }
    }
    if (i>=MAXITR) writelog("iteration divergent i=%d\n",i);
    
    free(v); free(H); free(var);
    
    return 0;
}
/* single-point positioning ----------------------------------------------------
* compute receiver position and clock bias by single-point positioning
* with pseudorange observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs,int n,const nav_t *nav,const opt_t *opt,sol_t *sol)
{
    double *rs,*dts,*var,*azel,*resp;
    int i,stat,vsat[MAXOBS]={0};
    
    if (n<=0) {writelog("no observation data\n"); return 0;}
    
    rs=mat(6,n);dts=mat(2,n);var=mat(1,n);azel=zeros(2,n);resp=mat(1,n);

    /* satellite positons, velocities and clocks */
    satposs(obs,n,nav,rs,dts,var);
    
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,rs,dts,var,nav,opt,sol,azel,vsat,resp);
    
    free(rs); free(dts); free(var); free(azel); free(resp);
    return stat;
}