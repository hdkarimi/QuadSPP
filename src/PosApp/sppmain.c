/* sppmain.c : main program compute single point position
*
* version : yyyy/mm/dd		ver		description
*			2018/06/12		1.0		create
*
* Copyright (C) 2018 by H.KARIMI, All rights reserved.
*-----------------------------------------------------------------------------*/
#include "gnsslib.h"

/* help text -----------------------------------------------------------------*/
static const char *help[]={
"",
" usage: spp [file.rnx] [file.sp3] [option...]",
"",
" Read RINEX OBS and SP3 files and compute receiver position for each epoch. ",
" For observation file, RINEX Version 3 is supported and for navigation file ",
" precise orbit and clock file, SP3c and SP3d are supported in this program. ",
" The extension of the RINEX OBS shall be .yyo or .rnx and the extension of  ",
" SP3 shall be .sp3 or .eph. Command line options are as follows ([]:default)",
"",
" -hh          print help",
" -ts  ds ts   start date/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te  de te   end date/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti  tint    time interval (sec) [all]",
" -sys s[s...] nav system(s) (s=G:GPS,R:GLO,E:GAL,C:BDS) [G]",
" -sf  format  output position format (0=lat/lon/hgt or 1=x/y/z-ecef) [1]",
" -el  mask    elevation mask angle (deg) [15]",
" -gd  gdop    max GDOP for reject epoch solution [30]",
"",
"example: spp mmmmddd0.rnx aaawwwwd.sp3 -ts yyyy/mm/dd hh:mm:ss -el 20",
"example: spp mmmmddd0.rnx aaawwwwd.sp3 -sys C -sf 1 -el 10 -gd 30",
};
/* show message --------------------------------------------------------------*/
extern int showmsg(char *format, ...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}
/* spp main ------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    opt_t prcopt={	/* processing option default */
		SYS_GPS ,			/* navigation system */
		SOLF_XYZ,			/* solution format (SOLF_xxx) */
		15*D2R  ,			/* elevation mask angle (rad) */
		30      ,			/* reject threshold of gdop */
		"-sys=G",			/* rinex option */
		"\\t"				/* field separator */
	};		

    gtime_t ts={0},te={0};
    double tint=0.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59};
    int i,n,ret;
	char *infile[2],*p;

    for (i=1,n=0;i<argc;i++) {
        if      (!strcmp(argv[i],"-ts")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",es,es+1,es+2);
            sscanf(argv[++i],"%lf:%lf:%lf",es+3,es+4,es+5);
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",ee,ee+1,ee+2);
            sscanf(argv[++i],"%lf:%lf:%lf",ee+3,ee+4,ee+5);
            te=epoch2time(ee);
        }
        else if (!strcmp(argv[i],"-ti")&&i+1<argc) tint=atof(argv[++i]);
        else if (!strcmp(argv[i],"-sys")&&i+1<argc) {
			prcopt.navsys=SYS_NON;
			strcpy(prcopt.rnxopt,"-sys=");
			strcat(prcopt.rnxopt,argv[++i]);
            for (p=argv[i];*p;p++) {
                switch (*p) {
                    case 'G': prcopt.navsys|=SYS_GPS; break;
                    case 'R': prcopt.navsys|=SYS_GLO; break;
                    case 'E': prcopt.navsys|=SYS_GAL; break;
                    case 'C': prcopt.navsys|=SYS_BDS; break;
                }
            }
        }
        else if (!strcmp(argv[i],"-sf")&&i+1<argc) prcopt.posf   =atoi(argv[++i]);
	    else if (!strcmp(argv[i],"-el")&&i+1<argc) prcopt.elmin  =atof(argv[++i])*D2R;
		else if (!strcmp(argv[i],"-gd")&&i+1<argc) prcopt.maxgdop=atof(argv[++i]);
        else if (!strcmp(argv[i],"-hh")) printhelp();
        else if (n<2) infile[n++]=argv[i];
    }
    if (n<=0) {
        showmsg("error : no input file");
        return -2;
    }

    ret=postpos(ts,te,tint,&prcopt,infile);
    if (!ret) {fprintf(stderr,"%40s\r",""); showmsg("done!");}
    return ret;
}