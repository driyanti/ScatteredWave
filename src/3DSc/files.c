#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#define HEADERLENGTH 240
#define MAXSAMP 1024 

#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

FILE *fid;
static FILE *fid2;
char **environ;
extern int errno;

void makedir_(char *dir,int dirlen)
{ 
    char tempdir[256];
	strncpy(tempdir,dir,dirlen); tempdir[dirlen] = '\0';
	mkdir(tempdir,0777);
	return;
}

void cleandir_(char *dir,int dirlen)
{
	char tempdir[256],comm[256];
	strncpy(tempdir,dir,dirlen); tempdir[dirlen] = '\0';
	sprintf(comm,"rm -r %s",tempdir);
	system(comm);
	return;
}

void readheader_(int *isrc,int *nsamp,float *deltat, char *basnam, \
                 int basnamlen)
{
    char header[HEADERLENGTH] = {0},sufile[20],tempbasnam[256];
    unsigned short ns,dt;
    FILE *fid3;

	strncpy(tempbasnam,basnam,basnamlen);
	tempbasnam[basnamlen] = '\0';
    sprintf(sufile,"%s_%02d.su",tempbasnam,*isrc);
    fid3 = fopen(sufile,"rb"); 
    if (fid3 == NULL) 
    {
		fprintf(stderr,"File %s not present\n",sufile);
		exit(-1);
    }
    fread((void *)header,(size_t)1,(size_t)114,fid3);
    fread((void *)&ns,(size_t)2,(size_t)1,fid3);
    fread((void *)&dt,(size_t)2,(size_t)1,fid3);
    fread((void *)header,(size_t)1,(size_t)122,fid3);
    *nsamp = (int)ns;
    *deltat = (float)dt * 1e-6;
	fclose(fid3);
}
FILE *fid;
static FILE *fid2;
char **environ;
extern int errno;

void writesutrace_(int *isrc, int *nsamp, float *timetrace, 
  float *sx, float *sy, float *sdepth, float *gx, float *gy,
  float *dt, int *srclen, char *basnam)
{
/* Makes a header, and puts the trace behind the header */
	static int curr=0,tracl=0;
	int isx,isy,igx,igy,isz,fldr,len;
	unsigned short atsamp,dtms; 
	short scalco=-1000;
	char header[HEADERLENGTH] = {0},suoutfile[20],tempbasnam[256];
 	static char prevbasnam[256] = "";
	unsigned short trid=1,counit=1;
	float scale=1000.0;
	
	len=40;	
	strncpy(tempbasnam,basnam,len);
	tempbasnam[len] = '\0';
	if (*srclen == 4) 
		sprintf(suoutfile,"%s_%04d.su",tempbasnam,*isrc);
	else
		sprintf(suoutfile,"%s_%02d.su",tempbasnam,*isrc);

	if (curr!=*isrc || strcmp(prevbasnam,tempbasnam))
	{
		curr = *isrc;
		strcpy(prevbasnam,tempbasnam);
		remove(suoutfile); 
		tracl=0;
	}
	fid = fopen(suoutfile,"a+");
	if (fid==NULL) printf("writing not possible, errno=%d\n",errno);

	tracl++;
	fldr = *isrc;
	atsamp = (unsigned short)(*nsamp); 
	dtms = (unsigned short)(*dt * 1e+6);
	isx = (int)(*sx * scale);
	isy = (int)(*sy * scale);
	igx = (int)(*gx * scale);
	igy = (int)(*gy * scale);
	isz = (int)(*sdepth * scale);

	fwrite(&tracl   ,(size_t)4,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)4  ,fid);
	fwrite(&fldr    ,(size_t)4,(size_t)1  ,fid);
	fwrite(&tracl   ,(size_t)4,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)12 ,fid); 
    fwrite(&trid    ,(size_t)2,(size_t)1  ,fid);
    fwrite(header   ,(size_t)1,(size_t)18 ,fid);
    fwrite(&isz     ,(size_t)4,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)16 ,fid);
	fwrite(&scalco  ,(size_t)2,(size_t)1  ,fid);
	fwrite(&scalco  ,(size_t)2,(size_t)1  ,fid);
	fwrite(&isx     ,(size_t)4,(size_t)1  ,fid);
	fwrite(&isy     ,(size_t)4,(size_t)1  ,fid);
	fwrite(&igx     ,(size_t)4,(size_t)1  ,fid);
	fwrite(&igy     ,(size_t)4,(size_t)1  ,fid);
    fwrite(&counit  ,(size_t)2,(size_t)1  ,fid); 
    fwrite(header   ,(size_t)1,(size_t)24 ,fid); 
	fwrite(&atsamp  ,(size_t)2,(size_t)1  ,fid);
	fwrite(&dtms    ,(size_t)2,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)122,fid); 
	fwrite(timetrace,(size_t)4,(size_t)(*nsamp),fid);
	fclose(fid);
}

void writesugrid_(int *isrc, int *nsamp, float *timetrace, 
  float *sx, float *sy, float *sdepth, float *gx, float *gy,
  float *dx, int *srclen, char *basnam)
{
/* Makes a header, and puts the trace behind the header */
	static int curr=0,tracl=0;
	int isx,isy,igx,igy,isz,fldr,len;
	unsigned short atsamp,dxs; 
	short scalco=-1000;
	char header[HEADERLENGTH] = {0},suoutfile[60],tempbasnam[256];
 	static char prevbasnam[256] = "";
	unsigned short trid=1,counit=1;
	float scale=1000.0;
	
	len=40;	
	strncpy(tempbasnam,basnam,len);
	tempbasnam[len] = '\0';
	if (*srclen == 4) 
		sprintf(suoutfile,"%s_%04d.su",tempbasnam,*isrc);
	else
		sprintf(suoutfile,"%s_%02d.su",tempbasnam,*isrc);

	if (curr!=*isrc || strcmp(prevbasnam,tempbasnam))
	{
		curr = *isrc;
		strcpy(prevbasnam,tempbasnam);
		remove(suoutfile); 
		tracl=0;
	}
	fid = fopen(suoutfile,"a+");
	if (fid==NULL) printf("writing not possible, errno=%d\n",errno);

	tracl++;
	fldr = *isrc;
	atsamp = (unsigned short)(*nsamp); 
	dxs  = (unsigned short)(*dx * scale);
	isx = (int)(*sx * scale);
	isy = (int)(*sy * scale);
	igx = (int)(*gx * scale);
	igy = (int)(*gy * scale);
	isz = (int)(*sdepth * scale);

	fwrite(&tracl   ,(size_t)4,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)4  ,fid);
	fwrite(&fldr    ,(size_t)4,(size_t)1  ,fid);
	fwrite(&tracl   ,(size_t)4,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)12 ,fid); 
    fwrite(&trid    ,(size_t)2,(size_t)1  ,fid);
    fwrite(header   ,(size_t)1,(size_t)18 ,fid);
    fwrite(&isz     ,(size_t)4,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)16 ,fid);
	fwrite(&scalco  ,(size_t)2,(size_t)1  ,fid);
	fwrite(&scalco  ,(size_t)2,(size_t)1  ,fid);
	fwrite(&isx     ,(size_t)4,(size_t)1  ,fid);
	fwrite(&isy     ,(size_t)4,(size_t)1  ,fid);
	fwrite(&igx     ,(size_t)4,(size_t)1  ,fid);
	fwrite(&igy     ,(size_t)4,(size_t)1  ,fid);
    fwrite(&counit  ,(size_t)2,(size_t)1  ,fid); 
    fwrite(header   ,(size_t)1,(size_t)24 ,fid); 
	fwrite(&atsamp  ,(size_t)2,(size_t)1  ,fid);
	fwrite(&dxs     ,(size_t)2,(size_t)1  ,fid);
	fwrite(header   ,(size_t)1,(size_t)122,fid); 
	fwrite(timetrace,(size_t)4,(size_t)(*nsamp),fid);
	fclose(fid);
}

