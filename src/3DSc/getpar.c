/* getpar(file,name,val): routines to get a parameter from a file 
 *
 * contains:
 *  getparint(const char *name,int *val)
 *  getparfloat(const char *name,float *val)
 *  getparbool(const char *name,int *val)
 *  getparstring(const char *name,char *val)
 *  getpararray(const char *name,int *vals,int *lenght_of_vals)
 * 
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#define REGELLEN 256

static FILE *fp;
char line[REGELLEN],chix[REGELLEN];

int openfile(char *hfile,int len)
{ 
	char file[REGELLEN];
	strncpy(file,hfile,REGELLEN);
	file[len] = '\0';
	fp = fopen(file,"r");
    if (fp==NULL)
		fprintf(stderr,"File %s could not be getparred\n",file);
	return fp==NULL ? 0 : 1;
}

int seekline(char *name2,int len)
{
	char s[REGELLEN],s2[REGELLEN],name[REGELLEN];
	strncpy(name,name2,REGELLEN);
	name[len] = '\0';
	if (fp==NULL) return 0;
	while(fgets(s,256,fp)!=NULL)
	{
		sscanf(s," %s",s2);
		if(!strncmp(s2,name,len)) 
		{
			fclose(fp);
			strcpy(line,s);
			return 1;
		}
	}
	fclose(fp);
	return 0;
}


void getparint_(char *file,char *name,int *val,int len1,int len2)
{
	int hval;
	if (!openfile(file,len1)) return;
	if (seekline(name,len2)) 
		strcpy(chix,strchr(line,'='));
	else
		return;
	if (strlen(chix)<1) return;
	if (sscanf(chix+1,"%d",&hval)>0)
		*val = hval;
	printf("variable %1$.*2$s, value %3$d\n",name,len2,*val);
	return;
}
	
void getparbool_(char *file,char *name,int *val,int len1,int len2)
{
	char hval;
	if (!openfile(file,len1)) return;
	if (seekline(name,len2)) 
		strcpy(chix,strchr(line,'='));
	else
	{
/*		printf("Variable %s not found in file %s\n",name,file); */
		return;
	}
	if (strlen(chix)<1) return;
	sscanf(chix+1," %c",&hval);
	*val = (tolower(hval)=='y' || tolower(hval)=='j' || hval=='1' || tolower(hval)=='t');
	printf("variable %1$.*2$s, value %3$s\n",name,len2,(*val ? "true" : "false") );
	return;
}
		
void getparfloat_(char *file,char *name,float *val,int len1,int len2)
{
	float hval;
	if (!openfile(file,len1)) return;
	if (seekline(name,len2)) 
		strcpy(chix,strchr(line,'='));
	else
		return;
	if (strlen(chix)<1) return;
	if (sscanf(chix+1,"%f",&hval)>0)
		*val = hval;
	printf("variable %1$.*2$s, value %3$f\n",name,len2,*val);
	return;
}

void getparstring_(char *file,char *name,char *val,int *slen,int len1,int len2,int len3)
{
	if (!openfile(file,len1)) return;
	if (seekline(name,len2)) 
		strcpy(chix,strchr(line,'='));
	else	
		return;
	if (strlen(chix)<1) return;
	strncpy(val,chix+1,len3);
	*slen = (int)strlen(val);
	return;
}


void getpararray_(char *file,char *name,int *vals,int *arrlen,int len1,int len2) 
{
	char u[REGELLEN],*u1,*u2,subs[REGELLEN];
	int currnum,prevnum,ix,i;
	for(i=0;i<*arrlen;vals[i++]=0); /* initialize to zero */
	if (!openfile(file,len1)) return; /* file not found */
	if (!seekline(name,len2)) return; /* variable not found */
	strncpy(subs,strchr(line,'='),REGELLEN);
	if (!strcmp(subs,"")) return;
	strcpy(line,subs+1);
	printf("variable %1$.*2$s, value %3$s\n",name,len2,line);
/* Parsing of line. Allowed is: 1,2,3; 1-3,4 */
	strcpy(u,line);
   	while (1)
	{
		u2 = strchr(u,','); /* search for ',' */
		if (u2==NULL)  /* if not found: end of line */
			u2 = strchr(u,'\0');
	  	ix = (int)(u2-u); 
		/* Parse substring */
		strcpy(subs,u);
		subs[ix] = '\0'; 
		/* First: check out if there is a number */
		if (sscanf(subs,"%d",&prevnum)>0)
		{
			u1 = strchr(subs,'-');
			if (!(u1!=NULL && sscanf(u1+1,"%d",&currnum)>0))
				currnum = prevnum;
		 	for (i=prevnum;i<=currnum;i++)
			  if (i>0 && i<=*arrlen)
		 		vals[i-1] = 1; /* FORTRAN starts with 1, C with 0 :-( */
		}
		if (strlen(u2)>0) 	
			strcpy(u,u2+1);
		else
			return;
	}
}
