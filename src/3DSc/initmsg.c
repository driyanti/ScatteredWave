#include <stdio.h>
#include <stdlib.h>
#define GETMSG "/home/mf/ditzel/Tools/getmsg"

/* Sets standard output of this process and its children to its
 * argument 'filnam'. Should be called from FORTRAN
 * 
 * Fabian Ernst, 220496
 */

char hostnam[20],newfilnam[256];
char *argarray[] = {"log_file",(char *)0};

void initmsg_(char *filnam,int len)
{
	int mytid,numt,logtid;
	strncpy(newfilnam,filnam,256);
	newfilnam[len] = '\0';
	printf("getmsg: starting %s\n",GETMSG); 
    printf("getmsg: redirecting node output to file %s\n",newfilnam);
	while(len>0 && isspace(newfilnam[len-1])) newfilnam[len-1] = '\0';
	argarray[0] = (char *)malloc((len+4)*sizeof(char));
        strcpy(argarray[0],newfilnam);
	strcpy(hostnam,getenv("HOST"));
	mytid = pvm_mytid();
  	numt = pvm_spawn(GETMSG,argarray,1,hostnam,1,&logtid); 
	if (numt<1) 
	{
         fprintf(stderr,"getmsg: could not spawn output-process\n");
		 fprintf(stderr,"tid = %d\n",logtid);
         return;
	}
	pvm_setopt(10,logtid);  /* set output of this task */
	pvm_setopt(11,55);      /* output msgid ( != 66)   */
	pvm_setopt(4,logtid);   /* set output of children  */
	return;
}

void  getname_(char *machinename)
{
     strcpy(machinename,getenv("HOST"));
}

