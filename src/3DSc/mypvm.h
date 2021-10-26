      integer   MAXPARTS
      parameter(MAXPARTS=10)

      integer   nproc,proctid,KILLMSG,MYMSGTAG,
     .          MYCOMMTAG,MYSRCMSGTAG,MYINTERMSGTAG,
     .          MYRESMSGTAG,MYFINMSGTAG,parentid
      dimension proctid(MAXPARTS)

      parameter (KILLMSG=666,MYMSGTAG=11,
     .          MYRESMSGTAG=12,MYCOMMTAG=13,MYSRCMSGTAG=14,
     .          MYFINMSGTAG=15,MYINTERMSGTAG=16)

      integer   COMM_STOP,COMM_CALC,COMM_FINNODE
      parameter (COMM_FINNODE=60,COMM_CALC=61,COMM_STOP=62)
      
      
      integer   BOOLEAN
      parameter (BOOLEAN=INTEGER4)
      
      common    /myp1/ nproc,proctid,parentid
