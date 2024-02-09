cdoc********************************************************************
c* name:  main_demult
c*
c* descr:  
c*
c* NOTES:
c*
c* HIST:   Version Date       Author   			Comment
c*         1.0     28 Jul 08  Christina Dwi Riyanti     Initial version
c*
cdoc********************************************************************

c***********************************************************************
c     declare variables
c***********************************************************************
      program main_dmult
        
      integer     kk,ns,ntrac
      integer     nf
      integer     nf1,nf2
      integer     nkmax
      integer     nsi,ntrc,nrec,noff
      integer     ii,jj,ij,i,j,k,l
      integer     nlay,nsamp,nwav,iom
      integer     nx,nmax
      integer     maxtrc,it,maxit
      real        df,fnum
      real        tol
      real        dw,pi
      real        f1,f2,dt0   
      real        dx,dt,x3r,x3s,eps,x1s
      real        xoff  
      real        rmax  
      real*8      TIME_1
      real*8      TIME_2
      real*8      rt,rcond
      parameter   (maxtrc = 512)  
      parameter   (nmax=50000)
      
      real        ttap(2048)
      real        xtap(2048)
      real        wave(2048),etap(2048)
      real        velo(100),rho(100),depth(100),t(2048)
      real        pupdat(maxtrc,1024)
      real        prdat(maxtrc,1024)
      real        vdodat(maxtrc,1024)
      real        temp1(nmax,1024),temp2(nmax,1024)
      complex     ci,ces,imat
      complex     cpup(2048),cvdo(2048)
      complex     cwave(2048)
      complex     cwave2(2048),ct(2048)
      complex     cpupdat(nmax,2048)
      complex     cvdodat(nmax,2048)

      integer     maxtr,ntot
      parameter   (maxtr=189)      
      integer     lda,ipvt(maxtr)
      complex*16  som,som1,som2
      complex*16  cvzdat(maxtr,maxtr)
      complex*16  matdat(maxtr,maxtr)
      complex*16  cpudat(maxtr,maxtr)
      complex*16  tmp(maxtr,maxtr)
      complex*16  bvect(maxtr),zz(maxtr),det(2)
      complex     cpmdat(maxtr,1024) 
       
      DATA lda/maxtr/

c***********************************************************************
c     read the geometry from the ascii file                     
c***********************************************************************

      open (10,file='geo4',status='unknown')
      read(10,*)nlay 
      do i=1,nlay
         read(10,*)depth(i),velo(i), rho(i)
      enddo
      close(10)

      CALL CPU_TIME(TIME_1)
 
c**********************************************************************C
c     define constants       
c**********************************************************************C

         pi = 4. * atan (1.)
         dt = 0.004
         nf = 2048 
      nsamp = 1024 
         ns = 1024 
         nx = 1024
        x3s = 20.
        x3r = 80.
        x1s = 50.
         dx = 20. 
         df = 1./(nf*dt)
         f1 = 2.
         f2 = 80.
        ci  = cmplx(0.,1.)
        c0  = cmplx(0.,0.)
        eps = 0.025 
        stab = 0.05
        ntrac= 189   
        ntot = ntrac*ntrac
        dt0 = 0.0
        nf1 = int (f1/df)
        nf2 = int (f2/df)
        dw  = 2.*pi*df

c**********************************************************************C
c     calculate the time-taper    
c**********************************************************************C

      call timetap(etap,nsamp,dt,eps)

c**********************************************************************C
c     calculate the spacial tapers    
c**********************************************************************C

      call gen_taper(ttap,nf,1,200,nsamp-400,nsamp,2.0)
      call gen_taper(xtap,nx,1,30,ntrac-30,ntrac,2.0)
c**********************************************************************C
c     calculate the wavelet spectrum
c**********************************************************************C

      open(15,file='wave.save',status='unknown')

      read(15,*)nwav

      do i=1,nwav
         read(15,*)wave(i) 
      enddo

      call vclr(cwave,1,nf)        
      
      do i=1,nwav  
        cwave(i)=cmplx (wave(i), 0.)
      enddo  
      
      call pfft(cwave,nf,-dt)

c... use cwave2 below to get an artifact-free vdown field ...
c... use cwave if you want the best pdemult result ...

      do iom=1,nf/2
        cwave2(iom) = cexp(cmplx(eps,(iom-1)*dw*dt0))*cwave(iom) 
      enddo
      
      OPEN(unit=70,file='Psup.b',access='direct',
     .     form='unformatted',recl=4*ns*ntot)

      read(70,rec=1)((temp1(i,kk),kk=1,ns),i=1,ntot)
      
      close(70)

      OPEN(unit=80,file='Vdo.b',access='direct',
     &     form='unformatted',recl=4*ns*ntot)

      read(80,rec=1)((temp2(i,kk),kk=1,ns),i=1,ntot)

      close(80)

C**********************************************************************C
C     go to frequency   
C**********************************************************************C

c     ... clear arrays ...

      call vclr(cpupdat,1,2*nx)
      call vclr(cvdodat,1,2*nx)

      do jj=1,ntot

         call vclr(cpup,1,2*nx)
         call vclr(cvdo,1,2*nx)

c        ... make complex ...
         do ii=1,ns
            cpup(ii) = cmplx(temp1(jj,ii),0.)
            cvdo(ii) = cmplx(temp2(jj,ii),0.)
         enddo
           	    
c        ... transform to frequency domain ...
         call pfft(cpup,nf,-dt)	 
         call pfft(cvdo,nf,-dt)

c        ... store result in data array ...
         do iom = 1, nf/2+1
            cpupdat(jj,iom) = cpup(iom)
            cvdodat(jj,iom) = cvdo(iom)
         enddo
      enddo

      write(*,*)'ready with the frequency-transformation ...'

C**********************************************************************C
C     start of frequency loop 
C**********************************************************************C
        
       it  = 0   
       do iom = 1,nf2
         WRITE(*,*) '- FREQ ',iom,' OF ',nf2

c......... calculate constants .......

        ces = cmplx(eps,(iom-1)*dw)
c         ces = cmplx(eps,iom*dw)
 
C*** Define vector B for a single shot( shot nr. 45 for example)

        do kk = 1,ntrac 
           do j = 1,ntrac
	     jj = (kk-1)*ntrac+j
 	     cpudat(j,kk) = cpupdat(jj,iom)*cwave2(iom)
	   end do
	     bvect(kk) = cpudat(45,kk)
	 end do  

C*** Define matrix A=Vdo(ii,jj)  

         do kk = 1,ntrac 
	   do j = 1,ntrac
             jj = (kk-1)*ntrac+j
 	     cvzdat(kk,j) = cvdodat(jj,iom)*2.*ces*rho(1)
	     matdat(j,kk) = cvzdat(j,kk)
	   end do
	 end do  

C*** calculate (A*conjg(A) + stab*I) 

         DO 120, J = 1, ntrac
            DO 110, I = 1, ntrac
                som = dcmplx(0.0,0.0)
	         if (I .eq. J) then
		  imat = cmplx(1.0,0.0)
		 else
		  imat = cmplx(0.0,0.0)
		 endif
              DO 100, L = 1, ntrac
                 som = som + cvzdat(I,L)*DCONJG(cvzdat(J,L))
  100         CONTINUE
                  cpudat( I, J ) = som + stab**2*imat
  110       CONTINUE
  120    CONTINUE
     
C***  calculate AA = (A*conjg(A) + stab*I)*(-1) 
            
      CALL zgeco(cpudat,lda,ntrac,ipvt,rcond,zz)

      CALL zgedi(cpudat,lda,ntrac,ipvt,det,zz,01)
      
c      WRITE(*,*)'rcond',rcond
	
c      rt = 1.0 + rcond	
	
c      IF (rt .EQ. 1.0) then 
c        WRITE(*,99)
c      ENDIF

C*** compute BB=conjg(A)*AA 
        
       DO  J = 1, ntrac
         DO I = 1, ntrac
            som1 = dcmplx(0.0,0.0)
           DO L = 1, ntrac
            som1 = som1 + DCONJG(cvzdat(L,I))*cpudat( L,J )
           END DO
            tmp( I, J ) = som1
         END DO
       END DO

C*** compute pdm= B*BB 

       DO I = 1, ntrac
            som2 = dcmplx(0.0,0.0)
         DO L = 1, ntrac
            som2 = som2 + bvect(L)*tmp(I,L)
         END DO
            cpmdat(I,iom) = som2*xtap(I)
        END DO

99    FORMAT(40H matrix is singular to working precision)

      END DO

c     ... end of frequency loop ...

      write(*,*)'ready with the model...'
     
c**********************************************************************C
c     go back to time and write files to disk 
c**********************************************************************C                
            
      do kk=1,ntrac

c        ... scale array to zero ...
         call vclr(ct,1,nf)        

c        ... set negative frequencies ...
         do iom = 1, nf2
            ct(iom) = cpmdat(kk,iom)
         enddo

c        ... transform back to time ...
         call pfft(ct,nf,df)

c        .... apply taper ...
         do i = 1,ns 
            prdat(kk,i) = 2. * real(ct(i))*ttap(i)
         enddo 

c        ... write out to file .........................

      enddo   

      write(*,*)'ready with the demultiple for a single source...'

      CALL CPU_TIME(TIME_2)
      
      WRITE(*,*) 'CPU_TIME=',TIME_2-TIME_1

      OPEN(unit=85,file='cprm_45i.b',access='direct',
     &     form='unformatted',recl=4*ns*ntrac)

      write(85,rec=1)((real(prdat(i,kk)),kk=1,ns),i=1,ntrac)
      
      close(85)

c***********************************************************************c
c     end of program    
c***********************************************************************c

      return 
      end  
c***********************************************************************
c     close all files

c***********************************************************************
c     some general utilities
c***********************************************************************

C====================================================================C
      subroutine timetap(t,n,dt,eps)

      implicit none

      integer i,n
      real    eps,dt,t(n)

c     .... calculate the time taper ....
      do i = 1,n
         t(i) = exp ((i-1)*dt*eps)
      enddo

C     ... end of subroutine ...
      return
      end
c***********************************************************************
      subroutine vclr(ct,l,n)
      
      integer l,n
      complex ct(n)
      
      do i=1,n
       ct(i) = cmplx(0.0,0.0)
      end do
      
      return
      end     
c***********************************************************************
      subroutine pfft(cx,n,delta)

c     ... implicit none ...

      integer j,i,N,m,istep,l
      complex cx(N),cw,ctemp
      real signi,sc,delta,arg

c     ... start of subroutine ...
      j=1
      if (delta.lt.0.) then
          signi     = -1.
      else
          signi     = 1.
      endif
      sc  = abs(delta)
      do 630 i=1,N
      if(i.gt.j) goto 610
      ctemp=cx(j)*sc
      cx(j)=cx(i)*sc
      cx(i)=ctemp
610   m=N/2
620   if(j.le.m) goto 630
      j=j-m
      m=m/2
      if(m.ge.1) goto 620
630   j=j+m
      l=1
640   istep=2*l
      do 650 m=1,l
c     carg=cmpN(0.,1.)*(3.141592653*signi*(m-1))/l
      arg=3.141592653*signi*(m-1)/l
      cw=cmplx(cos(arg),sin(arg))
      do 650 i=m,N,istep
      ctemp=cw*cx(i+l)
      cx(i+l)=cx(i)-ctemp
650   cx(i)=cx(i)+ctemp
      l=istep
      if(l.lt.N) goto 640

c     ... end of subroutine ...
      return
      end
      
c***********************************************************************
      subroutine gen_taper(taper,lentap,itap1,itap2,itap3,itap4,power)
      
cc*******create cos**power taper from itap1 to itap2 and itap3 to itap4
 
      implicit none

c***********************************************************************
c     subroutine variables
c***********************************************************************

      integer itap1,itap2,itap3,itap4,lentap,ix
      real power
      real taper(*)

c***********************************************************************
c     local variables
c***********************************************************************

      real tapfac,pi
      character*80 text

c***********************************************************************
c     define pi
c***********************************************************************

      pi=3.141592654

c***********************************************************************
c     check input variables
c***********************************************************************

      if (lentap.le.0) then
         write(*,*)'length taper < 1'
      endif
      if ( itap2.lt.itap1 .or. itap3.lt.itap2 .or. itap4.lt.itap3 ) then
         write(*,*)'gen_taper: taper values should increase'
         write(*,*)'wrong taper values!'
      endif

c***********************************************************************
c     create taper, start with zero part until itap1
c***********************************************************************

      do 100 ix=1,itap1-1
         taper(ix)=0.
  100 continue

c***********************************************************************
c     cosine edge from itap1 to itap2
c***********************************************************************

      tapfac=0.5*pi/float(max(1,itap2-itap1+1))
      do 200 ix=itap1,itap2
         taper(ix)=abs(cos(float(ix-itap2)*tapfac))**power
  200 continue

c***********************************************************************
c     flat part with value 1 from itap2 to itap3
c***********************************************************************

      do 300 ix=itap2+1,itap3-1
         taper(ix)=1.
  300 continue

c***********************************************************************
c     cosine edge from itap3 to itap4
c***********************************************************************

      tapfac=0.5*pi/float(max(1,itap4-itap3+1))
      do 400 ix=itap3,itap4
         taper(ix)=abs(cos(float(ix-itap3)*tapfac))**power
  400 continue

c***********************************************************************
c     zero beyond itap4
c***********************************************************************

      do 500 ix=itap4+1,lentap
         taper(ix)=0.
  500 continue

c***********************************************************************
c     end of subroutine GEN_TAPER
c***********************************************************************

      return
      end

C*********************************************************************C
      SUBROUTINE INVGAU(A,N,EPS)
C*********************************************************************C
C								                C
C    PURPOSE							                C
C    TO INVERT A SQUARE MATRIX A				                C
C								                C
C    USAGE							                C
C    CALL INVGAU(A,N,EPS,IIF)					                C
C								                C
C    DESCRIPTION OF PARAMETERS					                C
C    A	  (I) ARRAY OF DIMENSION (N,N)				                C
C	  (O) THE RESULTANT INVERSE MATRIX			                C
C    N	  (I) NUMBER OF EQUATIONS    =< 30			                C
C    EPS  (I) REAL TOLERANCE VALUE				               C
C								               C
C    METHOD							               C
C    GAUSS ELIMINATION WITH COLUMN PIVOTTING			                C
C								               C
C*********************************************************************C
      REAL	  EPS,PIV
      COMPLEX     A(N,N),B(N),C(N),T,R
      INTEGER	  IIF,I,J,K,N,P(N),Q(N)

      IFF=0
C** SEARCH FOR LARGEST ELEMENT
      DO 19 K=1,N
      PIV=0D0
      DO 4 I=K,N
      DO 4 J=K,N
      T=A(I,J)
      IF (PIV.GE.ABS(T)) GOTO 4
      P(K)=I
      Q(K)=J
      PIV=ABS(T)
   4  CONTINUE
      IF (PIV.GE.EPS) GOTO 6
      IIF=1
      GOTO 100
C** INTERCHANGE COLUMNS
   6  I=P(K)
      IF (K.EQ.I) GOTO 10
      DO 8 J=1,N
      T=A(I,J)
      A(I,J)=A(K,J)
   8  A(K,J)=T
C** INTERCHANGE ROWS
  10  J=Q(K)
      IF (K.EQ.J) GOTO 14
      DO 12 I=1,N
      T=A(I,J)
      A(I,J)=A(I,K)
  12  A(I,K)=T
C** DIVIDE BY PIVOT
  14  R=1D0/A(K,K)
      DO 16 J=1,N
      IF (J.EQ.K) THEN
      B(J)=R
      C(J)=CMPLX(1D0,0D0)
      ELSE
      B(J)=-A(K,J)*R
      C(J)=A(J,K)
      END IF
      A(K,J)=CMPLX(0D0,0D0)
  16  A(J,K)=CMPLX(0D0,0D0)
      DO 18 J=1,N
      DO 18 I=1,N
      A(I,J)=A(I,J)+C(I)*B(J)
  18  CONTINUE
  19  CONTINUE
C** DETERMINE INVERSE
      DO 36 K=N,1,-1
      J=P(K)
      IF (K.EQ.J) GOTO 22
      DO 20 I=1,N
      T=A(I,J)
      A(I,J)=A(I,K)
  20  A(I,K)=T
  22	I=Q(K)
      IF (K.EQ.I) GOTO 36
      DO 24 J=1,N
      T=A(I,J)
      A(I,J)=A(K,J)
  24  A(K,J)=T
  36  CONTINUE
 100  IF (IFF.NE.0)  STOP 'SINGULAR MATRIX'
      RETURN
      END


