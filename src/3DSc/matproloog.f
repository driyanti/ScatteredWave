      SUBROUTINE matrix_proloog(pabs,omega,Su,Sd,srclevel)

      implicit none
      
      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
!
! pabs = sqrt(p1^2 + p2^2)
!      
      integer	layer,ii,jj
      real*8	omega,srclevel,pabs
      complex*16 Td(3,3),Ru(3,3),Rd(3,3),Tu(3,3),
     1		Su(3,3),Sd(3,3),c0
      complex*16  Rfs(2,2)
      
      c0 = dcmplx(0.d0,0.d0)

      DO ii=1,3
        DO jj=1,3
          Td(ii,jj) = c0
	  Tu(ii,jj) = c0
	  Rd(ii,jj) = c0
	  Ru(ii,jj) = c0
        END DO
      END DO

      CALL mkRfs(pabs,Rfs)
      
      matRu(1,1,0) = Rfs(1,1)
      matRu(1,2,0) = Rfs(1,2)
      matRu(1,3,0) = c0
      matRu(2,1,0) = Rfs(2,1)
      matRu(2,2,0) = Rfs(2,2)
      matRu(2,3,0) = c0
      matRu(3,1,0) = c0
      matRu(3,2,0) = c0
      matRu(3,3,0) = dcmplx(1.0,0.0)

      DO layer = 1,nlay
        CALL mkTdRuRdTu(pabs,Td,Ru,Rd,Tu,layer)
        CALL TWODtoTHREED(Td,matTd,layer)
        CALL TWODtoTHREED(Ru,matRu,layer)
        CALL TWODtoTHREED(Rd,matRd,layer)
        CALL TWODtoTHREED(Tu,matTu,layer)	
      END DO
      
      CALL mkmatSu(omega,Su,srclevel)      
      CALL mkmatSd(omega,Sd,srclevel)
      
      END SUBROUTINE

c---------------------------------------------------------------------------
c      
c---------------------------------------------------------------------------

      SUBROUTINE mkRfs(pabs,Rfs)
      
      IMPLICIT NONE

      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'
c
c     This subroutine computes the reflection matrix of the free surface
c     for an upgoing wave in a solid.
c
      real*8 	 pabs,psum
      complex*16 Rfs(2,2),denom,v
      
      
      psum = pabs**2
      v = 2.d0*psum - 1.d0/beta(1)**2

!      denom   = 4.d0*psum*dcmplx(qalpha(1))*dcmplx(qbeta(1)) + v**2
      denom   = 4.d0*psum*(psum+qalpha(1)*qbeta(1)-1.d0/beta(1)**2) + 
     .	    	  	1.d0/beta(1)**4
      
      Rfs(1,1) = dcmplx(1.d0,0.d0) - 2.d0*v**2/denom
      Rfs(1,2) = (4.d0*ci*sqrt(psum)*v*sqrt(qalpha(1)*qbeta(1)))/denom
      Rfs(2,1) = Rfs(1,2)
      Rfs(2,2) = Rfs(1,1)

      END SUBROUTINE
c----------------------------------------------------------------------------
c
c----------------------------------------------------------------------------

      SUBROUTINE mkmatSu(omega,matSu,srclevel)
      
      IMPLICIT NONE
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'

      real*8	 omega,srclevel
      complex*16 Emat(3,3),matSu(3,3),hlpmat1(3,3)
c
c one-layer model (layer on top of half space)
c
c      EYE(1,1) = cmplx(1.d0,0.d0)
c      EYE(1,2) = cmplx(0.d0,0.d0)
c      EYE(1,3) = cmplx(0.d0,0.d0)
c      EYE(2,1) = cmplx(0.d0,0.d0)
c      EYE(2,2) = cmplx(1.d0,0.d0)
c      EYE(2,3) = cmplx(0.d0,0.d0)
c      EYE(3,1) = cmplx(0.d0,0.d0)
c      EYE(3,2) = cmplx(0.d0,0.d0)
c      EYE(3,3) = cmplx(1.d0,0.d0)
 
      CALL THREEDtoTWOD(matRu,matSu,srclayer-1)

c      IF (srclayer.gt.1) THEN
c        CALL mkE(omega,h(srclayer-1),Emat,srclayer-1)
c	CALL MxM(Emat,matSu,hlpmat1,3)
c	CALL MxM(hlpmat1,Emat,hlpmat2,3)
c	CALL THREEDtoTWOD(matRd,hlpmat1,srclayer-1)
c	CALL MxM(hlpmat1,hlpmat2,hlpmat3,3)
c        CALL matsub(EYE,hlpmat3,hlpmat1)
c	CALL invert(hlpmat1,hlpmat2,3)
c	CALL THREEDtoTWOD(matTu,hlpmat1,srclayer-1)
c	CALL MxM(hlpmat2,hlpmat1,hlpmat3,3)
c	CALL MxM(Emat,hlpmat3,hlpmat2,3)
c	CALL MxM(matSu,hlpmat2,hlpmat3,3)
c	CALL MxM(Emat,hlpmat3,hlpmat2,3)
c	CALL THREEDtoTWOD(matTd,hlpmat1,srclayer-1)
c	CALL MxM(hlpmat1,hlpmat2,hlpmat3,3)
c	CALL THREEDtoTWOD(matRu,hlpmat1,srclayer-1)
c	CALL summat(hlpmat1,hlpmat3,matSu)
c      END IF

cc      CALL mkE(omega,srcpos(3)-z(srclayer-1),Emat,srclayer)

      CALL mkE(omega,srclevel,Emat,srclayer)

      CALL MxM(Emat,matSu,hlpmat1,3)

      CALL MxM(hlpmat1,Emat,matSu,3) 
      
      END SUBROUTINE

c----------------------------------------------------------------------------
c
c----------------------------------------------------------------------------

      SUBROUTINE mkmatSd(omega,matSd,srclevel)

      IMPLICIT NONE  
          
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'
      
      real*8	 omega,srclevel
      complex*16 Emat(3,3),matSd(3,3),hlpmat1(3,3)
c
c one-layer model
c
      
c      IF (srclayer.eq.nlay+1) THEN
c        matSd(1,1) = cmplx(0.d0,0.d0)
c        matSd(1,2) = cmplx(0.d0,0.d0)
c        matSd(1,3) = cmplx(0.d0,0.d0)
c        matSd(2,1) = cmplx(0.d0,0.d0)
c        matSd(2,2) = cmplx(0.d0,0.d0)
c        matSd(2,3) = cmplx(0.d0,0.d0)
c        matSd(3,1) = cmplx(0.d0,0.d0)
c        matSd(3,2) = cmplx(0.d0,0.d0)
c        matSd(3,3) = cmplx(0.d0,0.d0)
c      ELSE

	CALL THREEDtoTWOD(matRd,matSd,srclayer)

c      END IF

!      CALL mkE(omega,z(srclayer)-srclevel,Emat,srclayer)
      CALL mkE(omega,z(1)-srclevel,Emat,srclayer)
      
      CALL MxM(Emat,matSd,hlpmat1,3)

      CALL MxM(hlpmat1,Emat,matSd,3) 

      END SUBROUTINE

ccccccccccccccccccccccccccccc
      SUBROUTINE THREEDtoTWOD(inmat1,outmat,layer)

      implicit none
      
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      
      integer  layer
      complex*16 inmat1(3,3,0:maxlay),outmat(3,3)
      
      outmat(1,1) = inmat1(1,1,layer)
      outmat(1,2) = inmat1(1,2,layer)
      outmat(1,3) = inmat1(1,3,layer)
      outmat(2,1) = inmat1(2,1,layer)
      outmat(2,2) = inmat1(2,2,layer)
      outmat(2,3) = inmat1(2,3,layer)
      outmat(3,1) = inmat1(3,1,layer)
      outmat(3,2) = inmat1(3,2,layer)
      outmat(3,3) = inmat1(3,3,layer)
      
      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE TWODtoTHREED(inmat,outmat,layer)

      implicit none
      
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'

      integer  layer
      complex*16 inmat(3,3),outmat(3,3,0:maxlay)
      
      outmat(1,1,layer) = inmat(1,1)
      outmat(1,2,layer) = inmat(1,2)
      outmat(1,3,layer) = inmat(1,3)
      outmat(2,1,layer) = inmat(2,1)
      outmat(2,2,layer) = inmat(2,2)
      outmat(2,3,layer) = inmat(2,3)
      outmat(3,1,layer) = inmat(3,1)
      outmat(3,2,layer) = inmat(3,2)
      outmat(3,3,layer) = inmat(3,3)
      
      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE matsub(inmat1,inmat2,outmat)

      implicit none
      
      complex*16 inmat1(3,3),inmat2(3,3),outmat(3,3)
      
      outmat(1,1) = inmat1(1,1) - inmat2(1,1)
      outmat(1,2) = inmat1(1,2) - inmat2(1,2)
      outmat(1,3) = inmat1(1,3) - inmat2(1,3)
      outmat(2,1) = inmat1(2,1) - inmat2(2,1)
      outmat(2,2) = inmat1(2,2) - inmat2(2,2)
      outmat(2,3) = inmat1(2,3) - inmat2(2,3)
      outmat(3,1) = inmat1(3,1) - inmat2(3,1)
      outmat(3,2) = inmat1(3,2) - inmat2(3,2)
      outmat(3,3) = inmat1(3,3) - inmat2(3,3)
      
      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE vecsub(invec1,invec2,outvec)
      implicit none
      
      complex*16 invec1(3),invec2(3),outvec(3)
      
      outvec(1) = invec1(1) - invec2(1) 
      outvec(2) = invec1(2) - invec2(2) 
      outvec(3) = invec1(3) - invec2(3) 
     
      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE build6x6(inmat1,inmat2,outmat)

      implicit none

      complex*16 inmat1(4,4),inmat2(2,2),outmat(6,6),c0
      
      c0 = dcmplx(0.d0,0.d0)

      outmat(1,1) = inmat1(1,1)
      outmat(1,2) = inmat1(1,2)
      outmat(1,3) = c0
      outmat(1,4) = inmat1(1,3)
      outmat(1,5) = inmat1(1,4)
      outmat(1,6) = c0

      outmat(2,1) = inmat1(2,1)
      outmat(2,2) = inmat1(2,2)
      outmat(2,3) = c0
      outmat(2,4) = inmat1(2,3)
      outmat(2,5) = inmat1(2,4)
      outmat(2,6) = c0

      outmat(3,1) = c0
      outmat(3,2) = c0
      outmat(3,3) = inmat2(1,1)
      outmat(3,4) = c0
      outmat(3,5) = c0
      outmat(3,6) = inmat2(1,2)

      outmat(4,1) = inmat1(3,1)
      outmat(4,2) = inmat1(3,2)
      outmat(4,3) = c0
      outmat(4,4) = inmat1(3,3)
      outmat(4,5) = inmat1(3,4)
      outmat(4,6) = c0

      outmat(5,1) = inmat1(4,1)
      outmat(5,2) = inmat1(4,2)
      outmat(5,3) = c0
      outmat(5,4) = inmat1(4,3)
      outmat(5,5) = inmat1(4,4)
      outmat(5,6) = c0

      outmat(6,1) = c0
      outmat(6,2) = c0
      outmat(6,3) = inmat2(2,1)
      outmat(6,4) = c0
      outmat(6,5) = c0
      outmat(6,6) = inmat2(2,2)

      END SUBROUTINE

c=================================================================
c
c=================================================================
      SUBROUTINE MxM(inmat1,inmat2,outmat,dim)
      implicit none

      integer  dim
      complex*16 inmat1(dim,dim),inmat2(dim,dim),outmat(dim,dim)
!
! single precision
!      
!      CALL CGEMM('N','N',dim,dim,dim,cmplx(1.d0,0.d0),inmat1,dim,
!     &		inmat2,dim,cmplx(0.d0,0.d0),outmat,dim)
!
! double precision
!
      CALL ZGEMM('N','N',dim,dim,dim,dcmplx(1.d0,0.d0),inmat1,dim,
     &		inmat2,dim,dcmplx(0.d0,0.d0),outmat,dim)

      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE matxvec(inmat,invec,outvec,dim)

      implicit none

      integer  dim
      complex*16 inmat(dim,dim),invec(dim),outvec(dim)
!
! single precision
!      
!      CALL CGEMV('N',dim,dim,cmplx(1.d0,0.d0),inmat,dim,
!     &		invec,1,cmplx(0.d0,0.d0),outvec,1)
!
! double precision
!
      CALL ZGEMV('N',dim,dim,dcmplx(1.d0,0.d0),inmat,dim,
     &		invec,1,dcmplx(0.d0,0.d0),outvec,1)

      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE summat(inmat1,inmat2,outmat)

      implicit none
      
      complex*16 inmat1(3,3),inmat2(3,3),outmat(3,3)
      
      outmat(1,1) = inmat1(1,1) + inmat2(1,1)
      outmat(1,2) = inmat1(1,2) + inmat2(1,2)
      outmat(1,3) = inmat1(1,3) + inmat2(1,3)
      outmat(2,1) = inmat1(2,1) + inmat2(2,1)
      outmat(2,2) = inmat1(2,2) + inmat2(2,2)
      outmat(2,3) = inmat1(2,3) + inmat2(2,3)
      outmat(3,1) = inmat1(3,1) + inmat2(3,1)
      outmat(3,2) = inmat1(3,2) + inmat2(3,2)
      outmat(3,3) = inmat1(3,3) + inmat2(3,3)
      
      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE submat(inmat1,inmat2,outmat)
      implicit none
      
      complex*16 inmat1(3,3),inmat2(3,3),outmat(3,3)
      
      outmat(1,1) = inmat1(1,1) - inmat2(1,1)
      outmat(1,2) = inmat1(1,2) - inmat2(1,2)
      outmat(1,3) = inmat1(1,3) - inmat2(1,3)
      outmat(2,1) = inmat1(2,1) - inmat2(2,1)
      outmat(2,2) = inmat1(2,2) - inmat2(2,2)
      outmat(2,3) = inmat1(2,3) - inmat2(2,3)
      outmat(3,1) = inmat1(3,1) - inmat2(3,1)
      outmat(3,2) = inmat1(3,2) - inmat2(3,2)
      outmat(3,3) = inmat1(3,3) - inmat2(3,3)
      
      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE invert(inmat,outmat,dim)
      implicit none
      
      integer ii,jj,info,dim,ipiv(dim)
      complex*16 inmat(dim,dim),outmat(dim,dim),work(dim*dim)
      
      DO ii=1,dim
        DO jj=1,dim
	  IF (ii.eq.jj) THEN 
	    outmat(ii,jj)=dcmplx(1.d0,0.d0)
	  ELSE
            outmat(ii,jj)=dcmplx(0.d0,0.d0)
	  END IF
	END DO
      END DO   
      
c
c single precision
c
!      CALL  CGESV(dim,dim,inmat,dim,ipiv,outmat,dim,info)
c
c double precision
c
cc      CALL ZGETRF( dim, dim, inmat, dim, IPIV, INFO )
cc      CALL ZGETRI( dim, inmat, dim, IPIV, WORK, dim*dim, INFO )
cc      CALL  ZGESV(dim,dim,inmat,dim,ipiv,outmat,dim,info)

       IF (info.ne.0) THEN
         STOP 'false inverse operation'
       END IF
      

      END SUBROUTINE
c=================================================================
c
c=================================================================
      SUBROUTINE vecsum(invec1,invec2,outvec)
      implicit none
      
      complex*16 invec1(3),invec2(3),outvec(3)
      
      outvec(1) = invec1(1) + invec2(1) 
      outvec(2) = invec1(2) + invec2(2) 
      outvec(3) = invec1(3) + invec2(3) 
     
      END SUBROUTINE

c     **********************************************
c     begin subroutine split
c     **********************************************
      
      SUBROUTINE split(D,mu,nu,md,nd)

      implicit none

c
c     This subroutine splits up the matrix D into its submatrices:
c         / Mu  Md \
c     D=  |        |
c         \ Nu  Nd /

c     Subroutines used: none
c     This subroutine is called from: mkint

c     Global variables
c     Name      I/O   Type     Description
c     ------------------------------------------------------------------------
c     D          I    c(4,4)   The eigenvector matrix to be split up
c     mu         O    c(2,2)   The upper left 2x2 submatrix of D
c     md         O    c(2,2)   The upper right 2x2 submatrix of D
c     nu         O    c(2,2)   The lower left 2x2 submatrix of D
c     nd         O    c(2,2)   The lower right 2x2 submatrix of D

      complex*16 D(4,4),mu(2,2),nu(2,2),md(2,2),nd(2,2)

      mu(1,1) = D(1,1)
      mu(1,2) = D(1,2)
      mu(2,1) = D(2,1)
      mu(2,2) = D(2,2)
      md(1,1) = D(1,3)
      md(1,2) = D(1,4)
      md(2,1) = D(2,3)
      md(2,2) = D(2,4)
      nu(1,1) = D(3,1)
      nu(1,2) = D(3,2)
      nu(2,1) = D(4,1)
      nu(2,2) = D(4,2)
      nd(1,1) = D(3,3)
      nd(1,2) = D(3,4)
      nd(2,1) = D(4,3)
      nd(2,2) = D(4,4)
      
      END SUBROUTINE
        
c     ***********************************************
c     Einde subroutine split, begin subroutine minvar
c     ***********************************************

      SUBROUTINE minvar(mu,nu,md,nd,invar)
      
      implicit none

c     This subroutine computes a matrix-invariant <M_U,M_D> as defined in 
c     eq. 2.58 of Brian Kennett's book "Seismic wave-propagation in 
c     stratified media" [1983].

c     Subroutines used: transp, Mmin and MxM
c     This subroutine is called from: mkint

c     Global variables
c     Name      I/O   Type     Description
c     ------------------------------------------------------------------------
c     mu         I    c(2,2)   M_U
c     md         I    c(2,2)   M_D
c     nu         I    c(2,2)   N_U
c     nd         I    c(2,2)   N_D
c     invar      O    c(2,2)   The matrix-invariant to be computed

      complex*16 mu(2,2),nu(2,2),md(2,2),nd(2,2),invar(2,2)
      complex*16 mut(2,2),nut(2,2),dum1(2,2),dum2(2,2)

      CALL transp(mu,mut,2)
      CALL transp(nu,nut,2)

      CALL MxM(mut,nd,dum1,2)
      CALL MxM(nut,md,dum2,2)

      CALL Mmin(dum1,dum2,invar)

      END

c=============================================================
c-----------------------------------------------------------------
c=============================================================

      SUBROUTINE transp(A,AT,n)
      implicit none
c
c     This subroutine computes the transpose of a square matrix A.
c     The result is called AT
c
c     Global variables
	integer  n
	complex*16 A(n,n),AT(n,n)

	IF (n.eq.2) THEN
	  AT(1,1) = A(1,1)
	  AT(1,2) = A(2,1)
	  AT(2,1) = A(1,2)
	  AT(2,2) = A(2,2)
	ELSE 
	  AT(1,1) = A(1,1)
	  AT(1,2) = A(2,1)
	  AT(1,3) = A(3,1)
	  AT(2,1) = A(1,2)
	  AT(2,2) = A(2,2)
	  AT(2,3) = A(3,2)
	  AT(3,1) = A(1,3)
	  AT(3,2) = A(2,3)
	  AT(3,3) = A(3,3)
	END IF
	END SUBROUTINE
c----------------------------------------------------------------
      SUBROUTINE Mmin(inmat1,inmat2,outmat)
      implicit none
      
      complex*16 inmat1(2,2),inmat2(2,2),outmat(2,2)
      
      outmat(1,1) = inmat1(1,1) - inmat2(1,1)
      outmat(1,2) = inmat1(1,2) - inmat2(1,2)
      outmat(2,1) = inmat1(2,1) - inmat2(2,1)
      outmat(2,2) = inmat1(2,2) - inmat2(2,2)
      
      END SUBROUTINE
c=================================================================
      SUBROUTINE TransposeMat(A,n)
      implicit none
c
c     This subroutine computes the transpose of a square matrix A.
c     The result is called AT
c
c     Global variables
	integer n,i,j
	complex*16 A(n,n),AT(n,n)

	  DO i=1,n
	    DO j=1,n
	     AT(i,j) = A(j,i)
	    END DO
	  END DO
	  DO i=1,n
	    DO j=1,n
	     A(i,j) = AT(i,j)
	    END DO
	  END DO	  

	END SUBROUTINE

ccccccccccccccccccccccccccccccccccccccccccc
     
      SUBROUTINE mkTdRuRdTu(pabs,Td,Ru,Rd,Tu,layer)

      implicit none
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'
      
      integer	layer
      real*8   	pabs

      complex*16 Rd(3,3),Ru(3,3),Td(3,3),Tu(3,3)
      complex*16 Rd_sh,Ru_sh,Td_sh,Tu_sh
      complex*16 D1(4,4),D2(4,4),invar1(2,2),invar2(2,2),invar3(2,2),
     1		dum1(2,2),mu1(2,2),md1(2,2),nu1(2,2),nd1(2,2),
     2        	mu2(2,2),md2(2,2),nu2(2,2),nd2(2,2)

c
c  SH-part
c     
      CALL mkmatSH(Td_sh,Ru_sh,Rd_sh,Tu_sh,layer)
c
c  P-SV
c
      IF ((alpha(layer).eq.alpha(layer+1)).and.
     .	    	  (beta(layer).eq.beta(layer+1)) ) THEN
       
      Td(1,1) = dcmplx(1.d0,0.d0)
      Td(1,2) = dcmplx(0.d0,0.d0)
      Td(1,3) = dcmplx(0.d0,0.d0)
      Td(2,1) = dcmplx(0.d0,0.d0)
      Td(2,2) = dcmplx(1.d0,0.d0)
      Td(2,3) = dcmplx(0.d0,0.d0)
      Td(3,1) = dcmplx(0.d0,0.d0)
      Td(3,2) = dcmplx(0.d0,0.d0)
      Td(3,3) = dcmplx(1.d0,0.d0)

      Rd(1,1) = dcmplx(0.d0,0.d0)
      Rd(1,2) = dcmplx(0.d0,0.d0)
      Rd(1,3) = dcmplx(0.d0,0.d0)
      Rd(2,1) = dcmplx(0.d0,0.d0)
      Rd(2,2) = dcmplx(0.d0,0.d0)
      Rd(2,3) = dcmplx(0.d0,0.d0)
      Rd(3,1) = dcmplx(0.d0,0.d0)
      Rd(3,2) = dcmplx(0.d0,0.d0)
      Rd(3,3) = dcmplx(0.d0,0.d0)

      Ru(1,1) = dcmplx(0.d0,0.d0)
      Ru(1,2) = dcmplx(0.d0,0.d0)
      Ru(1,3) = dcmplx(0.d0,0.d0)
      Ru(2,1) = dcmplx(0.d0,0.d0)
      Ru(2,2) = dcmplx(0.d0,0.d0)
      Ru(2,3) = dcmplx(0.d0,0.d0)
      Ru(3,1) = dcmplx(0.d0,0.d0)
      Ru(3,2) = dcmplx(0.d0,0.d0)
      Ru(3,3) = dcmplx(0.d0,0.d0)

      ELSE
      
      CALL mkD(pabs,D1,layer)
      CALL split(D1,mu1,nu1,md1,nd1)
      CALL mkD(pabs,D2,layer+1)

      CALL split(D2,mu2,nu2,md2,nd2)
      CALL minvar(md1,nd1,md2,nd2,invar1)
      CALL minvar(mu1,nu1,md2,nd2,dum1)
      CALL invert(dum1,invar2,2)
      CALL minvar(mu1,nu1,mu2,nu2,invar3)

      CALL MxM(invar1,invar2,dum1,2)

      Rd(1,1) = -dum1(1,1)
      Rd(1,2) = -dum1(1,2)
      Rd(1,3) = dcmplx(0.d0,0.d0)
      Rd(2,1) = -dum1(2,1)
      Rd(2,2) = -dum1(2,2)
      Rd(2,3) = dcmplx(0.d0,0.d0)
      Rd(3,1) = dcmplx(0.d0,0.d0)
      Rd(3,2) = dcmplx(0.d0,0.d0)
      Rd(3,3) = Rd_sh

      Td(1,1) = ci*invar2(1,1)
      Td(1,2) = ci*invar2(1,2)
      Td(1,3) = dcmplx(0.d0,0.d0)
      Td(2,1) = ci*invar2(2,1)
      Td(2,2) = ci*invar2(2,2)
      Td(2,3) = dcmplx(0.d0,0.d0)
      Td(3,1) = dcmplx(0.d0,0.d0)
      Td(3,2) = dcmplx(0.d0,0.d0)
      Td(3,3) = Td_sh

      CALL MxM(invar2,invar3,dum1,2)

      Ru(1,1) = -dum1(1,1)
      Ru(1,2) = -dum1(1,2)
      Ru(1,3) = dcmplx(0.d0,0.d0)
      Ru(2,1) = -dum1(2,1)
      Ru(2,2) = -dum1(2,2)
      Ru(2,3) = dcmplx(0.d0,0.d0)
      Ru(3,1) = dcmplx(0.d0,0.d0)
      Ru(3,2) = dcmplx(0.d0,0.d0)
      Ru(3,3) = Ru_sh
      
      ENDIF
      
      CALL transp(Td,Tu,3)
      
      END SUBROUTINE
c------------------------------------------------------------------
c
c------------------------------------------------------------------
      SUBROUTINE mkmatSH(Td,Ru,Rd,Tu,layer)

      implicit none

      INCLUDE 'bounds.h'
      INCLUDE 'param.h'

      integer	 layer
      complex*16 Rd,Ru,Td,Tu

	Rd = (rho(layer)*beta(layer)**2*qbeta(layer)-
     1			rho(layer+1)*beta(layer+1)**2*qbeta(layer+1))/
     2		(rho(layer)*beta(layer)**2*qbeta(layer)+
     3     		rho(layer+1)*beta(layer+1)**2*qbeta(layer+1))
	Td = 2.d0*sqrt(rho(layer)*beta(layer)**2*
     1			rho(layer+1)*beta(layer+1)**2*
     1				qbeta(layer)*qbeta(layer+1))/
     2		     (rho(layer)*beta(layer)**2*qbeta(layer)+
     3			rho(layer+1)*beta(layer+1)**2*qbeta(layer+1))
	Ru = -1.d0*Rd 
	Tu = Td		
      END SUBROUTINE

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
