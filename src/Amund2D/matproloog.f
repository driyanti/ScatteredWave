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
      CALL ZGEMM('N','C',dim,dim,dim,dcmplx(1.d0,0.d0),inmat1,dim,
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
