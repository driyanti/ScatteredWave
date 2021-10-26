      
      SUBROUTINE mkD(pabs,matD,layer)

      implicit none
      
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'
c
c     This subroutine computes the matrix D as defined in eqs. 3.36 and
c     3.37 of Brian Kennet's "Seismic wave propagation in stratfied media"
c
      integer	   layer
      real*8 	   pabs
      complex*16   matD(4,4),ealpha,ebeta, factor
      
      
      factor = sqrt(qalpha(layer))*sqrt(2.0*rho(layer))
      ealpha = 1.d0/factor
      factor = sqrt(qbeta(layer))*sqrt(2.0*rho(layer))
      ebeta  = 1.d0/factor
      
      matD(1,1) = -ci*sqrt(qalpha(layer))/sqrt(2.*rho(layer)) 
      	    	  !-ealpha*ci*qalpha(layer)
      matD(1,2) = ebeta*pabs
      matD(1,3) = -matD(1,1)
      matD(1,4) = matD(1,2)
      matD(2,1) = ealpha*pabs
      matD(2,2) = -ci*sqrt(qbeta(layer))/sqrt(2.*rho(layer))
      	    	  !-ebeta*ci*qbeta(layer)
      matD(2,3) = matD(2,1)
      matD(2,4) = -matD(2,2)

      matD(3,1) = -sqrt(rho(layer)/2.)*
     .	    	  beta(layer)**2/sqrt(qalpha(layer))*
     .	    	  (qbeta(layer)**2 - pabs**2)
      	          !ealpha*rho(layer)*(2.0*beta(layer)**2*pabs**2-1.0)
      matD(3,2) = -ci*sqrt(2.*rho(layer))*sqrt(qbeta(layer))*
     .	    	  	pabs*beta(layer)**2
      matD(3,3) = matD(3,1)
      matD(3,4) = -matD(3,2)
      matD(4,1) = -ci*sqrt(2.*rho(layer))*sqrt(qalpha(layer))*
     .	    	  	pabs*beta(layer)**2     
      matD(4,2) = -sqrt(rho(layer)/2.)*
     .	    	  beta(layer)**2/sqrt(qbeta(layer))*
     .	    	  (qbeta(layer)**2 - pabs**2)
      	    	  !ebeta*rho(layer)*(2.0*beta(layer)**2*pabs**2-1.0)
      matD(4,3) = -matD(4,1)
      matD(4,4) = matD(4,2)

      END SUBROUTINE
        
c     *********************************************
c     Einde subroutine mkD, begin subroutine mkDinv
c     *********************************************

      SUBROUTINE mkDinv(matD,Dinv)

      implicit none

c     This subroutine computes the inverse of matrix D as defined in 
c     eqs. 3.36, 3.37 and 3.40 of Brian Kennet's "Seismic wave propagation 
c     in stratfied media".

      complex*16 matD(4,4),Dinv(4,4)

      complex*16 ci
      ci = cmplx(0.d0,1.d0)

      Dinv(1,1) =-ci*matD(3,3)
      Dinv(1,2) =-ci*matD(4,3)
      Dinv(1,3) = ci*matD(1,3)
      Dinv(1,4) = ci*matD(2,3)
      Dinv(2,1) =-ci*matD(3,4) 
      Dinv(2,2) =-ci*matD(4,4)
      Dinv(2,3) = ci*matD(1,4)
      Dinv(2,4) = ci*matD(2,4)
      Dinv(3,1) = ci*matD(3,1)
      Dinv(3,2) = ci*matD(4,1)
      Dinv(3,3) =-ci*matD(1,1)
      Dinv(3,4) =-ci*matD(2,1)
      Dinv(4,1) = ci*matD(3,2)
      Dinv(4,2) = ci*matD(4,2)
      Dinv(4,3) =-ci*matD(1,2)
      Dinv(4,4) =-ci*matD(2,2)

      END SUBROUTINE
