
      SUBROUTINE mkDH(DH,layer)
      
      implicit none
      
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
c
c     This subroutine computes the matrix D as defined in eqs. 3.36 and
c     3.38 of Brian Kennet's "Seismic wave propagation in stratfied media"
c
      integer 	layer
      complex*16 	DH(2,2),ebeta,ci,factor
      
      ci 	    = cmplx(0.d0,1.d0)
      factor    = sqrt(2.d0*rho(layer)*qbeta(layer))
      ebeta     = 1.d0/factor

      DH(1,1) = 1.d0/(factor*beta(layer))
      DH(1,2) = 1.d0/(factor*beta(layer))
      DH(2,1) =-ci*beta(layer)*sqrt(rho(layer)*qbeta(layer))/sqrt(2.d0)
      DH(2,2) = ci*beta(layer)*sqrt(rho(layer)*qbeta(layer))/sqrt(2.d0)

      END SUBROUTINE
        
c----------------------------------------------------------------------
      SUBROUTINE mkDinvh(DH,DHinv)

      implicit none

c     This subroutine computes the inverse of matrix Dh as defined in 
c     eqs. 3.36, 3.38 and 3.40 of Brian Kennet's "Seismic wave propagation 
c     in stratfied media".


      complex*16 DH(2,2),DHinv(2,2),ci

      ci = cmplx(0.d0,1.d0)

      DHinv(1,1) =-ci*DH(2,2)
      DHinv(1,2) = ci*DH(1,2)
      DHinv(2,1) = ci*DH(2,1) 
      DHinv(2,2) =-ci*DH(1,1)

      END SUBROUTINE
      
