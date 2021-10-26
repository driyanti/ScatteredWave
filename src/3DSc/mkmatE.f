      SUBROUTINE mkE(omega,deltaz,E,layer)

      implicit none
      
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'
      
c     This subroutine computes the 2x2 propagation matrix E

      integer	 layer
      real*8     omega,deltaz,check1,check2,check,maxexp
      complex*16 E(3,3),c0
      
      check1 = -dble(ci*qalpha(layer)*omega)*abs(deltaz)
      check2 = -dble(ci*qbeta(layer)*omega)*abs(deltaz)  
      check  = min(check1,check2)
                  
      maxexp = 30.0
      
      c0 = cmplx(0.d0,0.d0)
      
      E(1,1) = c0
      E(1,2) = c0
      E(1,3) = c0
      E(2,1) = c0
      E(2,2) = c0
      E(2,3) = c0
      E(3,1) = c0
      E(3,2) = c0
      E(3,3) = c0

      IF (deltaz.eq.0.0) THEN
        E(1,1) = cmplx(1.d0,0.d0)
        E(2,2) = cmplx(1.d0,0.d0)
        E(3,3) = cmplx(1.d0,0.d0)
      ELSE
        IF (check.lt.maxexp) THEN
          E(1,1) = exp(ci*omega*qalpha(layer)*deltaz)
          E(2,2) = exp(ci*omega*qbeta(layer)*deltaz)
          E(3,3) = E(2,2)
        ENDIF
      ENDIF

      END SUBROUTINE
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
