      SUBROUTINE initial
      
      implicit none

      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
      INCLUDE 'record.h'
      
c
c     local variables
c

      integer layer
      real*8    dumalpha,dumbeta,dumral
      
      ci = dcmplx(0.d0,1.d0)
      
      IF (nsamp.eq.2) THEN
        logns = 1
      ELSE IF (nsamp.eq.4) THEN
        logns = 2
      ELSE IF (nsamp.eq.8) THEN
        logns = 3
      ELSE IF (nsamp.eq.16) THEN
        logns = 4
      ELSE IF (nsamp.eq.32) THEN
        logns = 5
      ELSE IF (nsamp.eq.64) THEN
        logns = 6
      ELSE IF (nsamp.eq.128) THEN
        logns = 7
      ELSE IF (nsamp.eq.256) THEN
        logns = 8
      ELSE IF (nsamp.eq.512) THEN
        logns = 9
      ELSE IF (nsamp.eq.1024) THEN
        logns = 10
      ELSE IF (nsamp.eq.2048) THEN
        logns = 11
      ELSE IF (nsamp.eq.4096) THEN
        logns = 12
      ELSE IF (nsamp.eq.8192) THEN
        logns = 13
c      ELSE IF (nsamp.eq.16384) THEN
c        logns = 14
c      ELSE IF (nsamp.eq.32768) THEN
c        logns = 15
c      ELSE IF (nsamp.eq.65536) THEN
c        logns = 16
c      ELSE IF (nsamp.eq.131072) THEN
c        logns = 17
c      ELSE IF (nsamp.eq.262144) THEN
c        logns = 18
c      ELSE IF (nsamp.eq.524288) THEN
c        logns = 19
c      ELSE IF (nsamp.eq.1048576) THEN
c        logns = 20
      ELSE
        write(*,*) '(input.in) nsamp should be a power of 2. '
        write(*,*) 'Your input was: nsamp = ',nsamp
      END IF

      IF (nsamp.gt.maxtsamp) STOP 'Too many samples. (Init.f)'
c
c test for dt and omega domain
c
      IF (fmax .gt. (nsamp/2-1)/(nsamp*dt)) THEN
        write(*,*) 'Error: violation omega domain'
        STOP
      END IF
c
c test for validity of p-domain
c      
      IF ((int(tot_inc_pmin/tot_inc_deltap)+1).gt.maxksamp) THEN
        write(*,*)'p_max too large or delta_p too small(Init.f)'
        write(*,*)'nk_total too small, param.h'
	STOP 
      END IF

cccccccccccccc

      DO layer = 1,nlay+1                 
c        G(layer) = G(layer) - ci*2.d0*D(layer)/100.d0*dble(G(layer))
c        beta(layer) = sqrt(G(layer)/rho(layer))
        beta(layer) = G(layer)*cmplx(1d0,-D(layer))
        beta(layer) = abs(dble(beta(layer))) - 
     1			ci*abs(dble(-ci*beta(layer)))
       if (debug .eq. 1) then
           write(10,*) 'beta(layer) = ', beta(layer)
           write(10,*) 'rho(layer) = ', rho(layer)
       endif
        alpha(layer) = K(layer)*cmplx(1d0,-D(layer)) 
c        alpha(layer) = sqrt((K(layer)+4.d0/3.d0*G(layer))/rho(layer))
        alpha(layer) = abs(dble(alpha(layer))) - 
     1			ci*abs(dble(-ci*alpha(layer)))
        if (debug .eq. 1) then
           write(10,*) 'alpha(layer) = ', alpha(layer)
        endif
        
cc-new----
c      DO layer = 1,nlay+1                 
c        G(layer) = G(layer) - ci*2.d0*D(layer)/100.d0*dble(G(layer))
cc        beta(layer) = sqrt(G(layer)/rho(layer))
c        beta(layer) = abs(dble(beta(layer))) - 
c     1			ci*abs(dble(-ci*beta(layer))) 
c        write(10,*) 'beta(layer) = ', beta(layer)
c      
c        alpha(layer) = sqrt((K(layer)+4.d0/3.d0*G(layer))/rho(layer))
c        alpha(layer) = abs(dble(alpha(layer))) - 
c     1			ci*abs(dble(-ci*alpha(layer)))
c        write(10,*) 'alpha(layer) = ', alpha(layer)
cccccccccccccAuke's versioncccccccccc
        dumalpha = dble(alpha(layer))
        dumbeta  = dble(beta(layer))      
      
        CALL calcr(dumalpha,dumbeta,dumral)
        if (debug .eq. 1) then
          write(10,*) 'rayleigh(layer) = ', dumral
        endif
      END DO
      domega = 2.d0*pi/(dt*nsamp)
c
c determination layer position
c      
      z(0) = 0.d0
      DO layer=1,nlay
        z(layer) = z(layer-1) + h(layer)
      END DO
      z(nlay+1) = z(nlay) + 1000.
c
c determination of srclayer and rcvlayer
c
      DO layer=1,nlay+1
        IF (rcvpos(3).lt.z(layer)) THEN
	  rcvlayer = layer
	  IF (rcvpos(3).eq.z(layer-1)) write(*,*) 'receiver at layer'
	  GOTO 401
	END IF
      END DO
401   IF (rcvlayer.eq.nlay+1) write(*,*) 'receiver in half space'
      
      DO layer=1,nlay+1
        IF (srcpos(3).lt.z(layer)) THEN
	  srclayer = layer
	  IF (srcpos(3).eq.z(layer-1)) write(*,*)'source at layer'
	  GOTO 402
	END IF
      END DO
402   IF (srclayer.eq.nlay+1) write(*,*) 'source in half space'

      END SUBROUTINE
c
c-------------------------------------------------------------------------------
c      
      subroutine CALCR(cp,cs,cr)

      real*8 cp,cs,cr,x0,eta,xk,fxk,ff,fac,xh

      x0=0.0
      eta=(cs/cp)**2
      xk=ff(x0,eta)
      fxk=xk
2020  if (.not. (abs(fxk) .gt. 1E-5)) go to 2030
      xh=ff(xk,eta)/fac(xk,eta)
      xk=xk-xh
      fxk=ff(xk,eta)
      go to 2020
2030  if ((xk .lt. 0.) .or. (xk .gt. 1.)) then
      write(*,*) 'Error in calculation of cR'
      endif
      cr=cs*sqrt(xk)
      end
c
c-------------------------------------------------------------------------------
c
      real*8 function ff(x,a)

      real*8 x,a
      ff=x**3 - 8.*x**2 + 8.*x*(3.-2.*a) + 16.*(a-1.)

      end function
c
c-------------------------------------------------------------------------------      
c
      real*8 function fac(x,a)

      real*8 x,a
      fac=3.*x**2 - 16.*x + 8.*(3.-2.*a)

      end function
