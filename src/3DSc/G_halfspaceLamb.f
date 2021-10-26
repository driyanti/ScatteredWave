      SUBROUTINE mkgreen_halfspace(green_halfspace,omega,p1,p2
     . ,zrcv,zsrc)

cSUBROUTINE mkgreen_halfspace(green_halfspace,
c     .			omega,layer,p1,p2,zrcv,zsrc,sinc_char)

      implicit none

      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
c      INCLUDE 'embank.h'
      
      character   sinc_char
      integer     layer
      real*8      zsrc,zrcv,p1,p2,omega,psum,zdist
      complex*16  green_halfspace(dim,dim),c1,Kks,DD_plus,DD_min
      complex*16  matP(2*dim),matS(2*dim),matPP(2*dim),matSS(2*dim),
     .		  matPS(2*dim),matSP(2*dim)
                                                       
      c1=cmplx(0.0,1.0)     
                        
      psum  =  p1**2 + p2**2
      Kks   =  1.d0/beta(1)**2 - 2.d0*psum
      DD_plus = Kks**2 + 4.d0*qalpha(1)*qbeta(1)*psum
      DD_min  = Kks**2 - 4.d0*qalpha(1)*qbeta(1)*psum
      
      zdist = zrcv-zsrc
      
      matP(1) = p1**2	/qalpha(1)/2.d0
     .		*exp(ci*omega*qalpha(1)*abs(zrcv-zsrc))
      matP(2) = p1*p2	/qalpha(1)/2.d0
     .		*exp(ci*omega*qalpha(1)*abs(zrcv-zsrc))
      matP(3) = p1*qalpha(1)*sign(1,zsrc-zrcv)/qalpha(1)/2.d0
     .		*exp(ci*omega*qalpha(1)*abs(zrcv-zsrc))
      matP(4) = p2**2/qalpha(1)/2.d0
     .		*exp(ci*omega*qalpha(1)*abs(zrcv-zsrc))
      matP(5) = p2*qalpha(1)*sign(1,zsrc-zrcv)/qalpha(1)/2.d0
     .		*exp(ci*omega*qalpha(1)*abs(zrcv-zsrc))
      matP(6) = qalpha(1)**2/qalpha(1)/2.d0
     .		*exp(ci*omega*qalpha(1)*abs(zrcv-zsrc))
      
      matS(1) = (1.d0/beta(1)**2 - p1**2)/qbeta(1)/2.d0
     .		*exp(ci*omega*qbeta(1)*abs(zrcv-zsrc))
      matS(2) = -p1*p2/qbeta(1)/2.d0
     .		*exp(ci*omega*qbeta(1)*abs(zrcv-zsrc))
      matS(3) = -p1*qbeta(1)*sign(1,zsrc-zrcv)/qbeta(1)/2.d0
     .		*exp(ci*omega*qbeta(1)*abs(zrcv-zsrc))
      matS(4) = (1.d0/beta(1)**2 - p2**2)/qbeta(1)/2.d0
     .		*exp(ci*omega*qbeta(1)*abs(zrcv-zsrc))
      matS(5) = -p2*qbeta(1)*sign(1,zsrc-zrcv)/qbeta(1)/2.d0
     .		*exp(ci*omega*qbeta(1)*abs(zrcv-zsrc))
      matS(6) = psum/qbeta(1)/2.d0
     .		*exp(ci*omega*qbeta(1)*abs(zrcv-zsrc))

      matPP(1) = -p1**2*DD_min/(qalpha(1)*DD_plus*2.d0)
     .		*exp(ci*omega*qalpha(1)*(zrcv+zsrc))
      matPP(2) = -p1*p2*DD_min/(qalpha(1)*DD_plus*2.d0)
     .		*exp(ci*omega*qalpha(1)*(zrcv+zsrc))
      matPP(3) = -p1*qalpha(1)*DD_min/(qalpha(1)*DD_plus*2.d0)
     .		*exp(ci*omega*qalpha(1)*(zrcv+zsrc))
      matPP(4) = -p2**2*DD_min/(qalpha(1)*DD_plus*2.d0)
     .		*exp(ci*omega*qalpha(1)*(zrcv+zsrc))
      matPP(5) = -p2*qalpha(1)*DD_min/(qalpha(1)*DD_plus*2.d0)
     .		*exp(ci*omega*qalpha(1)*(zrcv+zsrc))
      matPP(6) = qalpha(1)**2*DD_min/(qalpha(1)*DD_plus*2.d0)
     .		*exp(ci*omega*qalpha(1)*(zrcv+zsrc))
       
      matSS(1) = ((1./beta(1)**2-p1**2)*DD_plus-
     .				8.*qalpha(1)*qbeta(1)**3*p1**2)
     .		*exp(ci*omega*qbeta(1)*(zrcv+zsrc))/qbeta(1)/DD_plus/2.d0
      matSS(2) = -p1*p2*(DD_plus+8.*qalpha(1)*qbeta(1)**3)
     .		*exp(ci*omega*qbeta(1)*(zrcv+zsrc))/qbeta(1)/DD_plus/2.d0
      matSS(3) = -p1*qbeta(1)*DD_min
     .		*exp(ci*omega*qbeta(1)*(zrcv+zsrc))/qbeta(1)/DD_plus/2.d0
      matSS(4) = ((1./beta(1)**2-p2**2)*DD_plus-
     .				8.*qalpha(1)*qbeta(1)**3*p2**2)
     .		*exp(ci*omega*qbeta(1)*(zrcv+zsrc))/qbeta(1)/DD_plus/2.d0
      matSS(5) = -p2*qbeta(1)*DD_min
     .		*exp(ci*omega*qbeta(1)*(zrcv+zsrc))/qbeta(1)/DD_plus/2.d0
      matSS(6) = -psum*DD_min
     .		*exp(ci*omega*qbeta(1)*(zrcv+zsrc))/qbeta(1)/DD_plus/2.d0

      matPS(1) = p1**2*qbeta(1)*2.*Kks/DD_plus
     .		*exp(ci*omega*(qbeta(1)*zrcv+qalpha(1)*zsrc))
      matPS(2) = p1*p2*qbeta(1)*2.*Kks/DD_plus
     .		*exp(ci*omega*(qbeta(1)*zrcv+qalpha(1)*zsrc))
      matPS(3) = p1*qalpha(1)*qbeta(1)*2.*Kks/DD_plus
     .		*exp(ci*omega*(qbeta(1)*zrcv+qalpha(1)*zsrc))
      matPS(4) = p2**2*qbeta(1)*2.*Kks/DD_plus
     .		*exp(ci*omega*(qbeta(1)*zrcv+qalpha(1)*zsrc))
      matPS(5) = p2*qalpha(1)*qbeta(1)*2.*Kks/DD_plus
     .		*exp(ci*omega*(qbeta(1)*zrcv+qalpha(1)*zsrc))
      matPS(6) = qalpha(1)*psum*2.*Kks/DD_plus
     .		*exp(ci*omega*(qbeta(1)*zrcv+qalpha(1)*zsrc))
      
      matSP(1) = p1**2*qbeta(1)*2.*Kks/DD_plus
     .			*exp(ci*omega*(qbeta(1)*zsrc+qalpha(1)*zrcv))
      matSP(2) = p1*p2*qbeta(1)*2.*Kks/DD_plus
     .			*exp(ci*omega*(qbeta(1)*zsrc+qalpha(1)*zrcv))
      matSP(3) = -p1*psum*2.*Kks/DD_plus
     .			*exp(ci*omega*(qbeta(1)*zsrc+qalpha(1)*zrcv))
      matSP(4) = p2**2*qbeta(1)*2.*Kks/DD_plus
     .			*exp(ci*omega*(qbeta(1)*zsrc+qalpha(1)*zrcv))
      matSP(5) = -p2*psum*2.*Kks/DD_plus
     .			*exp(ci*omega*(qbeta(1)*zsrc+qalpha(1)*zrcv))
      matSP(6) = qalpha(1)*psum*2.*Kks/DD_plus
     .			*exp(ci*omega*(qbeta(1)*zsrc+qalpha(1)*zrcv))

      green_halfspace(1,1) = ci*(matP(1)+matS(1)+matPP(1)+
     .		matSS(1)+matPS(1)+matSP(1))/rho(1)/omega
      green_halfspace(2,1) = ci*(matP(2)+matS(2)+matPP(2)+
     .		matSS(2)+matPS(2)+matSP(2))/rho(1)/omega
      green_halfspace(3,1) = -ci*(matP(3)+matS(3)+matPP(3)+
     .		matSS(3)+matPS(3)+matSP(3))/rho(1)/omega
      green_halfspace(1,2) = green_halfspace(2,1)
      green_halfspace(2,2) = ci*(matP(4)+matS(4)+matPP(4)+
     .		matSS(4)+matPS(4)+matSP(4))/rho(1)/omega
      green_halfspace(3,2) = -ci*(matP(5)+matS(5)+matPP(5)+
     .		matSS(5)+matPS(5)+matSP(5))/rho(1)/omega
      green_halfspace(1,3) = green_halfspace(3,1)
      green_halfspace(2,3) = green_halfspace(3,2)
      green_halfspace(3,3) = ci*(matP(6)+matS(6)+matPP(6)+
     .		matSS(6)+matPS(6)+matSP(6))/rho(1)/omega
                         
      end subroutine
