       SUBROUTINE mkgreen_halfspace(green_halfspace,omega,p1,p2,zrcv,
     .                              zsrc)


      implicit none

      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
      INCLUDE 'embank.h'
      
      real*8      zsrc,zrcv,p1,p2,omega,psum
      complex*16  green_halfspace(dim*dim),dz_green_halfspace(dim*dim),
     .            Kks,DD_plus,DD_min,matPS(dim*dim),matSP(dim*dim),
     .            matP(dim*dim),matS(dim*dim),
     .            matPP(dim*dim),matSS(dim*dim),hlpvar                                  
      alpha_embank = alpha(1)
      beta_embank  = beta(1)
      qalpha_embank = qalpha(1)
      qbeta_embank = qbeta(1)
            
      psum  	    = p1**2 + p2**2
      hlpvar        = sqrt(-psum + 1.d0/alpha_embank**2)
      qalpha_embank = cmplx(abs(dble(hlpvar)),dble(-ci*hlpvar))      
      hlpvar        = sqrt(-psum + 1.d0/beta_embank**2)
      qbeta_embank  = cmplx(abs(dble(hlpvar)),dble(-ci*hlpvar))
      Kks   	    = 1.d0/beta_embank**2 - 2.d0*psum
      DD_plus       = Kks**2 + 4.d0*qalpha_embank*qbeta_embank*psum
      DD_min        = Kks**2 - 4.d0*qalpha_embank*qbeta_embank*psum
      
c
c initialize q_alpha and q_beta with appropriate signs
c sign of q_alpha and q_beta should correspond to the fourier transform
c see De Hoop (waves in elastic media, dissipative medium and our fourier trnsf)
c
      
    
      matP(1) = p1**2   /qalpha_embank/2.d0
     .          *exp(ci*omega*qalpha_embank*abs(zrcv-zsrc))
      matP(2) = p1*p2   /qalpha_embank/2.d0
     .          *exp(ci*omega*qalpha_embank*abs(zrcv-zsrc))
      matP(3) = p1*qalpha_embank*sign(1.d0,zsrc-zrcv)/qalpha_embank/2.d0
     .          *exp(ci*omega*qalpha_embank*abs(zrcv-zsrc))
      matP(4) = matP(2)
      matP(5) = p2**2/qalpha_embank/2.d0
     .          *exp(ci*omega*qalpha_embank*abs(zrcv-zsrc))
      matP(6) = p2*qalpha_embank*sign(1.d0,zsrc-zrcv)/qalpha_embank/2.d0
     .          *exp(ci*omega*qalpha_embank*abs(zrcv-zsrc))
      matP(7) = matP(3)
      matP(8) = matP(6)
      matP(9) = qalpha_embank**2/qalpha_embank/2.d0
     .          *exp(ci*omega*qalpha_embank*abs(zrcv-zsrc))


            
      matS(1) = (1.d0/beta_embank**2 - p1**2)/qbeta_embank/2.d0
     .          *exp(ci*omega*qbeta_embank*abs(zrcv-zsrc))
      matS(2) = -p1*p2/qbeta_embank/2.d0
     .          *exp(ci*omega*qbeta_embank*abs(zrcv-zsrc))
      matS(3) = -p1*qbeta_embank*sign(1.d0,zsrc-zrcv)/qbeta_embank/2.d0
     .          *exp(ci*omega*qbeta_embank*abs(zrcv-zsrc))
      matS(4) = matS(2)
      matS(5) = (1.d0/beta_embank**2 - p2**2)/qbeta_embank/2.d0
     .          *exp(ci*omega*qbeta_embank*abs(zrcv-zsrc))
      matS(6) = -p2*qbeta_embank*sign(1.d0,zsrc-zrcv)/qbeta_embank/2.d0
     .          *exp(ci*omega*qbeta_embank*abs(zrcv-zsrc))
      matS(7) = matS(3)
      matS(8) = matS(6)
      matS(9) = psum/qbeta_embank/2.d0
     .          *exp(ci*omega*qbeta_embank*abs(zrcv-zsrc))


      matPP(1) = -p1**2*DD_min/(qalpha_embank*DD_plus*2.d0)
     .          *exp(ci*omega*qalpha_embank*(zrcv+zsrc))
      matPP(2) = -p1*p2*DD_min/(qalpha_embank*DD_plus*2.d0)
     .          *exp(ci*omega*qalpha_embank*(zrcv+zsrc))
      matPP(3) = -p1*qalpha_embank*DD_min/(qalpha_embank*DD_plus*2.d0)
     .          *exp(ci*omega*qalpha_embank*(zrcv+zsrc))
      matPP(4) = matPP(2)      
      matPP(5) = -p2**2*DD_min/(qalpha_embank*DD_plus*2.d0)
     .          *exp(ci*omega*qalpha_embank*(zrcv+zsrc))
      matPP(6) = -p2*qalpha_embank*DD_min/(qalpha_embank*DD_plus*2.d0)
     .          *exp(ci*omega*qalpha_embank*(zrcv+zsrc))
      matPP(7) = -matPP(3)
      matPP(8) = -matPP(6)
      matPP(9) = qalpha_embank**2*DD_min/(qalpha_embank*DD_plus*2.d0)
     .          *exp(ci*omega*qalpha_embank*(zrcv+zsrc))
       
      matSS(1) = ((1./beta_embank**2-p1**2)*DD_plus-
     .                          8.*qalpha_embank*qbeta_embank**3*p1**2)
     .*exp(ci*omega*qbeta_embank*(zrcv+zsrc))/qbeta_embank/DD_plus/2.0
      matSS(2) = -p1*p2*(DD_plus+8.*qalpha_embank*qbeta_embank**3)
     .*exp(ci*omega*qbeta_embank*(zrcv+zsrc))/qbeta_embank/DD_plus/2.0       
      matSS(3) = -p1*qbeta_embank*DD_min
     .*exp(ci*omega*qbeta_embank*(zrcv+zsrc))/qbeta_embank/DD_plus/2.0
      
      matSS(4) = matSS(2)
      matSS(5) = ((1./beta_embank**2-p2**2)*DD_plus-
     .                          8.*qalpha_embank*qbeta_embank**3*p2**2)
     .*exp(ci*omega*qbeta_embank*(zrcv+zsrc))/qbeta_embank/DD_plus/2.0
      matSS(6) = -p2*qbeta_embank*DD_min
     .*exp(ci*omega*qbeta_embank*(zrcv+zsrc))/qbeta_embank/DD_plus/2.0
      
      matSS(7) = -matSS(3)
      matSS(8) = -matSS(6)
      
      matSS(9) = -psum*DD_min
     .*exp(ci*omega*qbeta_embank*(zrcv+zsrc))/qbeta_embank/DD_plus/2.d0

      matPS(1) = p1**2*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      matPS(2) = p1*p2*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      matPS(3) = p1*qalpha_embank*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      matPS(4) = matPS(2)
      matPS(5) = p2**2*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      matPS(6) = p2*qalpha_embank*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      matPS(7) = p1*psum*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      matPS(8) = p2*psum*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      matPS(9) = qalpha_embank*psum*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zrcv+qalpha_embank*zsrc))
      
      matSP(1) = p1**2*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
      matSP(2) = p1*p2*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
      matSP(3) = -p1*psum*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
      matSP(4) = matSP(2)
      matSP(5) = p2**2*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
      matSP(6) = -p2*psum*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
      matSP(7) = -p1*qalpha_embank*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
      matSP(8) = -p2*qalpha_embank*qbeta_embank*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
      matSP(9) = qalpha_embank*psum*2.*Kks/DD_plus
     .          *exp(ci*omega*(qbeta_embank*zsrc+qalpha_embank*zrcv))
!
! direction x3: extra minus sign, per index 
! 13 = -
! 23 = -
! 33 = - * - = +
!
      green_halfspace(1) = -ci*(matP(1)+matS(1)+matPP(1)+
     .          matSS(1)+matPS(1)+matSP(1))/rho(1)/omega
      green_halfspace(2) = -ci*(matP(2)+matS(2)+matPP(2)+
     .          matSS(2)+matPS(2)+matSP(2))/rho(1)/omega
      green_halfspace(3) = -ci*(matP(3)+matS(3)+matPP(3)+
     .          matSS(3)+matPS(3)+matSP(3))/rho(1)/omega
      green_halfspace(4) = -ci*(matP(4)+matS(4)+matPP(4)+
     .          matSS(4)+matPS(4)+matSP(4))/rho(1)/omega
      green_halfspace(5) = -ci*(matP(5)+matS(5)+matPP(5)+
     .          matSS(5)+matPS(5)+matSP(5))/rho(1)/omega
      green_halfspace(6) = -ci*(matP(6)+matS(6)+matPP(6)+
     .          matSS(6)+matPS(6)+matSP(6))/rho(1)/omega
      green_halfspace(7) = -ci*(matP(7)+matS(7)+matPP(7)+
     .          matSS(7)+matPS(7)+matSP(7))/rho(1)/omega
      green_halfspace(8) = -ci*(matP(8)+matS(8)+matPP(8)+
     .          matSS(8)+matPS(8)+matSP(8))/rho(1)/omega
      green_halfspace(9) = -ci*(matP(9)+matS(9)+matPP(9)+
     .          matSS(9)+matPS(9)+matSP(9))/rho(1)/omega
!
! direction x3: extra minus sign, per index 
! 13 = -
! 23 = -
! 33 = - * - = +
! plus extra - due to differentiation with respect to x3
!
      dz_green_halfspace(1) = (
     .          qalpha_embank*sign(1.d0,zrcv-zsrc)*matP(1)+
     .           qbeta_embank*sign(1.d0,zrcv-zsrc)*matS(1)+
     .          qalpha_embank*matPP(1)+
     .           qbeta_embank*matSS(1)+
     .           qbeta_embank*matPS(1)+
     .          qalpha_embank*matSP(1))/rho(1)
      dz_green_halfspace(2) = (
     .          qalpha_embank*sign(1.d0,zrcv-zsrc)*matP(2)+
     .           qbeta_embank*sign(1.d0,zrcv-zsrc)*matS(2)+
     .          qalpha_embank*matPP(2)+
     .           qbeta_embank*matSS(2)+
     .           qbeta_embank*matPS(2)+
     .          qalpha_embank*matSP(2))/rho(1)
      dz_green_halfspace(3) = -(
     .          qalpha_embank*sign(1.d0,zrcv-zsrc)*matP(3)+
     .           qbeta_embank*sign(1.d0,zrcv-zsrc)*matS(3)+
     .          qalpha_embank*matPP(3)+
     .           qbeta_embank*matSS(3)+
     .           qbeta_embank*matPS(3)+
     .          qalpha_embank*matSP(3))/rho(1)
      dz_green_halfspace(4) = (
     .          qalpha_embank*sign(1.d0,zrcv-zsrc)*matP(4)+
     .           qbeta_embank*sign(1.d0,zrcv-zsrc)*matS(4)+
     .          qalpha_embank*matPP(4)+
     .           qbeta_embank*matSS(4)+
     .           qbeta_embank*matPS(4)+
     .          qalpha_embank*matSP(4))/rho(1)
      dz_green_halfspace(5) = -(
     .          qalpha_embank*sign(1.d0,zrcv-zsrc)*matP(5)+
     .           qbeta_embank*sign(1.d0,zrcv-zsrc)*matS(5)+
     .          qalpha_embank*matPP(5)+
     .           qbeta_embank*matSS(5)+
     .           qbeta_embank*matPS(5)+
     .          qalpha_embank*matSP(5))/rho(1)
      dz_green_halfspace(6) = (
     .          qalpha_embank*sign(1.d0,zrcv-zsrc)*matP(6)+
     .           qbeta_embank*sign(1.d0,zrcv-zsrc)*matS(6)+
     .          qalpha_embank*matPP(6)+
     .           qbeta_embank*matSS(6)+
     .           qbeta_embank*matPS(6)+
     .          qalpha_embank*matSP(6))/rho(1)
                         
      end
