      SUBROUTINE mkpdom(omega,pommax,pommin)
      
c
c author: Auke Ditzel
c year:   1999
c
c subroutine to calculate p_y (integration parameter)
c cosine tapering function is applied on the Greens function (slowness domain)
c to get rid of edge effects
c size of pdom depends on omega (frequency)
c
c next step: Filon integration module to speed up integration procedure
c
      implicit none

      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'

      real*8      p_y(2*maxksamp+1),ptaper(2*maxksamp+1),
     1            omegamax,omegamin,omega,p_max,factor,pommax,pommin
      common/pvalues/p_y,ptaper

      integer counter

      factor   = 1.d0
      omegamin = (int(2.d0*pi*fmin/domega)+1)*domega
      omegamax = (int(2.d0*pi*fmax/domega)+2)*domega
      
      p_max =   pommin + 
     1  (pommin-pommax)/(omegamax-omegamin)*(omegamin-omega)

      nk_totpos = int(p_max / deltap) + 1
      nk_total = 2*nk_totpos + 1
      
      DO counter=1,nk_total
         p_y(counter) = dble(counter-nk_totpos-1)*deltap
      END DO
            
      DO counter=1,nk_total
         ptaper(counter) = factor
      END DO

      DO counter=1,int(nk_total/16.d0)
        ptaper(counter) = 0.d0
      END DO      

      DO counter=int(15.d0/16.d0*nk_total)+1,nk_total
        ptaper(counter) = 0.d0
      END DO

      DO counter=int(6.d0/8.d0*nk_total)+1,int(15.d0/16.d0*nk_total)
         ptaper(counter) =factor* 
     1   (cos(
     2    dble(counter-6.d0/8.d0*nk_total-1)/
     3    dble(nk_total*15.d0/16.d0-6.d0/8.d0*nk_total-1)*
     4    pi/2.d0))**2
      END DO
      
      DO counter=int(nk_total/16.d0)+1,int(nk_total/8.d0)
        ptaper(counter) =factor* 
     1   (sin(dble(counter-1-nk_total/16.d0)/
     3    dble(nk_total/8.d0-1-nk_total/16.d0)*pi/2.d0))**2
      END DO
      
      END SUBROUTINE
