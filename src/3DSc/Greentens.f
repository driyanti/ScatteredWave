      SUBROUTINE greenstensor(omegacnt,greens_p,greens_h,rcvlevel,
     . srclevel)
      
c
c author:   	  Christina Dwiriyanti
c last changes:   December 2000
c
c subroutine to calculate Greens function in slowness domain
c p_x is vast (p_x = p1)
c p_y is the integration parameter
c 
c September 2000: implementation force direction (mksource)
c 
      implicit none

      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
      
      integer	  kcnt,wavecnt,layer,wavecntx,ii,omegacnt
      
      real*8 	  p_y(2*maxksamp+1),ptaper(2*maxksamp+1),omega
                             
      real*8      p1,p2,pabs,psum,k_x(ndim),k_y(ndim),epsilon
      
      real*8      length_k,delta_k,rcvlevel,srclevel,p_x(2*maxksamp+1)
      
      
      complex*16  matD(4,4),matDh(2,2),Su(3,3),Sd(3,3),psrc,wlin(ns),
     .		  v_src_min(6),v_src_plus(6),v_rcv(6),b_rcv(6),hlpvar,
     .            green_halfspace(dim*dim)
     
      complex  greens_p(dim,dim),greens_h(dim,dim)

      common/pvalues/p_y,p_x,ptaper
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c v_src_min 	  	wavefield above src level
c
c v_src_plus	  	wavefield below src level
c
c matDcombi 	  	combination of matrices D
c     	    	  	(P-SV) and (SH)
c 
c SU,SD     	  	scattering matrices (Kennett)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c do NOT calculate for omega==0.0
c
c      nk_totpos = 0
      
      omega = 2d0*pi*float(omegacnt)
       
      IF (omega.eq.0.d0) GOTO 750
      
      length_k  = 2d0*pi
      delta_k = length_k/float(ndim)
       

c length_k_1=length_k_2 denotes the length interval in the wavenumber domain c(k_x,k_y)
      
      do ii=1,ndim/2+1
	k_x(ii) = float(ii-1)*delta_k
	k_y(ii) = k_x(ii)
      end do
       
      do ii=ndim/2+2,ndim
	k_x(ii) = float(ii-1)*delta_k-length_k
	k_y(ii) = k_x(ii)	
      end do
      
      CALL waveform(wlin)
                   
		               
          epsilon = delta_k/1000d0
          p1 = epsilon
          p2 = epsilon
        
	  p1=0d0
	  p2=0d0
	  
	  psum = p1**2 + p2**2
	  pabs = sqrt(psum)
          
	  psrc= cmplx(1.0,0.0)

c          
c
c initialize q_alpha and q_beta with appropriate signs
c sign of q_alpha and q_beta should correspond to the fourier transform
c see De Hoop (waves in elastic media, dissipative medium and our fourier trnsf)
c
	
        DO layer = 1,nlay+1
          hlpvar        = sqrt(-psum + 1.d0/alpha(layer)**2)
	  qalpha(layer) = cmplx(abs(dble(hlpvar)),dble(-ci*hlpvar)) 
          hlpvar        = sqrt(-psum + 1.d0/beta(layer)**2)
	  qbeta(layer)  = cmplx(abs(dble(hlpvar)),dble(-ci*hlpvar)) 
        END DO	  ! end layer loop
        

c        
c construct matDcombi matrix to be used later
c
 	CALL mkD(pabs,matD,rcvlayer)          
          	
        CALL mkDh(matDh,rcvlayer)

c
c matrices SU and SD are constructed
c matrices RD, RU, TD, TU are constructed as well
c for each layer
c
        CALL  mkgreen_halfspace(green_halfspace,omega,p1,p2
     .     ,rcvlevel,srclevel)

          greens_h(1,1) = wlin(omegacnt)*psrc*green_halfspace(1)
          greens_h(1,2) = wlin(omegacnt)*psrc*green_halfspace(2)
          greens_h(1,3) = wlin(omegacnt)*psrc*green_halfspace(3)
          greens_h(2,1) = wlin(omegacnt)*psrc*green_halfspace(4)
          greens_h(2,2) = wlin(omegacnt)*psrc*green_halfspace(5)
          greens_h(2,3) = wlin(omegacnt)*psrc*green_halfspace(6)
          greens_h(3,1) = wlin(omegacnt)*psrc*green_halfspace(7)
          greens_h(3,2) = wlin(omegacnt)*psrc*green_halfspace(8)
          greens_h(3,3) = wlin(omegacnt)*psrc*green_halfspace(9)
                  
	CALL matrix_proloog(pabs,omega,Su,Sd,srclevel)
	        
	DO kcnt = 1,dim	! for each force direction
	    
            
             epsilon = delta_k/1000d0
             p1 = epsilon
             p2 = epsilon
             
	     p1=0d0
	     p2=0d0

             psrc = exp(-ci*omega*(p1*srcpos(1)+p2*srcpos(2)))

             psrc= cmplx(1.0,0.0)
	     
c at source level the field in terms of up- and downgoing waves 
c is constructed.
c output: v_src_min and v_src_plus represent field above and below 
c the source level
c kcnt : direction force
c
	  
	  CALL mksource(p1,p2,omega,Su,Sd,v_src_min,v_src_plus,kcnt)
                        
          CALL calcfield(pabs,omega,v_src_min,v_src_plus,v_rcv,
     .         rcvlevel,srclevel)
     
c
c DISPLACEMENTS
c

     	    b_rcv(1) = matD(1,3)*(v_rcv(4) - v_rcv(1)) +
     .                   matD(1,4)*(v_rcv(2) + v_rcv(5))

	    b_rcv(2) = matD(2,1)*(v_rcv(4) + v_rcv(1)) +
     .                   matD(2,4)*(v_rcv(5) - v_rcv(2))

	    b_rcv(3) = matDh(1,1)*(v_rcv(3) + v_rcv(6))

c
c STRESSES
c
	
!	  b_rcv(4) = matD(3,1)*(v_rcv(4) + v_rcv(1)) +
!     .	  	     matD(3,4)*(v_rcv(5) - v_rcv(2))

!	  b_rcv(5) = matD(4,3)*(v_rcv(4) - v_rcv(1)) +
!     .	  	     matD(4,4)*(v_rcv(5) + v_rcv(2))

!	  b_rcv(6) = matDh(2,2)*(-v_rcv(3) + v_rcv(6))

c
c Transformation (U,V,W) to (u,v,w) (slowness domain) 
c

          IF ((p2.eq.0.d0).and.(p1.eq.0.d0)) THEN
            
     
	      greens_p(1,kcnt) = wlin(omegacnt)*psrc*
     .              (ci*(b_rcv(2)-b_rcv(3))/sqrt(2d0))
     
              greens_p(2,kcnt) = wlin(omegacnt)*psrc*
     .              (ci*(b_rcv(2)+b_rcv(3))/sqrt(2d0))
     
              greens_p(3,kcnt) = wlin(omegacnt)*(b_rcv(1))*psrc

	    ELSE !NOT ((p2.eq.0.d0).and.(p1.eq.0.d0))

	      greens_p(1,kcnt) = wlin(omegacnt)*psrc*
     .              (ci*(p1/pabs)*b_rcv(2)-ci*(p2/pabs)*b_rcv(3))
     
              greens_p(2,kcnt) = wlin(omegacnt)*psrc*
     .              (ci*(p2/pabs)*b_rcv(2)+ci*(p1/pabs)*b_rcv(3))
     
              greens_p(3,kcnt) = wlin(omegacnt)*(b_rcv(1))*psrc
     
	END IF
                      
	END DO    ! end kcnt loop
      
	      
750   END SUBROUTINE

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine waveform(wlin)
       
       include 'bounds.h'
       
       integer ii,jj,fqmax
       
       double precision w(ns,ns),wr(ns),wi(ns),df(ns)
       
       complex*16 wlin(ns),ci
       
       ci=cmplx(0.0,1.0)
              
       open(unit=9,file='wlin256.out',status='old')
         
       read(9,*)((w(ii,jj),jj=1,2),ii=1,ns)       
       
       fqmax = 47

       do ii=1,ns
         df(ii) = w(ii,1)
	 wr(ii) = w(ii,2)
         wi(ii) = 0.0           
       end do
       
       call fftpot(2,slog,ns,wr,wi)
       
       do ii=1,fqmax
        wlin(ii)=cmplx(wr(ii),wi(ii))
       end do
       
       do ii=(fqmax+1),ns
        wlin(ii)=ci*0.0
       end do	
              
       close(9)   
           
       return
       end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    	subroutine fftpot (kdft,mlog,ndim,rre,rim)
c
cccc    implicit double precision (a-h,o-z)
        double precision rre(ndim),rim(ndim)
        double precision cosh,freq,pi,rimh,rreh,sinh,wi,wih,wr,wrh
     0  integer index,k,kdft,mh,mhlk,ml,mlk,mlog,
     1          nd,ndh,ndim,ndimm1,nh,nhdim
c
c
        pi=4d0*datan(1d0)
	nhdim=ndim/2
	ndimm1=ndim-1
	index=1
	
	do 10 nd=1,ndimm1
	  if (nd.ge.index) goto 4
	  rreh=rre(index)
	  rimh=rim(index)
	  rre(index)=rre(nd)
	  rim(index)=rim(nd)
	  rre(nd)=rreh
	  rim(nd)=rimh
4	  k=nhdim
6	  if (k.ge.index) goto 8
	  index=index-k
	  k=k/2
	  goto 6
8	  index=index+k
10	continue
	
	do 40 ml=1,mlog
	  mlk=2**ml
	  mhlk=mlk/2
	  wr=1.
	  wi=0.
	  freq=pi/mhlk
	  cosh=dcos(freq)
	  sinh=dsin(freq)
	  do 30 mh=1,mhlk
	     do 20 nd=mh,ndim,mlk
		ndh=nd+mhlk
		rreh=wr*rre(ndh)-wi*rim(ndh)
		rimh=wi*rre(ndh)+wr*rim(ndh)
		rre(ndh)=rre(nd)-rreh
		rim(ndh)=rim(nd)-rimh
		rre(nd)=rre(nd)+rreh
		rim(nd)=rim(nd)+rimh
20	  continue
	  wrh=wr
	  wih=wi
	  wr=wrh*cosh-wih*sinh
	  wi=wrh*sinh+wih*cosh
30	  continue
40	  continue

	  do 60 nd=1,ndim
	  rre(nd)=rre(nd)/dsqrt(dfloat(ndim))
	  rim(nd)=rim(nd)/dsqrt(dfloat(ndim))
60	  continue
	  if (kdft.eq.2) goto 64
          do 50 nh=2,nhdim
	  index=ndim-nh+2
	  rreh=rre(index)
	  rimh=rim(index)
	  rre(index)=rre(nh)
          rim(index)=rim(nh)
          rre(nh)=rreh
          rim(nh)=rimh
50        continue
64        return
          end
	  
ccccccccccccccccccccccccccccccccccccccccccc
