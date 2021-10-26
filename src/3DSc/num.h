      integer   nsamp,logns,nroscils
      real*8      domega,dt,rcvpos(dim),srcpos(dim),
     1	     velo_src,oscilfreq(maxoscils),freqstart,freqend,pi,fmax,fmin
      
      complex*16 ci

      common/timesampl/logns,nsamp,dt,domega
      common/movingsrc/rcvpos,srcpos,velo_src,oscilfreq,nroscils
      common/wavelet/freqend,freqstart,fmin,fmax
      common/hulpvars/ci,pi

