	 subroutine fftprocess(omegacnt,cnt,kcnt,greens
     1   ,rcvlevel,srclevel)
         
         implicit none

         INCLUDE 'bounds.h'
         INCLUDE 'num.h'
         INCLUDE 'param.h'
      
         real*8    fr(ndim),fi(ndim),rcvlevel,srclevel
         
         complex   greens(ndim,ndim,dim,dim)
         
         integer   wavecntx,wavecnt,cnt,kcnt,omegacnt
          
         do  wavecntx=1,ndim	   
	  do wavecnt=1,ndim 
             fr(wavecnt)=real(greens(wavecntx,wavecnt,cnt,kcnt))
	     fi(wavecnt)=imag(greens(wavecntx,wavecnt,cnt,kcnt))
          end do
    	   
	    call fftpot(2,mlog,ndim,fr,fi)
	 	   
          do wavecnt = 1,ndim
	     greens(wavecntx,wavecnt,cnt,kcnt) = 
     .       cmplx(fr(wavecnt),fi(wavecnt))
	  end do	  	 
	 end do
        
	
	 do  wavecnt = 1,max_rcv	   
	  do wavecntx = 1,ndim 	
             fr(wavecntx) = real(greens(wavecntx,wavecnt,cnt,kcnt))
	     fi(wavecntx) = imag(greens(wavecntx,wavecnt,cnt,kcnt))
	  end do
    	   
	   call fftpot(2,mlog,ndim,fr,fi)
	 
	  do wavecntx = 1,max_rcv
	 
	   greens(wavecntx,wavecnt,cnt,kcnt) = 1d0/float(ndim)*           
     .      cmplx(fr(wavecntx),fi(wavecntx))
	    
	  end do	  	 
	 end do
         return
         end subroutine
