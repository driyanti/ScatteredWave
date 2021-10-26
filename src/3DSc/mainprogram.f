      PROGRAM GREENTENSOR
      
      IMPLICIT NONE
      
      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
                  
      REAL*8   hr(ns),hi(ns)
                  
      COMPLEX  greens_p(dim,dim),greens_h(dim,dim)
            
      COMPLEX  gp(max_rcv, max_rcv, ns,dim)
      real*8   gn(max_rcv, max_rcv, ns)      
      
      INTEGER  ii,jj,freqcnt,cnt,nf
      
      REAL  gpr(max_rcv,ns), gr(max_rcv,max_rcv,ns) 

      !max number of freq 
      nf = 47
      
      ! Reading parameter
      CALL readpar
      
      ! Initial
      CALL initial

      write(*,*) "Now, computation of the greentensor, G11, G12, G13 will start..."
    
      DO freqcnt = 1,ns		! loop for omeg
        IF ( freqcnt .lt. (nf+1)) THEN
        if (debug .eq. 1) then
           print*,'start',freqcnt
        endif        	   	  	   
        CALL greenstensor(freqcnt,greens_p,greens_h,rcvpos(3),
     .                  srcpos(3))

          DO cnt = 1,dim               ! loop for direction
            DO ii=1,max_rcv
              DO jj=1,max_rcv
	              gp(ii,jj,freqcnt,cnt) = greens_p(cnt,3)
              END DO    ! end loop  p2
            END DO   ! end loop  p1
          END DO
                     	    
        ELSE

          DO cnt=1,dim       
             DO ii=1,max_rcv
               DO jj=1,max_rcv
	                gp(ii,jj,freqcnt,cnt) = cmplx(0.0,0.0)
               END DO 
             END DO
           END DO ! end loop cnt	
      END IF
      
      END DO   ! end loop freqcnt
                 
      DO cnt=1,dim       
         DO ii=1,max_rcv
	     DO jj=1,max_rcv
	 
             DO freqcnt=1,ns	  
                hr(freqcnt) = real(gp(ii,jj,freqcnt,cnt))
	              hi(freqcnt) = imag(gp(ii,jj,freqcnt,cnt))	  
	       END DO
	 
	       CALL fftpot(1,slog,ns,hr,hi)  	    
	 
	      DO freqcnt=1,ns	           
	         gp(ii,jj,freqcnt,cnt) =  cmplx(hr(freqcnt),hi(freqcnt))
           ! Only the vertical direction
               if (cnt .eq. 3) then
                   gn(ii,jj,freqcnt) = real(gp(ii,jj,freqcnt,3))
               endif    
	      END DO
         
	     END DO

         END DO
      END DO

      write(*,*) "write the output now"
      open(unit=70,file='G11.out',
     .access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)	    	    	    	      

      open(unit=75,file='G12.out',
     .access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)	    	    

      open(unit=80,file='G13.out',
     .access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)	    	 	      

      write(70,rec=1)((((real(gp(ii,jj,freqcnt,1))),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
     
      write(75,rec=1)((((real(gp(ii,jj,freqcnt,2))),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
     
      write(80,rec=1)((((real(gp(ii,jj,freqcnt,2))),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
      
     write(*,*) "Done my friend!"
      CLOSE(70)            
      CLOSE(75)
      CLOSE(80)

340   format(3(e12.6,2x))
800   format(128(e12.6,2x))
  
      RETURN
      END PROGRAM     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 