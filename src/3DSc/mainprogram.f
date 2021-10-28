      PROGRAM mainprogram

c Author      : Christina Dwiriyanti
c Last changes: December 2002
c Language    : Fortran 77 
c Description : Main progroam to calculate the Greens tensor in the time domain
c               Reference: Brian Kennet's "Seismic wave propagation in stratfied media"
c 
      IMPLICIT NONE
      
      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
                  
      REAL*8   hr(ns),hi(ns)
                  
      COMPLEX  greens_p(dim,dim),greens_h(dim,dim)
            
      COMPLEX  gp(max_rcv, max_rcv, ns,dim)
      
      INTEGER  ii,jj,freqcnt,cnt,nf
      
      write(*,*)'Reading parameters'
      CALL readpar
      nf=fmax
      write(*,*)'Done reading parameters'
      CALL initial
       
      write(*,*)'Start to calculate the 3d green function in the freq'           
        
      DO freqcnt = 1,ns		! loop for omeg	
        IF ( freqcnt .lt. (nf+1)) THEN
          CALL greenstensor(freqcnt,greens_p,greens_h,rcvpos(3),
     .                      srcpos(3))
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
      ! Back to time domain; apply fft           
      DO cnt=1,dim       
        DO ii=1,max_rcv
	        DO jj=1,max_rcv
             DO freqcnt=1,ns	  
                hr(freqcnt) = real(gp(ii,jj,freqcnt,cnt))
	              hi(freqcnt) = imag(gp(ii,jj,freqcnt,cnt))	  
	           END DO
	           CALL fftpot(1,slog,ns,hr,hi)
	           DO freqcnt=1,ns	           
	            gp(ii,jj,freqcnt,cnt) = cmplx(hr(freqcnt),hi(freqcnt))
	           END DO
	        END DO
        END DO
        if(debug .eq. 1) then
          write(*,*)'Index freq-',freqcnt, gp(64,64,128,3)
        endif

      END DO

      write(*,*)'Finish..'

      write(*,*)'Ok, now writing the output of the green function'
      open(unit=70,file='inc1.out',
     .access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)	    	    	    	      

      open(unit=75,file='inc2.out',
     .access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv) 	 	      

      write(70,rec=1)(((real(gp(ii,jj,freqcnt,1)),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
     
      write(75,rec=1)(((real(gp(ii,jj,freqcnt,2)),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
     
      open(unit=80,file='inc3.out',
     . access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)

       write(80,rec=1)(((real(gp(ii,jj,freqcnt,3)),freqcnt=1,ns),
     & jj=1,max_rcv),ii=1,max_rcv)

      CLOSE(70)            
      CLOSE(75)
      CLOSE(80)

340   format(3(e12.6,2x))
800   format(128(e12.6,2x))
  
      RETURN
      END PROGRAM     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
