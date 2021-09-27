program main 
implicit none

!finalval.f\
!resultval.f\
!         solveqn.f\
!         BHmxv.f\
!         Gmxm.f\
!         BGmxm.f\
!         Hpmxv.f\
!         Gfree.f\
!         Gfree_polarz.f\
!         Greentens.f\
!         MxMGreen.f\
!         matgreen.f\
!         kxkyvalue.f\
!         mksource.f \
!         pmksource.f \
!         mkpdom.f \
!         mkmatD.f \
!         mkmatDh.f\
!         Besselfunc.f\
!         mkmatE.f \ 
!         readpar.f \
!         initial.f\
!         zgemm.f\
!         zgeco.f\
!	   zgesl.f\
!         zgemv.f \
!         xerbla.f \
!         lsame.f \
!         matproloog.f \
!         calcfield.f \
!         Gcalcfield.f \
!         initmsg.c \
!         files.c \
!         getpar.c
!# - Include files
!INCLUDE = bounds.h num.h param.h path.h record.h
       IMPLICIT NONE

       INCLUDE 'bounds.h'        
       INCLUDE 'num.h'        
       INCLUDE 'param.h'
	
       INTEGER  freqcnt,cnt,kcnt,ii,jj,ij,kndat,nf        
       REAL*8   he(ns),hi(ns),omega               
       COMPLEX  u_sc(max_rcv,max_rcv,ns,dim),u_k(maxksamp,dim,ns),
     .           solv(maxksamp),  
     .           greens(max_rcv,max_rcv,maxksamp,dim,dim)
            
        
 	 CALL readpar        

        CALL initial        
        
        nf = 47

        DO freqcnt = 1,ns
        
        IF (freqcnt .lt. (int(fmax)+1)) THEN

        WRITE(*,*) 'freqk', freqcnt				
        
        omega = 2d0*pi*float(freqcnt)        

        CALL solveeqnval(freqcnt,solv)
          
        DO kndat=1,ndat        

         DO kcnt = (kndat-1)*dim+1,kndat*dim 
            u_k(kndat,kcnt-(kndat-1)*dim,freqcnt) = solv(kcnt)                      
         END DO
        END DO
  
        ELSE

        DO kndat=1,ndat
         DO kcnt = (kndat-1)*dim+1,kndat*dim                            
            u_k(kndat,kcnt-(kndat-1)*dim,freqcnt) = cmplx(0.0,0.0)
         END DO
        END DO
                
        END IF              
                
        END DO !finish looping freqcnt
  
                                   
        DO freqcnt=1,ns
         
          IF (freqcnt .lt. (int(fmax)+1)) THEN          
               omega = 2d0*pi*float(freqcnt)
         
               CALL resultval(freqcnt,greens)  
                                   
               DO ii=1,max_rcv
                   DO jj=1,max_rcv
                     DO cnt=1,dim
                        u_sc(ii,jj,freqcnt,cnt)=cmplx(0.0,0.0)
                        DO ij=1,ndat         
                           DO kcnt=1,dim                 
                               u_sc(ii,jj,freqcnt,cnt) = u_sc(ii,jj,freqcnt,cnt) +
     .                         greens(ii,jj,ij,cnt,kcnt)*u_k(ij,kcnt,freqcnt)*
     .                         drho(ij)*dx*dy*dz*omega**2 
                           END DO
                        END DO       
                     END DO  
                  END DO
               END DO           

          ELSE

            DO cnt=1,dim
               DO ii=1,max_rcv
                  DO jj=1,max_rcv
                     u_sc(ii,jj,freqcnt,cnt)=cmplx(0.0,0.0)
                  END DO        
                END DO                
            END DO
                   
          END IF                
       
       END DO ! end freqcnt  

cc  Fast Fourier Transform in t-domain

       DO cnt = 1,dim  	
	  DO ii = 1,max_rcv
	    DO jj = 1,max_rcv	 
          
             DO freqcnt = 1,ns          
                he(freqcnt) = real(u_sc(ii,jj,freqcnt,cnt))
                hi(freqcnt) = imag(u_sc(ii,jj,freqcnt,cnt))
             END DO
       	
    	      CALL fftpot(1,slog,ns,he,hi)
 	
             DO freqcnt = 1,ns                    
               u_sc(ii,jj,freqcnt,cnt) = cmplx(he(freqcnt),hi(freqcnt))
            END DO
          
	   END DO	
      	  END DO       
       END DO	  
	                   
      OPEN(unit=70,file='uscatter1.out'
     .,access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)
       
      OPEN(unit=75,file='uscatter2.out'
     .,access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)

      OPEN(unit=80,file='uscatter3.out'
     .,access='direct',form='unformatted',recl=2*4*ns*max_rcv*max_rcv)

      WRITE(70,rec=1)((((REAL(u_sc(ii,jj,freqcnt,1))),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
    
      WRITE(75,rec=1)((((REAL(u_sc(ii,jj,freqcnt,2))),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
     
      WRITE(80,rec=1)((((REAL(u_sc(ii,jj,freqcnt,3))),freqcnt=1,ns),
     .                   ii=1,max_rcv),jj=1,max_rcv)
           
      CLOSE(70)
      CLOSE(75)
      CLOSE(80)
700   STOP
end
       
