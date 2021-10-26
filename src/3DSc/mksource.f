      SUBROUTINE mksource(p1,p2,omega,Su,Sd,v_src_min,v_src_plus,
     .						force_direction)
      implicit none

!
! extra: times omega
!
      
      INCLUDE 'bounds.h'
      INCLUDE 'force.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'
      
      integer   force_direction
      real*8	  p1,p2,omega,p,forcefactor
      complex*16  Sd(3,3),Su(3,3),F0(6),
     .		  sigma0(6),sigma_up(3),sigma_down(3),
     .		  matD(4,4),Dinv(4,4),Dh(2,2),Dinvh(2,2),Dcombi(6,6),
     .		  v_src_up_min(3),v_src_up_plus(3),
     .		  v_src_down_min(3),v_src_down_plus(3),
     .		  v_src_min(6),v_src_plus(6),
     .		  c0,EYE(3,3),hlpvec1(3),hlpvec2(3),
     .		  hlpmat1(3,3),hlpmat2(3,3)
     
      EYE(1,1) = cmplx(1.d0,0.d0)
      EYE(1,2) = cmplx(0.d0,0.d0)
      EYE(1,3) = cmplx(0.d0,0.d0)
      EYE(2,1) = cmplx(0.d0,0.d0)
      EYE(2,2) = cmplx(1.d0,0.d0)
      EYE(2,3) = cmplx(0.d0,0.d0)
      EYE(3,1) = cmplx(0.d0,0.d0)
      EYE(3,2) = cmplx(0.d0,0.d0)
      EYE(3,3) = cmplx(1.d0,0.d0)

      c0          = cmplx(0.d0,0.d0)
      p	          = sqrt(p1**2 + p2**2)
      forcefactor = forcexp

c	
c F0 = (F0(U),F0(V),F0(W),F0(P)=F3,F0(S)=FV,F0(T)=FH)^T
c	
	F0(1) = c0
	F0(2) = c0
	F0(3) = c0
	F0(4) = c0
	F0(5) = c0
	F0(6) = c0
	IF (force_direction.eq.1) THEN
c  
c  F_V = j(p_1*f_1 + p_2*f_2)/p =  j p_1 / p * f_1
c  F_H = j(p_1*f_2 - p_2*f_1)/p = -j p_2 / p * f_1
c
          IF (p.ne.0.0) THEN
	      F0(5) =  ci*p1/(omega*p)
	      F0(6) = -ci*p2/(omega*p)
          ELSE
	      F0(5) =  ci/(sqrt(2.d0)*omega)
	      F0(6) =  ci/(sqrt(2.d0)*omega)
	  ENDIF

c          IF (p.ne.0.0) THEN
c	      F0(5) =  ci*p1/(omega**2*p)
c	      F0(6) = -ci*p2/(omega**2*p)
c          ELSE
c	      F0(5) =  ci/(sqrt(2.d0)*omega**2)
c	      F0(6) = -ci/(sqrt(2.d0)*omega**2)
c	  ENDIF

!
! extra omega terms have been added to 
!	  
	ELSE
	  IF (force_direction.eq.2) THEN
c
c  F_V = j(p_1*f_1 + p_2*f_2)/p =  j p_2 / p * f_2
c  F_H = j(p_1*f_2 - p_2*f_1)/p =  j p_1 / p * f_2
c
      	    IF (p.ne.0.0) THEN
     	      F0(5) =  ci*p2/(omega*p)
    	      F0(6) =  ci*p1/(omega*p)
	    ELSE
     	      F0(5) =  ci/(sqrt(2.d0)*omega)
    	      F0(6) =  ci/(sqrt(2.d0)*omega)
	    END IF
c            IF (p.ne.0.0) THEN
c     	      F0(5) =  ci*p2/(omega**2*p)
c    	      F0(6) =  ci*p1/(omega**2*p)
c	    ELSE
c     	      F0(5) =  ci/(sqrt(2.d0)*omega**2)
c    	      F0(6) =  ci/(sqrt(2.d0)*omega**2)
c	    END IF
	  ELSE
     	    F0(4) =  -1.d0/omega 
!
! extra omega term in F0(4) unknown
! this minus sign has been added to fit Hankel function
! don't ask me why!!!!
!
	  ENDIF
	ENDIF
	
	F0(4) = forcefactor*F0(4)
	F0(5) = forcefactor*F0(5)
	F0(6) = forcefactor*F0(6)

	CALL mkD(p,matD,srclayer)
	CALL mkDinv(matD,Dinv)
	CALL mkDh(Dh,srclayer)
	CALL mkDinvh(Dh,Dinvh)
	
	CALL build6x6(Dinv,Dinvh,Dcombi)
	CALL matxvec(Dcombi,F0,sigma0,6)

	sigma_up(1)   = sigma0(1) 
	sigma_up(2)   = sigma0(2)
	sigma_up(3)   = sigma0(3)

	sigma_down(1) = sigma0(4) 
	sigma_down(2) = sigma0(5)
	sigma_down(3) = sigma0(6) 
	
	CALL MxM(Su,Sd,hlpmat1,3)
	CALL matsub(EYE,hlpmat1,hlpmat2)

	CALL invert(hlpmat2,hlpmat1,3)
	CALL matxvec(Su,sigma_up,hlpvec1,3)
	CALL vecsub(sigma_down,hlpvec1,hlpvec2)
	CALL matxvec(hlpmat1,hlpvec2,v_src_down_plus,3)
	CALL matxvec(Sd,v_src_down_plus,v_src_up_plus,3)

	CALL MxM(Sd,Su,hlpmat1,3)
	CALL matsub(EYE,hlpmat1,hlpmat2)
	CALL invert(hlpmat2,hlpmat1,3)
	CALL matxvec(Sd,sigma_down,hlpvec1,3)
	CALL vecsub(hlpvec1,sigma_up,hlpvec2)
	CALL matxvec(hlpmat1,hlpvec2,v_src_up_min,3)	
	CALL matxvec(Su,v_src_up_min,v_src_down_min,3)
			
	v_src_min(1)  = v_src_up_min(1)
	v_src_min(2)  = v_src_up_min(2)
	v_src_min(3)  = v_src_up_min(3)
	v_src_min(4)  = v_src_down_min(1)
	v_src_min(5)  = v_src_down_min(2)
	v_src_min(6)  = v_src_down_min(3)

	v_src_plus(1) = v_src_up_plus(1)
	v_src_plus(2) = v_src_up_plus(2)
	v_src_plus(3) = v_src_up_plus(3)
	v_src_plus(4) = v_src_down_plus(1)
	v_src_plus(5) = v_src_down_plus(2)
	v_src_plus(6) = v_src_down_plus(3)

	RETURN
        
        END SUBROUTINE	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
