      SUBROUTINE calcfield(pabs,omega,v_src_min,v_src_plus,v_rcv,
     .            rcvlevel,srclevel)

c
c author:         Auke Ditzel
c date:           October 2000
c
c this subroutine calculates the field initially given at source level,
c at the receiver level in terms of up- and downgoing waves, by propagating
c the field from source to receiver, making use og the propagator matrix Q
c
c so far we have only implemented the case of srclayer==rcvlayer
c also in different layers now August 2002: Auke Ditzel
c
      implicit none

      INCLUDE 'bounds.h'
      INCLUDE 'param.h'         ! srclayer & rcvlayer defined
      INCLUDE 'num.h'
cs      INCLUDE 'embank.h'
      
      integer    zsgtzr,hlplayer
      real*8     omega,srclevel,rcvlevel,pabs 
      complex*16 v_src_min(2*dim),v_src_plus(2*dim),v_rcv(2*dim),
     1           Eu(dim,dim),Ed(dim,dim),
     2           Td(dim,dim),Rd(dim,dim),Ru(dim,dim),Tu(dim,dim),
     3           vu(dim),vd(dim),dumu(dim),dumd(dim),
     4           hlpvec1(dim),hlpvec2(dim),hlpmat1(dim,dim)
      
      IF (srclevel.ge.rcvlevel) THEN
        zsgtzr = 1
      ELSE
        zsgtzr = 0
      END IF

      vu(1) = zsgtzr*v_src_min(1)+(1-zsgtzr)*v_src_plus(1)
      vu(2) = zsgtzr*v_src_min(2)+(1-zsgtzr)*v_src_plus(2)
      vu(3) = zsgtzr*v_src_min(3)+(1-zsgtzr)*v_src_plus(3)
      vd(1) = zsgtzr*v_src_min(4)+(1-zsgtzr)*v_src_plus(4)
      vd(2) = zsgtzr*v_src_min(5)+(1-zsgtzr)*v_src_plus(5)
      vd(3) = zsgtzr*v_src_min(6)+(1-zsgtzr)*v_src_plus(6)
        
      IF ((srclayer.eq.rcvlayer).and.(rcvlevel.ne.srclevel)) THEN
      
        CALL mkE(omega,rcvlevel-srclevel,Ed,rcvlayer)

        IF (abs(Ed(1,1)).gt.(0.0)) THEN
          Eu(1,1) = 1.d0/Ed(1,1)
          Eu(2,2) = 1.d0/Ed(2,2)
          Eu(3,3) = 1.d0/Ed(3,3)
        ELSE
          Eu(1,1) = cmplx(0.d0,0.d0)
          Eu(2,2) = cmplx(0.d0,0.d0)
          Eu(3,3) = cmplx(0.d0,0.d0)
        ENDIF     
        
        Eu(1,2) = cmplx(0.d0,0.d0)
        Eu(1,3) = cmplx(0.d0,0.d0)
        Eu(2,1) = cmplx(0.d0,0.d0)
        Eu(2,3) = cmplx(0.d0,0.d0)
        Eu(3,1) = cmplx(0.d0,0.d0)
        Eu(3,2) = cmplx(0.d0,0.d0)

        CALL matxvec(Eu,vu,dumu,3)
        CALL matxvec(Ed,vd,dumd,3)

        vu(1) = dumu(1)
        vu(2) = dumu(2)
        vu(3) = dumu(3)
        vd(1) = dumd(1)
        vd(2) = dumd(2)
        vd(3) = dumd(3)

      ELSE IF (srclayer.gt.rcvlayer) THEN
      
        hlplayer = srclayer
        
        DO WHILE (hlplayer.gt.rcvlayer) 
        
          IF (hlplayer.eq.srclayer) THEN          

            CALL mkE(omega,z(srclayer-1)-srclevel,Ed,srclayer)          

c            CALL mkE(omega,srclevel-z(hlplayer-1),Ed,hlplayer)          

          ELSE          
            CALL mkE(omega,-h(hlplayer),Ed,hlplayer)
          ENDIF
         
c         print*,'z', z(hlplayer-1),srclevel,hlplayer-1
          IF (abs(Ed(1,1)).gt.(0.0)) THEN
            Eu(1,1) = 1.d0/Ed(1,1)
            Eu(2,2) = 1.d0/Ed(2,2)
            Eu(3,3) = 1.d0/Ed(3,3)
          ELSE
            Eu(1,1) = cmplx(0.d0,0.d0)
            Eu(2,2) = cmplx(0.d0,0.d0)
            Eu(3,3) = cmplx(0.d0,0.d0)
          ENDIF
        
          Eu(1,2) = cmplx(0.d0,0.d0)
          Eu(1,3) = cmplx(0.d0,0.d0)
          Eu(2,1) = cmplx(0.d0,0.d0)
          Eu(2,3) = cmplx(0.d0,0.d0)
          Eu(3,1) = cmplx(0.d0,0.d0)
          Eu(3,2) = cmplx(0.d0,0.d0)

          CALL matxvec(Eu,vu,dumu,3)
          CALL matxvec(Ed,vd,dumd,3)

          vu(1) = dumu(1)
          vu(2) = dumu(2)
          vu(3) = dumu(3)
          vd(1) = dumd(1)
          vd(2) = dumd(2)
          vd(3) = dumd(3)
!
! now the field just below the (hlplayer-1)-th interface is known
!          
          CALL mkTdRuRdTu(pabs,Td,Ru,Rd,Tu,hlplayer-1)
c          CALL mkTdRuRdTu(pabs,Td,Ru,Rd,Tu,rcvlayer)

!
! compute v_D^-
!          
          CALL matxvec(Ru,vu,hlpvec1,3)
          CALL vecsub(vd,hlpvec1,hlpvec2)
          CALL invert(Td,hlpmat1,3)
          CALL matxvec(hlpmat1,hlpvec2,vd,3)

!
! compute v_U^-
!          

          CALL matxvec(Rd,vd,hlpvec1,3)
          CALL matxvec(Tu,vu,hlpvec2,3)
          CALL vecsum(hlpvec1,hlpvec2,vu)
          
          hlplayer = hlplayer-1 
        
        END DO
!
! go from interface to rcvlevel
!
c        CALL mkE(omega,rcvlevel-z(rcvlayer),Ed,rcvlayer)

        CALL mkE(omega,rcvlevel-z(rcvlayer-1),Ed,rcvlayer)
        
c	print*,'m',rcvlevel,z(rcvlayer-1)       
        
        IF (abs(Ed(1,1)).gt.(0.0)) THEN
          Eu(1,1) = 1.d0/Ed(1,1)
          Eu(2,2) = 1.d0/Ed(2,2)
          Eu(3,3) = 1.d0/Ed(3,3)
        ELSE
          Eu(1,1) = cmplx(0.d0,0.d0)
          Eu(2,2) = cmplx(0.d0,0.d0)
          Eu(3,3) = cmplx(0.d0,0.d0) 
        ENDIF
        
        Eu(1,2) = cmplx(0.d0,0.d0)
        Eu(1,3) = cmplx(0.d0,0.d0)
        Eu(2,1) = cmplx(0.d0,0.d0)
        Eu(2,3) = cmplx(0.d0,0.d0)
        Eu(3,1) = cmplx(0.d0,0.d0)
        Eu(3,2) = cmplx(0.d0,0.d0)

        CALL matxvec(Eu,vu,dumu,3)
        CALL matxvec(Ed,vd,dumd,3)

        vu(1) = dumu(1)
        vu(2) = dumu(2)
        vu(3) = dumu(3)
        vd(1) = dumd(1)
        vd(2) = dumd(2)
        vd(3) = dumd(3)
!        
! done
!              
       ELSE      ! srclayer.lt.rcvlayer
      
        hlplayer = srclayer
                   
        DO WHILE (hlplayer.lt.rcvlayer) 
        
          IF (hlplayer.eq.srclayer) THEN
            CALL mkE(omega,z(hlplayer)-srclevel,Ed,hlplayer)
          ELSE
            CALL mkE(omega,h(hlplayer),Ed,hlplayer)
          ENDIF

c        print*,'cel',z(hlplayer),srclevel
         
          IF (abs(Ed(1,1)).gt.(0.0)) THEN
          Eu(1,1) = 1.d0/Ed(1,1)
          Eu(2,2) = 1.d0/Ed(2,2)
          Eu(3,3) = 1.d0/Ed(3,3)
          ELSE
          Eu(1,1) = cmplx(0.d0,0.d0)
          Eu(2,2) = cmplx(0.d0,0.d0)
          Eu(3,3) = cmplx(0.d0,0.d0)
          ENDIF
        
          Eu(1,2) = cmplx(0.d0,0.d0)
          Eu(1,3) = cmplx(0.d0,0.d0)
          Eu(2,1) = cmplx(0.d0,0.d0)
          Eu(2,3) = cmplx(0.d0,0.d0)
          Eu(3,1) = cmplx(0.d0,0.d0)
          Eu(3,2) = cmplx(0.d0,0.d0)


          CALL matxvec(Eu,vu,dumu,3)
          CALL matxvec(Ed,vd,dumd,3)

          vu(1) = dumu(1)
          vu(2) = dumu(2)
          vu(3) = dumu(3)
          vd(1) = dumd(1)
          vd(2) = dumd(2)
          vd(3) = dumd(3)

!
! now the field just below the (hlplayer)-th interface is known
!         
          CALL mkTdRuRdTu(pabs,Td,Ru,Rd,Tu,hlplayer)          
!
! compute v_U^+
!          
          CALL matxvec(Rd,vd,hlpvec1,3)          
          CALL vecsub(vu,hlpvec1,hlpvec2)
          CALL invert(Tu,hlpmat1,3)
          CALL matxvec(hlpmat1,hlpvec2,vu,3)
!
! compute v_D^+
!          
          CALL matxvec(Ru,vu,hlpvec1,3)
          CALL matxvec(Td,vd,hlpvec2,3)
          CALL vecsum(hlpvec1,hlpvec2,vd)

          hlplayer = hlplayer + 1
        
        END DO
!
! go from interface to rcvlevel
!
c        CALL mkE(omega,rcvlevel-z(rcvlayer-1),Ed,rcvlayer)
        CALL mkE(omega,-rcvlevel+z(rcvlayer-1),Ed,rcvlayer)
        
        IF (abs(Ed(1,1)).gt.(0.0)) THEN
        
          Eu(1,1) = 1.d0/Ed(1,1)
          Eu(2,2) = 1.d0/Ed(2,2)
          Eu(3,3) = 1.d0/Ed(3,3)
        
        ELSE
        
          Eu(1,1) = cmplx(0.d0,0.d0)
          Eu(2,2) = cmplx(0.d0,0.d0)
          Eu(3,3) = cmplx(0.d0,0.d0) 
        
        ENDIF
       
          Eu(1,2) = cmplx(0.d0,0.d0)
          Eu(1,3) = cmplx(0.d0,0.d0)
          Eu(2,1) = cmplx(0.d0,0.d0)
          Eu(2,3) = cmplx(0.d0,0.d0)
          Eu(3,1) = cmplx(0.d0,0.d0)
          Eu(3,2) = cmplx(0.d0,0.d0)

        CALL matxvec(Eu,vu,dumu,3)
        CALL matxvec(Ed,vd,dumd,3)

          vu(1) = dumu(1)
          vu(2) = dumu(2)
          vu(3) = dumu(3)
          vd(1) = dumd(1)
          vd(2) = dumd(2)
          vd(3) = dumd(3)
!
! done
!                       
      END IF
      
      v_rcv(1) = vu(1)
      v_rcv(2) = vu(2)
      v_rcv(3) = vu(3)
      v_rcv(4) = vd(1)
      v_rcv(5) = vd(2)
      v_rcv(6) = vd(3)
      
      END SUBROUTINE
