     
      SUBROUTINE mkTdRuRdTu(pabs,Td,Ru,Rd,Tu,layer)

      implicit none
      INCLUDE 'bounds.h'
      INCLUDE 'param.h'
      INCLUDE 'num.h'
      
      integer	layer
      real*8   	pabs

      complex*16 Rd(3,3),Ru(3,3),Td(3,3),Tu(3,3)
      complex*16 Rd_sh,Ru_sh,Td_sh,Tu_sh
      complex*16 D1(4,4),D2(4,4),invar1(2,2),invar2(2,2),invar3(2,2),
     1		dum1(2,2),mu1(2,2),md1(2,2),nu1(2,2),nd1(2,2),
     2        	mu2(2,2),md2(2,2),nu2(2,2),nd2(2,2)

c
c  SH-part
c     
      CALL mkmatSH(Td_sh,Ru_sh,Rd_sh,Tu_sh,layer)
c
c  P-SV
c
      IF ((alpha(layer).eq.alpha(layer+1)).and.
     .	    	  (beta(layer).eq.beta(layer+1)) ) THEN
       
      Td(1,1) = dcmplx(1.d0,0.d0)
      Td(1,2) = dcmplx(0.d0,0.d0)
      Td(1,3) = dcmplx(0.d0,0.d0)
      Td(2,1) = dcmplx(0.d0,0.d0)
      Td(2,2) = dcmplx(1.d0,0.d0)
      Td(2,3) = dcmplx(0.d0,0.d0)
      Td(3,1) = dcmplx(0.d0,0.d0)
      Td(3,2) = dcmplx(0.d0,0.d0)
      Td(3,3) = dcmplx(1.d0,0.d0)

      Rd(1,1) = dcmplx(0.d0,0.d0)
      Rd(1,2) = dcmplx(0.d0,0.d0)
      Rd(1,3) = dcmplx(0.d0,0.d0)
      Rd(2,1) = dcmplx(0.d0,0.d0)
      Rd(2,2) = dcmplx(0.d0,0.d0)
      Rd(2,3) = dcmplx(0.d0,0.d0)
      Rd(3,1) = dcmplx(0.d0,0.d0)
      Rd(3,2) = dcmplx(0.d0,0.d0)
      Rd(3,3) = dcmplx(0.d0,0.d0)

      Rd(1,1) = dcmplx(0.d0,0.d0)
      Rd(1,2) = dcmplx(0.d0,0.d0)
      Rd(1,3) = dcmplx(0.d0,0.d0)
      Rd(2,1) = dcmplx(0.d0,0.d0)
      Rd(2,2) = dcmplx(0.d0,0.d0)
      Rd(2,3) = dcmplx(0.d0,0.d0)
      Rd(3,1) = dcmplx(0.d0,0.d0)
      Rd(3,2) = dcmplx(0.d0,0.d0)
      Rd(3,3) = dcmplx(0.d0,0.d0)

      ELSE
      
      CALL mkD(pabs,D1,layer)
      CALL split(D1,mu1,nu1,md1,nd1)
      CALL mkD(pabs,D2,layer+1)

      CALL split(D2,mu2,nu2,md2,nd2)
      CALL minvar(md1,nd1,md2,nd2,invar1)
      CALL minvar(mu1,nu1,md2,nd2,dum1)
      CALL invert(dum1,invar2,2)
      CALL minvar(mu1,nu1,mu2,nu2,invar3)

      CALL MxM(invar1,invar2,dum1,2)

      Rd(1,1) = -dum1(1,1)
      Rd(1,2) = -dum1(1,2)
      Rd(1,3) = dcmplx(0.d0,0.d0)
      Rd(2,1) = -dum1(2,1)
      Rd(2,2) = -dum1(2,2)
      Rd(2,3) = dcmplx(0.d0,0.d0)
      Rd(3,1) = dcmplx(0.d0,0.d0)
      Rd(3,2) = dcmplx(0.d0,0.d0)
      Rd(3,3) = Rd_sh

      Td(1,1) = ci*invar2(1,1)
      Td(1,2) = ci*invar2(1,2)
      Td(1,3) = dcmplx(0.d0,0.d0)
      Td(2,1) = ci*invar2(2,1)
      Td(2,2) = ci*invar2(2,2)
      Td(2,3) = dcmplx(0.d0,0.d0)
      Td(3,1) = dcmplx(0.d0,0.d0)
      Td(3,2) = dcmplx(0.d0,0.d0)
      Td(3,3) = Td_sh

      CALL MxM(invar2,invar3,dum1,2)

      Ru(1,1) = -dum1(1,1)
      Ru(1,2) = -dum1(1,2)
      Ru(1,3) = dcmplx(0.d0,0.d0)
      Ru(2,1) = -dum1(2,1)
      Ru(2,2) = -dum1(2,2)
      Ru(2,3) = dcmplx(0.d0,0.d0)
      Ru(3,1) = dcmplx(0.d0,0.d0)
      Ru(3,2) = dcmplx(0.d0,0.d0)
      Ru(3,3) = Ru_sh
      
      ENDIF
      
      print *, Td


      CALL transp(Td,Tu,3)
      
      
      print *, Tu
            
      END SUBROUTINE
