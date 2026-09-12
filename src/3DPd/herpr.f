      FUNCTION HERPR(X,Y,N)
      
      IMPLICIT       NONE
      INTEGER        I,N 
      COMPLEX*16     X(N),Y(N) 
      COMPLEX*16     XY,HERPR 
      
      XY = (0.D0,0.D0)
      do 200 I = 1,N
        XY = XY + CONJG(X(I))*Y(I)
  200 continue
      
      HERPR = XY
      
      RETURN
      END
      
