      function dnorm2(x,n)
      
      implicit       none
      integer        i,n 
      complex*16     x(n),y(n),xx  
      real*8         dnorm2  
      
      xx = (0.d0,0.d0)
      do i = 1,n
         xx = xx + conjg(x(i))*x(i)
      enddo
      
      dnorm2 = sqrt(realpart(xx))
      
      return
      end
      
