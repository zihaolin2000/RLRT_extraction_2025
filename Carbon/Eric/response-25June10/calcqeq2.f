      PROGRAM CALCQEQ2

      IMPLICIT none
 
      real*8 e,theta,q2,qf,mp/0.93827/,mp2,sin2,ep

      mp2 = mp*mp
      
      read(5,*) e,theta

      sin2 = sin(3.14159/180.0*theta/2.0)
      sin2 = sin2*sin2
      ep = mp*e/(2.0*e*sin2+mp)
      q2 = 4.0*e*ep*sin2
      qf = sqrt(q2)/0.1975
      
      write(6,1000) e,theta,sqrt(q2),qf

 1000 format(4f8.3)
      end


 









