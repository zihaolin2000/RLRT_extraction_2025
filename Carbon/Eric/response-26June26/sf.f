      
      SUBROUTINE SF(w2,q2,F1p,FLp,F2p,F1n,FLn,F2n)
CCCC   Converts reduced cross sections to structure functions for protons and neutrons  CCCCC 

      IMPLICIT none

      real*8 w2,q2,x,sigtp,siglp,sigtn,sigln,f1p,f2p,fLp
      real*8 f1n,f2n,fLn,pi,pi2,alpha,mp,mp2

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599 
      x = q2/(q2+w2-mp2)

   
      call rescsp(w2,q2,sigTp,sigLp)
      call rescsn(w2,q2,sigTn,sigLn)

      f1p = sigTp/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f1n = sigTn/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLp = sigLp*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLn = sigLn*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f2p = (2.*x*f1p+fLp)/(1.+4.*mp2*x*x/q2)
      f2n = (2.*x*f1n+fLn)/(1.+4.*mp2*x*x/q2)

      return
      end

      
CCC   -----------------
