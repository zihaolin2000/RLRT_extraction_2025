      SUBROUTINE NUC12SF(Z,A,w2,q2,state,F1,F2)
CCCCC  Calculates the nuclear elastic and excited state cross sections for 12C    CCCCC
CCCCC  This version smears in nu with a gaussian of width defined below           CCCCC 
CCCCC  Cross section is dsig/domega/dE in ub                                      CCCCC
      
      IMPLICIT NONE
      REAL*8 q2,q2f,e,ep,th,signe,ff,fft2,cs,W1,W2,tau,mott,nuccs,q2t
      REAL*8 sin2,cos2,tan2,recoil,width,norm,nu,nuex,nuel,fs,omega
      INTEGER i,j,IZ,IA,state,F1,F2
      REAL*8 Z,A,mp/0.93827/,hbarc/1973.28/,alpha/7.29735E-03/
      REAL*8 pi/3.14159/,radcon/0.0174533/,mA
      REAL*8 exc(10)/0.0,0.00444,0.00765,0.00964,0.0161,0.0185,0.020,     
     &            0.0230, 0.031, 0.020/                       !!! excited state energies in GeV
      REAL*8 x,tot,epnuc
      
c     width = 0.00225            !!! in GeV
      width = 0.0030
      if(state.EQ.8) width = 0.006
      if(state.EQ.9) width = 0.0075
      if(state.EQ.10) width = 0.007
      norm = width*sqrt(pi) 

      nu  = (w2-mp*mp+q2)/2.0/mp
      
      nuel = 0.0

      F1 = 0.0
      F2 = 0.0
      
      nuex =  nuel+exc(state)
c      omega = max(0.0,nu-nuex)
c     Q2 = 4.*e*(e-nuex)*sin2  !!!! Q^2 at the E' for the given state
      
      q2f = q2/0.1975/0.1975

      if(state.EQ.18) then   !!! multiple excitations in this region so energy appears to shift  !!!
        exc(18) = min(0.01930+0.00021*sqrt(q2f),0.01955)
      endif
     
      
      fs = exp(-1.0*(nu-nuex)**2./width/width)  !!! Gaussian probability based on distance from peak
      fs = fs/norm              !!! normalized
      if(((nu-nuex)/width).GT.4.0) fs = 0.0
 
c      write(6,*) omega, omega+width,log(omega), (log(e-nuel)-log(width)) 

      FF = 0.0
      if(state.LE.9) call NUCFFS12C(A,Z,Q2,state,FF) !!! longitudinal FFs
        

CCCC  Now calculate the total transverse FF^2

      FFT2 = 0.0
      if(state.EQ.10)
     &     FFT2 = 2.0E-4+0.054*q2-0.5017*q2*q2+1.7807*q2*q2*q2
     &           -3.0511*q2**4.0+ 2.535*q2**5.0-0.8198*q2**6.0

      
CCCC
      
      if(q2.GE.0.3) FF = 0.0

c      if(state.EQ.8) FF = FF/sqrt(3.0)
c      if(state.EQ.7) FF = FF/sqrt(2.0)
      
c      write(6,*) "Nuccs12cs: ",state,e,ep,th,q2,FF
      
      W2  = (Z*FF)**2+Z*Z*FFT2
      W1 = 0.0
      W1 = Z*Z*FFT2
                     

      F2 = nu*W2
      F1 = mp*W1  
      
c      nuccs = fs

      if(state.EQ.9) then
         F2 = 0.2*F2      
      elseif(state.EQ.8) then
         F2 = 0.75*F2
      elseif(state.EQ.7) then
         F2 = 1.0*F2
      elseif(state.EQ.6) then
         F2 = 1.0*F2
      endif
      F2 = fs*F2
      F1 = fs*F1

 
 1001 format(4F9.3,3E11.3)
      
      return

      end
