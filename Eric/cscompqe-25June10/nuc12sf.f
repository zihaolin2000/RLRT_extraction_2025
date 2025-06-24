      SUBROUTINE NUC12SF(Z,A,nu,q2p,state,F1,FL)
CCCCC  Calculates the nuclear elastic and excited state cross sections for 12C    CCCCC
CCCCC  This version smears in nu with a gaussian of width defined below           CCCCC 
CCCCC  Cross section is dsig/domega/dE in ub                                      CCCCC
      
      IMPLICIT NONE
      REAL*8 q2,q2p,q2f,e,ep,th,signe,ff,fft,cs,W1,W2,WL,tau,mott,nuccs
      REAL*8 sin2,cos2,tan2,recoil,width,norm,nu,nuex,nuel,fs,omega
      REAL*8 qv2,q2s,q2vs,smwid
      INTEGER i,j,IZ,IA,state,ST
      REAL*8 Z,A,mp/0.93827/,hbarc/1973.28/,alpha/7.29735E-03/
      REAL*8 pi/3.14159/,radcon/0.0174533/,mA
      REAL*8 exc(21)/0.0,0.00444,0.00765,0.00964,0.01084,0.0137,0.0151,    !!! excited state energies in GeV
     &     0.0161,0.0183,0.020,0.0230,0.0315,0.042,0.0151,0.0161,0.0166,
     &     0.0181,0.0193,0.0206,0.0235,0.0315 /
      REAL*8 wid(21)/0.00002,0.00002,0.00002,0.00002,0.00002,0.00125,       !!! excited state widths in GeV
     &     0.00002,0.00002,0.00002,0.0002,0.00475,0.009,0.012,0.00002,
     &     0.00002,0.00002,0.0002,0.00035,0.00015,0.004,0.009  /
            
      
      REAL*8 x,tot,epnuc,F1,F2,FL

      q2 = q2p
      q2 = max(q2,0.0D0)
      qv2 = q2+nu*nu
      
c      smwid = 0.0035            !!! GeV, Use for Barreau data
      smwid = 0.00048  !!! Use for Yamaguchi data

      width = sqrt(smwid*smwid+wid(state)*wid(state))
      
      norm = width*sqrt(pi) 

c      nuel = 0.0
c      nuel = ((0.931494*A)**2+q2p-mp*mp)/2./mp
      nuel = q2/2./(0.931494*A)
      
      F1 = 0.0
      F2 = 0.0
      FL = 0.0
     
      x = q2/2.0/mp/nu
      

      q2s = q2
      
      q2f = q2/0.1975/0.1975

      
      
      if(state.EQ.18) then 
         exc(18) = min(0.0194+0.00016*sqrt(q2f),0.01955)
      endif

      nuex =  nuel+exc(state)
      
    
      
c      write(6,*) e,
      
      
      fs = exp(-1.0*(nu-nuex)**2./width/width)  !!! Gaussian probability based on distance from peak
      fs = fs/norm              !!! normalized
      if(((nu-nuex)/width).GT.5.0) fs = 0.0

      FF = 0.0
      FFT = 0.0
      if(state.LE.13) then
        call NUCFFS12C(A,Z,qv2,state,FF) !!! longitudinal FFs
      else
        call NUCFFS12CT(A,Z,qv2,state,FFT) !!! transverse FFs
      endif  

      
c      W2  = (Z*FF)**2+Z*Z*FFT*FFT

      W1 = 0.5*(Z*FFT)**2*fs

      WL = (Z*FF)**2*fs
      

c      F2 = nu*W2*fs
      F1 = mp*W1
      FL = q2*q2/qv2/nu*WL
c      nuccs = fs

c      write(6,*) state,q2s,q2vs,nuex

      
      
c      F2 = fs*F2
c      F1 = fs*F1

c      FL = (1.0+q2/nu/nu)*F2-2.0*x*F1

c      FL = fs*FL
      
c      write(6,*) A,Z,state,q2s,nu,nuel,nuex,FF,FFT
      

 
 1001 format(4F9.3,3E11.3)
      
      return

      end
