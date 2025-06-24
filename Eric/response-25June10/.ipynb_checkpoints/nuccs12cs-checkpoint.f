      SUBROUTINE NUCCS12CS(Z,A,e,ep,th,state,nuccs)
CCCCC  Calculates the nuclear elastic and excited state cross sections for 12C    CCCCC
CCCCC  This version smears in nu with a gaussian of width defined below           CCCCC 
CCCCC  Cross section is dsig/domega/dE in ub                                      CCCCC
      
      IMPLICIT NONE
      REAL*8 q2,q2f,q2v,q2vs,q2s,e,ep,th,signe,ff,fft,cs,W1,W2,tau,mott
      REAL*8 nuccs,sin2,cos2,tan2,recoil,width,smwid,norm,nu,nuex,nuel
      REAL*8 fs,f2,omega,WL2,WT2
      INTEGER i,j,IZ,IA,state
      REAL*8 Z,A,mp/0.93827/,hbarc/1973.28/,alpha/7.29735E-03/
      REAL*8 pi/3.14159/,radcon/0.0174533/,mA
      REAL*8 exc(22)/0.0,0.00444,0.00765,0.00964,0.01084,0.0137,0.0151,    !!! excited state energies in GeV
     &     0.0161,0.0183,0.020,0.0230,0.0315,0.042,0.0151,0.0161,0.0166,
     &     0.0181,0.0193,0.0206,0.0235,0.0315,0.01271/
      REAL*8 wid(22)/0.00002,0.00002,0.00002,0.00002,0.00002,0.00125,       !!! excited state widths in GeV
     &     0.00002,0.00002,0.00002,0.0002,0.00475,0.009,0.012,0.00002,
     &     0.00002,0.00002,0.0002,0.00035,0.00015,0.004,0.009,0.0002 /
      
      REAL*8 x,tot,epnuc


c      smwid = 0.0035   !!! GeV, Use for Barreau data 
c      smwid = 0.00048  !!! Use for Yamaguchi data
c      smwid = 0.00085  !!! Use for Ryan data
c       smwid = 0.00018           !!! Bates
c       smwid = 0.00025   !!! Yamaguchi fine binning
c       smwid = 0.00046  !!! Use for LEDEX 0.685 GeV
c      smwid = 0.00041  !!! Use for LEDEX 0.362 GeV
c       smwid = 0.000625  !!! Use for Crannell 250 MeV  (0.25%)
c      smwid = 0.0015    !!!  Use for Crannell 600 MeV (0.25%)

       smwid = 0.001   !!!  New value for Jan05

c      smwid = 0.00025   !!! Yamaguchi hi res
      
      width = sqrt(smwid*smwid+wid(state)*wid(state))
      
      norm = width*sqrt(pi) 
      
      sin2 = dsin(radcon*th/2.0)
      sin2 = sin2*sin2
      cos2 = 1.0-sin2
      tan2 = sin2/cos2
      tau = q2/4.0/mp

      nu = e-ep
      q2 = 4.0*e*ep*sin2
      q2v = q2+nu*nu

c     epnuc = A*mp*e/1.0079/(mp*A/1.0079+2.0*e*sin2)

      epnuc = A*0.931494*e/(A*0.931494+2.0*e*sin2)  !!! Equivalent to above  !!!
      
      nuel = e-epnuc


      
c       mA = A*0.931494
c       epnuc = mA*e/(mA+2.0*e*sin2)
c       nuel = e-epnuc
     
           
      nuex =  nuel+exc(state)
c      omega = max(0.0,nu-nuex)
      q2f = q2/0.1975/0.1975 
      
      if(state.EQ.18) then 
        exc(18) = min(0.0193+0.00015*sqrt(q2f),0.0195)
        nuex =  nuel+exc(state)
      endif
      Q2s = 4.*e*(e-nuex)*sin2  !!!! Q^2 at the E' for the given state
      q2vs = q2s+nuex*nuex
      
      fs = exp(-1.0*(nu-nuex)**2./width/width)  !!! Gaussian probability based on distance from peak
      fs = fs/norm              !!! normalized
      if(((nu-nuex)/width).GT.4.0) fs = 0.0
 

      FF = 0.0
      FFT = 0.0
      if(state.LE.13) then
        call NUCFFS12C(A,Z,Q2Vs,state,FF) !!! longitudinal FFs
      else
        call NUCFFS12CT(A,Z,Q2Vs,state,FFT) !!! transverse FFs
      endif
      
      
c      if(q2.GE.0.3) then
c        FF = 0.0  
c        FFT = 0.0
c      endif
        

      
      WL2  = (Z*FF)**2
      WT2 = 0.0
      WT2 = (Z*FFT)**2
      

      mott = 0.3894e3*alpha**2/4.                               
      mott = mott*cos2/e/e/sin2/sin2
      recoil = A*mp/1.007276/(12.0*mp+e*(1.0-dcos(radcon*th)))
      
      nuccs = 1000.0*mott*recoil*(q2*q2/q2v/q2v*WL2+
     &          (q2/2.0/q2v+tan2)*WT2)                     

      nuccs = fs*nuccs

     
c      write(6,*) "TEST"
      
 1001 format(4F9.3,3E11.3)
      
      return

      end
