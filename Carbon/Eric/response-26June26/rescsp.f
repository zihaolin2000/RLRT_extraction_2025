CCC  Version 061525  -  Author:  M.E. Christy, christy@jlab.org             CCC
CCC  This fit version includes data from a number of JLab Hall experiments  CCC
CCC  as well as DIS data from SLAC (L. Whitlow) and photoproduction data    CCC
CCC  from DAPHNE and older data sets.                                       CCC
CCC  Subroutine to get Transverse and Longitudinal eP cross sections        CCC  
CCC  from fits cross sections over a range of epsilon.  The subroutine      CCC
CCC  resmod.f is required.  Units are in ub/Sr/Gev.                         CCC
CCC  
CCC   Region of applicability has been extended to cover the full JLab      CCC
CCC   11 GeV kinematic range of Q^2 < 30 GeV^2 and W^2 < 20                 CCC


      SUBROUTINE rescsp(W2,Q2,sigT,sigL)

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),xval(100)
      real*8 mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

      data xval/
     & 0.12290E+01,0.15245E+01,0.15053E+01,0.17130E+01,0.16809E+01,
     & 0.14382E+01,0.12604E+00,0.23500E+00,0.89515E-01,0.85287E-01,
     & 0.75000E-01,0.37589E+00,0.75995E+01,0.55341E+01,0.62953E+01,
     & 0.58604E+00,0.67857E+01,0.24949E+01,0.43242E+00,0.17928E+01,
     & 0.49760E+00,0.22988E+02,0.20914E-01,0.92496E+00,0.29301E-01,
     & 0.10046E+03,0.22301E+00,0.49346E+00,0.23833E+01,0.52660E+01,
     & 0.14292E+02,0.42200E+00,0.53650E-01,0.10436E+03,0.16219E+00,
     & 0.55858E+00,0.27000E+00,0.32575E+00,0.20218E+01,0.96124E-01,
     & 0.71437E+00,0.11438E+01,0.50419E+00,0.16961E+01,0.36714E+00,
     & 0.51125E+00,0.19923E+01,0.57000E+00,0.42071E+01,0.71386E+00,
     & 0.10000E+01,0.99887E+00,0.99292E+00,0.10078E+01,0.98970E+00,
     & 0.99921E+00,0.98542E+00,0.10247E+01,0.10450E+01,0.99335E+00,
     & 0.10120E+01,0.99404E+00,0.77754E+01,0.32542E+01,0.18685E+01,
     & 0.37473E+01,0.00000E+00,0.93168E+02,0.57939E+01,0.10210E-02,
     & 0.31525E+01,0.73161E+01,0.22053E+01,0.53151E+00,0.32473E+01,
     & 0.13766E+01,0.14511E+01,0.12058E+01,0.21044E-03,0.79047E+00,
     & 0.78882E+00,0.71543E+00,0.95268E-01,0.10405E+02,0.18309E+01,
     & 0.00000E+00,0.13763E+02,0.18380E+02,0.56863E-01,0.16637E-01,
     & -.20137E+00,0.59193E+01,-.19076E+00,0.40936E+00,0.54781E+03,
     & 0.11113E+03,0.00000E+00,0.00000E+00,0.00000E+00,0.71670E-01/
      
      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 

        if(i.LE.12) xvalL(i) = xval1(i)
        if(i.EQ.47.OR.i.EQ.48) xvalL(i) = xval1(i)
      enddo

       
      xb = q2/(w2+q2-mp2)

      call resmodp(1,w2,q2,xval1,sigT)
      call resmodp(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
c      write(6,*) w2,q2,F1,FL,R

      end


       




