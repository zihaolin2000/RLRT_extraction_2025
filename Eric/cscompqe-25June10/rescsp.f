CCC  Version 120521  -  Author:  M.E. Christy, christy@jlab.org             CCC
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
     & 0.12291E+01,0.15173E+01,0.15044E+01,0.17100E+01,0.16801E+01,
     & 0.14312E+01,0.12616E+00,0.23000E+00,0.92594E-01,0.90606E-01,
     & 0.75000E-01,0.35067E+00,0.75729E+01,0.56091E+01,0.94606E+01,
     & 0.20156E+01,0.66190E+01,0.41732E+00,0.23980E-01,0.53136E+01,
     & 0.63752E+00,0.11484E+02,0.69949E-01,0.26191E+01,0.53603E-01,
     & 0.65000E+02,0.15351E+00,0.20624E+01,0.23408E+01,0.16100E+02,
     & 0.62414E+02,0.17201E+01,0.23261E+00,0.65000E+02,0.23292E+01,
     & 0.14980E+01,0.23000E+00,0.63385E+00,0.19093E-01,0.61061E-01,
     & 0.29146E-02,0.54388E+00,0.77997E+00,0.28783E+00,0.10605E+01,
     & 0.69793E+00,0.20009E+01,0.57000E+00,0.41632E+01,0.38427E+00,
     & 0.10000E+01,0.99842E+00,0.98719E+00,0.10168E+01,0.98945E+00,
     & 0.99594E+00,0.98799E+00,0.10271E+01,0.10650E+01,0.97920E+00,
     & 0.10152E+01,0.99622E+00,0.81011E+01,0.10070E-02,0.14857E+01,
     & 0.33445E+01,0.31641E-09,0.69755E+02,0.55228E+01,0.14438E+00,
     & 0.60474E+01,0.65395E-07,0.14129E+01,0.58609E+00,0.36220E+01,
     & 0.92699E+00,0.14418E+01,0.86403E-02,0.10001E-03,0.75106E+00,
     & 0.76077E+00,0.42272E+00,0.55511E-11,0.52486E+00,0.58153E+00,
     & 0.15798E+01,0.50105E+00,0.89149E+02,0.72789E+00,0.24813E-01,
     & -.61906E+00,0.10000E+01,0.00000E+00,0.00000E+00,0.68158E+03,
     & 0.12429E+01,0.00000E+00,0.00000E+00,0.00000E+00,0.10000E-05 /
           
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


       




