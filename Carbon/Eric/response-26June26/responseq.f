      PROGRAM RESPONSEQ

      IMPLICIT NONE

      real*8 Z, A, Q2, W2, xb, qv, nu, dnu, F1, FL, RT, RL, RTE, RLE
      real*8 nuel, ex, RTQE, RLQE, RTIE, RLIE, RTNS, RLNS, RTTOT, RLTOT
      real*8 flNS, f1NS, fLt, f1t, mp/0.938273/
      integer i,j,type
      real*8 xvalc(100) /     
     & 0.18408E+00,0.69756E+01,0.17765E+00,0.12545E+02,0.84099E+00,
     & 0.12120E+01,0.29829E+01,0.15953E+01,0.64928E+00,-.34309E+01,
     & 0.58794E+02,0.17722E+00,0.77061E+01,0.10000E+01,0.21511E+00,
     & 0.50422E+00,0.26811E+00,0.24094E+00,0.37362E+01,0.18371E-01,
     & 0.80561E-01,0.15184E+00,0.22200E+00,0.00000E+00,0.22200E+00,
     & 0.46401E-01,0.66384E+00,-.35015E-01,0.22589E-01,0.86312E-01,
     & 0.91300E+00,0.13026E+00,0.16315E+00,0.10000E+01,0.10472E-01,
     & 0.69558E+01,0.21098E+01,0.00000E+00,0.24792E+02,0.10000E-01,
     & 0.10000E+01,0.96116E+00,0.10726E+01,0.99999E+00,0.93939E+00,
     & 0.10094E+01,0.96949E+00,0.10255E+01,0.98498E+00,0.99732E+00,
     & 0.10414E+01,0.10000E+01,0.10211E+01,0.10350E+01,0.10290E+01,
     & 0.95411E+00,0.11187E+01,0.11900E+01,0.11745E+01,0.10000E+01,
     & 0.10012E+01,0.10687E+01,0.10269E+01,0.11001E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01 /
      
      
      A = 12.0
      Z = 6.0

      read(5,*) qv
      
      
      dnu = 0.0005
!       dnu = 0.005

      nu = 0.0
      
      do i=1,7000
        nu = nu + dnu
        if (nu>=qv) exit

        q2 = qv*qv-nu*nu
        nuel = q2/2./(0.931494*A)
        w2 = mp*mp+2.0*mp*nu-q2
        xb = q2/2.0/mp/nu

        ex = nu-nuel
        
        type = 1
        call csfitcomp(w2,q2,A,Z,XVALC,type,f1,fL) !!!  total response
        fL = 2.0*xb*fL
        RTTOT = 2.0/mp*F1/1000.0
        RLTOT = qv*qv/q2/2.0/mp/xb*FL/1000.

        type = 2
        call csfitcomp(w2,q2,A,Z,XVALC,type,f1,fL) !!!  QE response
        fL = 2.0*xb*fL
        RTQE = 2.0/mp*F1/1000.0
        RLQE = qv*qv/q2/2.0/mp/xb*FL/1000.0
        
        type = 3
        call csfitcomp(w2,q2,A,Z,XVALC,type,f1,fL) !!!  IE response
        fL = 2.0*xb*fL
        RTIE = 2.0/mp*F1/1000.0
        RLIE =  qv*qv/q2/2.0/mp/xb*FL/1000.0
        
        type = 4
        call csfitcomp(w2,q2,A,Z,XVALC,type,f1,fL) !!!  TE response
        fL = 2.0*xb*fL  
        RTE = 2.0/mp*F1/1000.0
        RLE = 0.0


c        write(6,*) RLTOT,RLIE+RLQE
        
        fLNS = 0.0
        f1NS = 0.0
        do j=2,22
           call nuc12sf(Z,A,nu,q2,j,f1t,fLt)

          fLNS = fLNS + fLt
          f1NS = f1NS + f1t      
        enddo
        RTNS = 2.0/mp*F1NS/1000.0 
        RLNS =  qv*qv/q2/2.0/mp/xb*FLNS/1000.0
        if(RLNS.LE.1E-40) RLNS = 0.0
        if(RTNS.LE.1E-40) RTNS = 0.0

        if(ex.LE.0.012) then  !!! Only needed for plotting purposes
c           RTNS = RTNS/6.0
c           RLNS = RLNS/6.0
        endif
        
        RLTOT = RLTOT+RLNS
        RTTOT = RTTOT+RTNS
        
        if(q2.GT.0.0) 
     &       write(6,2000) qv,q2,ex,nu,RTTOT,RLTOT,RTQE,RLQE,RTIE,RLIE,
     &                      RTE,RLE,RTNS,RLNS          
 
        
      enddo

 2000  format(4f9.5,10E11.3)
      

      return
      end


      
      
      

      
     
CCC-----------------

      
