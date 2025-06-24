      PROGRAM response_q2_w2

      IMPLICIT NONE

      real*8 Z, A, Q2, W2, xb, qv, nu, dnu, F1, FL, RT, RL, RTE, RLE
      real*8 nuel, ex, RTQE, RLQE, RTIE, RLIE, RTNS, RLNS, RTTOT, RLTOT
      real*8 flNS, f1NS, fLt, f1t, mp/0.938273/
      integer i,j,type
      integer io_status, arg_status, unit
      character(len=30) filename
      real*8 xvalc(45) /     
c     & 0.94187E-01,0.10367E+02,0.13523E+00,0.68783E+01,0.76964E+00,
c     & 0.74171E+00,0.20269E+01,0.20269E+01,0.69603E+00,-.41246E+01,
c     & 0.98677E+00,0.97242E+00,0.10338E+01,0.98935E+00,0.10000E+01,
c     & 0.10041E+01,0.97018E+00,0.10073E+01,0.98889E+00,0.99529E+00,
c     & 0.10000E+01,0.99757E+00,0.10040E+01,0.10191E+01,0.10083E+01,
c     & 0.81949E+00,0.00000E+00,-.95810E+00,0.90514E+00,0.18136E+01,
c     & 0.21799E+01,0.21799E+01,0.27720E+01,0.47927E+00,0.22800E+00,
c     & 0.10133E-01,0.26823E+00,0.38365E-01,0.71421E-01,0.88152E-01,
c     & 0.74219E-01,0.82864E-01,0.29490E+00,0.66002E+00,0.81059E+01 /
c     & 0.94131E-01,0.10746E+02,0.13847E+00,0.69187E+01,0.78080E+00,
c     & 0.11348E+00,0.20269E+01,0.20737E+01,0.70874E+00,-.41560E+01,
c     & 0.99234E+00,0.98155E+00,0.10240E+01,0.99269E+00,0.10000E+01,
c     & 0.10058E+01,0.97857E+00,0.10044E+01,0.99274E+00,0.99636E+00,
c     & 0.10000E+01,0.99844E+00,0.10026E+01,0.10106E+01,0.10050E+01,
c     & 0.78561E+00,0.00000E+00,-.10365E+01,0.92719E+00,0.20435E+01,
c     & 0.19717E+01,0.19779E+01,0.28591E+01,0.50813E+00,0.22800E+00,
c     & 0.11486E-01,0.27182E+00,0.40453E-01,0.75158E-01,0.45838E-01,
c     & 0.68540E-02,0.78288E-01,0.29602E+00,0.57452E+00,-.13525E+00 /
c     & 0.86448E-01,0.12116E+02,0.12905E+00,0.68753E+01,0.76463E+00,
c     & 0.76437E-01,0.87115E+01,0.19440E+01,0.67375E+00,-.40009E+01,
c     & 0.99298E+00,0.98302E+00,0.10302E+01,0.10011E+01,0.10000E+01,
c     & 0.10073E+01,0.97422E+00,0.10058E+01,0.98887E+00,0.99433E+00,
c     & 0.10000E+01,0.99588E+00,0.10026E+01,0.10120E+01,0.10046E+01,
c     & 0.79586E+00,0.11295E-05,-.97614E+00,0.92788E+00,0.20215E+01,
c     & 0.24889E+01,0.24890E+01,0.31432E+01,0.74838E+00,0.22800E+00,
c     & 0.70982E-02,0.25749E+00,0.31475E-01,0.56800E-01,-.13645E+00,
c     & 0.36333E-01,0.75341E-01,0.27892E+00,0.15692E+00,0.10459E+00 /
     & 0.91648E-01,0.12714E+02,0.13380E+00,0.69068E+01,0.77023E+00,
     & 0.76437E-01,0.87115E+01,0.18976E+01,0.66472E+00,-.39215E+01,
     & 0.99320E+00,0.98312E+00,0.10302E+01,0.10009E+01,0.10000E+01,
     & 0.10070E+01,0.97472E+00,0.10059E+01,0.98892E+00,0.99434E+00,
     & 0.10000E+01,0.99596E+00,0.10028E+01,0.10122E+01,0.10045E+01,
     & 0.79845E+00,0.11295E-05,-.97071E+00,0.92502E+00,0.20146E+01,
     & 0.24416E+01,0.24499E+01,0.31154E+01,0.72998E+00,0.22800E+00,
     & 0.76502E-02,0.25718E+00,0.31429E-01,0.58780E-01,-.15059E+00,
     & 0.38790E-01,0.77051E-01,0.26795E+00,0.17673E+00,0.10451E-01 /
      
      
      A = 12.0
      Z = 6.0
      call get_command_argument(1, filename, arg_status)
      unit = 20
      
      open(UNIT=unit, FILE=filename, STATUS='old', IOSTAT=io_status)
      
      if (io_status/= 0) then
        print *, 'Unable to open file:',filename
        stop
      endif


      
      
      i = 0
      do 
        read(unit,*,IOSTAT=io_status) i, q2, w2
        if (io_status /= 0) exit
        
        nu = (w2+q2-mp*mp)/(2.0*mp)
        qv = sqrt(q2+nu*nu)
        nuel = q2/2./(0.931494*A)


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

c        if(ex.LE.0.012) then  !!! Only needed for plotting purposes
c           RTNS = RTNS/6.0
c           RLNS = RLNS/6.0
c        endif

        if(RLNS.LE.1E-40) RLNS = 0.0
        if(RTNS.LE.1E-40) RTNS = 0.0
        
        RLTOT = RLTOT+RLNS
        RTTOT = RTTOT+RTNS
        !write(6,2000) i, RTTOT, RLTOT
        write(6, 2000) i, RTTOT, RLTOT, (RTQE+RTIE+RTE), (RLQE+RLIE+RLE) 
        
c        if(q2.GT.0.0) 
c     &       write(6,2000) i,qv,q2,ex,nu,RTTOT,RLTOT,RTQE,RLQE,RTIE,RLIE,
c     &                      RTE,RLE,RTNS,RLNS        
        
c        if(q2.GT.0.0) 
c     &       write(6,2000) qv,q2,ex,nu,RTTOT,RLTOT,RTQE,RLQE,RTIE,RLIE,
c     &                      RTE,RLE,RTNS,RLNS   

c         if(q2.LE.0.0) 
c     &       write(6,2000) qv,q2,ex,nu,RTTOT,RLTOT,RTQE,RLQE,RTIE,RLIE,
c     &                      RTE,RLE,RTNS,RLNS   
c        i = i + 1

 
        
      enddo
c 2000  format(4f9.5,10E15.7)
c 2000  format(4f9.5,10E11.3)
 2000  format(1I5,4E15.7)      
c 2000  format(1I5,2E15.7) 
      close(UNIT=unit)
      return
      end


      
      
      
      
     
CCC-----------------

      
