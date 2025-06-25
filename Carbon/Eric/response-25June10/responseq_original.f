      PROGRAM RESPONSEq

      IMPLICIT NONE

      real*8 Z, A, Q2, W2, xb, qv, nu, dnu, F1, FL, RT, RL, RTE, RLE
      real*8 nuel, ex, RTQE, RLQE, RTIE, RLIE, RTNS, RLNS, RTTOT, RLTOT
      real*8 flNS, f1NS, fLt, f1t, mp/0.938273/
      integer i,j,type
      real*8 xvalc(60) /     
     & 0.83076E-01,0.20000E+02,0.14359E+00,0.75694E+01,0.78533E+00,
     & 0.19341E+00,0.87115E+01,0.16646E+01,0.47764E+00,-.33752E+01,
     & 0.10033E+01,0.10834E+01,0.10161E+01,0.99018E+00,0.98696E+00,
     & 0.10113E+01,0.97623E+00,0.10009E+01,0.98933E+00,0.99457E+00,
     & 0.10187E+01,0.99589E+00,0.10023E+01,0.10143E+01,0.99514E+00,
     & 0.10992E+01,0.20094E+02,-.22835E+00,0.34511E+01,0.25992E+01,
     & 0.10000E-03,0.59395E+01,0.27600E+01,0.32563E+00,0.24997E+00,
     & 0.13790E-01,0.26621E+00,0.30937E-01,0.68585E-01,0.27873E+00,
     & -.84553E-01,0.10137E+00,0.24202E+00,0.15765E-01,0.11633E-01,
     & 0.10231E+01,0.10986E+01,0.11014E+01,0.11632E+01,0.99299E+00,
     & 0.10005E+01,0.10000E+01,0.10402E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01 /
      
      
      A = 27.0
      Z = 13.0

      read(5,*) qv
      
      
      dnu = 0.00015
      nu = 0.0
      
      do i=1,3000
        nu = nu + dnu

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

        if(ex.LE.0.012) then  !!! Only needed for plotting purposes
           RTNS = RTNS/6.0
           RLNS = RLNS/6.0
        endif
        
c        RLTOT = RLTOT+RLNS
c        RTTOT = RTTOT+RTNS

        
        if(q2.GT.0.0) 
     &       write(6,2000) qv,q2,ex,nu,RTTOT,RLTOT,RTQE,RLQE,RTIE,RLIE,
     &                      RTE,RLE,RTNS,RLNS          
 
        
      enddo

 2000  format(4f9.5,10E11.3)
      

      return
      end


      
      
      
      
     
CCC-----------------

      
