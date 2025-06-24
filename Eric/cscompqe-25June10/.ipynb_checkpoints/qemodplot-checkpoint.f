      PROGRAM QEMODPLOT

      IMPLICIT none

 
      real*8 x,q2,w2,A,Z,f1mec,f1,f1qe,r,rqe,f2,f2qe,rq3,mp,mp2,rat
      real*8 sum,sum2,sumrat,w2min,w2max,e,eb,theta,ep,epmax,nu,q2c
      real*8 pi,pi2,alpha,flux,sigL,sigT,sigm,eps,fL,fLqe,kappa,de
      real*8 sigtmec, sigmec,sigqe,sigtqe,siglqe,sigie,ebcc,veff,foc
      real*8 epnuc,nuel,nuccs,nuccstot,ev,epv,q2v,w2v,epsv,fluxv
      real*8 sigmnonuc,sin2,cos2,tan2,radcon/0.0174533/
      real*8 sigtot,signonuc
      
      integer i,j,k,state,type
      logical thend/.false./
      LOGICAL GOODFIT/.true./,coulomb/.true./
      real*8 xvalc(45) /          
c     & 0.76003E-01,0.84447E+01,0.99274E-01,0.65625E+01,0.69897E+00,
c     & 0.70012E+00,0.20269E+01,0.19072E+01,0.69159E+00,-.39652E+01,
c     & 0.98981E+00,0.97391E+00,0.10337E+01,0.99053E+00,0.10000E+01,
c     & 0.10040E+01,0.97073E+00,0.10072E+01,0.98901E+00,0.99552E+00,
c     & 0.10000E+01,0.99763E+00,0.10039E+01,0.10190E+01,0.10084E+01,
c     & 0.81286E+00,0.00000E+00,-.96577E+00,0.89281E+00,0.18759E+01,
c     & 0.21825E+01,0.21823E+01,0.28339E+01,0.49261E+00,0.22500E+00,
c     & 0.85763E-02,0.26821E+00,0.39391E-01,0.73511E-01,0.12573E+00,
c     & 0.70756E-01,0.96280E-01,0.27874E+00,0.64755E+00,0.87647E-01 /
c     & 0.89939E-01,0.10201E+02,0.13086E+00,0.68400E+01,0.76771E+00,
c     & 0.75766E+00,0.20269E+01,0.20714E+01,0.70139E+00,-.41708E+01,
c     & 0.98622E+00,0.97207E+00,0.10339E+01,0.98952E+00,0.10000E+01,
c     & 0.10043E+01,0.96991E+00,0.10074E+01,0.98893E+00,0.99523E+00,
c     & 0.10000E+01,0.99760E+00,0.10040E+01,0.10190E+01,0.10083E+01,
c     & 0.81790E+00,0.00000E+00,-.96381E+00,0.91077E+00,0.18259E+01,
c     & 0.21899E+01,0.21906E+01,0.27930E+01,0.48048E+00,0.22800E+00,
c     & 0.10571E-01,0.26830E+00,0.38317E-01,0.71332E-01,0.75446E-01,
c     & 0.74997E-01,0.81582E-01,0.30230E+00,0.66323E+00,0.99925E+01 /
c     & 0.92048E-01,0.11466E+02,0.13563E+00,0.69188E+01,0.76870E+00,
c     & 0.69650E-01,0.20269E+01,0.20653E+01,0.70040E+00,-.41622E+01,
c     & 0.99256E+00,0.98228E+00,0.10238E+01,0.99270E+00,0.10000E+01,
c     & 0.10032E+01,0.97868E+00,0.10045E+01,0.99293E+00,0.99675E+00,
c     & 0.10000E+01,0.99840E+00,0.10025E+01,0.10110E+01,0.10052E+01,
c     & 0.80156E+00,0.00000E+00,-.98951E+00,0.89455E+00,0.19710E+01,
c     & 0.21378E+01,0.21390E+01,0.28917E+01,0.48810E+00,0.22800E+00,
c     & 0.10979E-01,0.27037E+00,0.40631E-01,0.74907E-01,-.29387E-01,
c     & 0.19512E-01,0.81638E-01,0.30699E+00,0.51130E+00,-.23686E-01 /
c     & 0.94449E-01,0.10821E+02,0.13888E+00,0.69224E+01,0.78078E+00,
c     & 0.11235E+00,0.20269E+01,0.20732E+01,0.70846E+00,-.41558E+01,
c     & 0.99233E+00,0.98152E+00,0.10240E+01,0.99265E+00,0.10000E+01,
c     & 0.10058E+01,0.97847E+00,0.10044E+01,0.99271E+00,0.99628E+00,
c     & 0.10000E+01,0.99844E+00,0.10026E+01,0.10105E+01,0.10050E+01,
c     & 0.78676E+00,0.00000E+00,-.10317E+01,0.92242E+00,0.20365E+01,
c     & 0.19625E+01,0.19652E+01,0.28451E+01,0.50889E+00,0.22800E+00,
c     & 0.11420E-01,0.27154E+00,0.40438E-01,0.75126E-01,0.44351E-01,
c     & 0.71092E-02,0.78423E-01,0.29565E+00,0.62607E+00,-.13227E+00 /
     & 0.73010E-01,0.10612E+02,0.11337E+00,0.67674E+01,0.74110E+00,
     & 0.76437E-01,0.87115E+01,0.18999E+01,0.75048E+00,-.38472E+01,
     & 0.99281E+00,0.98412E+00,0.10303E+01,0.10018E+01,0.10000E+01,
     & 0.10072E+01,0.97384E+00,0.10063E+01,0.98851E+00,0.99300E+00,
     & 0.10000E+01,0.99661E+00,0.10022E+01,0.10116E+01,0.10039E+01,
     & 0.78136E+00,0.11295E-05,-.10064E+01,0.93493E+00,0.22110E+01,
     & 0.26555E+01,0.26555E+01,0.34117E+01,0.76660E+00,0.21000E+00,
     & 0.65674E-02,0.25672E+00,0.31976E-01,0.56458E-01,-.86596E-01,
     & 0.27536E-01,0.75467E-01,0.31597E+00,0.99112E-01,0.96124E-02 /      
      
      
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1./137.

      A = 12.
      Z = 6.

      read(5,*) e,theta

      
      sin2 = dsin(radcon*theta/2.0)
      sin2 = sin2*sin2
      cos2 = 1.0-sin2
      tan2 = sin2/cos2
 
      nuel = e*e*sin2      
      nuel = nuel/(8.0/1.00797*mp+e*2)
      epnuc = e-nuel
      epmax = epnuc-0.000

      de = 0.001
      if(e.LT.0.8) then
         de = 0.0005
      elseif(e.LT.1.4) then
         de = 0.0007
      elseif(e.LT.1.7) then
        de = 0.0009
      elseif(e.LT.3.0) then
        de = 0.0025
      endif   
      do i=1,2000
        ep = epmax-i*de
        q2 = 4.*e*ep*sin2
        w2 =  mp2+2.*mp*nu-q2
        if(coulomb) then
         call vcoul(A,Z,veff)
         foc = 1.0D0 + veff/e
         ev =  e + veff     
         epv = ep + veff
        endif
         
        nu = e-ep
        q2v = 4.*ev*epv*sin2
        epsv = 1./(1. + 2.*(nu*nu+q2v)/q2v*tan2)
        w2v = mp2+2.*mp*nu-q2v
        x = q2v/(w2v-mp2+q2v)
        kappa = abs(w2v-mp2)/2./mp
        fluxv = alpha*kappa/(2.*pi2*q2v)*epv/ev/(1.-epsv)

        type = 1
        do type=1,4
         call csfitcomp(w2v,q2v,A,Z,xvalc,type,sigt,sigL)      
         sigm = fluxv*(sigt+epsv*sigl)
         sigm =  0.3894e3*8.0d0*pi2*alpha/abs(w2v-mp2)*sigm
         sigm = sigm*foc*foc
         sigm = sigm/A
         if(type.eq.1) then
            sigtot = sigm
         elseif(type.eq.2) then
           sigqe = sigm
         elseif(type.eq.3) then
           sigie = sigm
         elseif(type.eq.4) then
           sigmec = sigm
         endif
        enddo  

 

        nuccstot = 0.0
        do k=2,21
           call nuccs12cs(Z,A,ev,epv,theta,k,nuccs)
c           if(q2.GT.0.3) nuccs = 0.0
           nuccstot = nuccstot+nuccs/1000.0
           nuccstot = nuccstot
        enddo
        nuccstot = nuccstot/A
        signonuc = sigtot
        sigtot = sigtot + nuccstot
        
         

         if(ep.GT.0.01.AND.w2.LT.40.) then
            write(6,2000) ep,theta,nu-nuel,w2,q2,sigtot,sigqe,sigie,
     &         sigmec,nuccstot,signonuc 
     
         endif

      enddo


 2000 format(5f9.4,6f12.4)
      end


 









