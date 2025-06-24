      SUBROUTINE CSFITCOMP(w2,q2,A,Z,XVALC,type,sigt,sigL)
      IMPLICIT none

      real*8 e,ep,th,q2,w2,x,cs,flux,kappa2,sin2,tan2,csmod
      real*8 f1,f2,fl,f1qe,f2qe,flqe,f1mec,f2mec,fLmec,f1i,f2i,fLi
      real*8 r,rqe,sigt,sigl,sigm
      real*8 alpha,pi,pi2,mp,mp2,res,veff,foc,Z,A,xvalc(45)
      real*8 psip,psimax,psimin,fy1,fy2,int1,int2,rat,f1t,f2t,fLt,dpsi
      integer i,j,k,ntot,nbins,type
      LOGICAL GOODFIT/.true./  
      character*40 filename

      psimin = -2.3
      psimax = 5.0
      nbins = 220

      mp = .938272
      mp2 = mp*mp

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

      x = q2/abs(w2-mp2+q2)
      
      dpsi = (psimax - psimin)/float(nbins)

      kappa2 = (1.+4.*x*x*mp2/q2)   

      int1 = 0.0D0
      int2 = 0.0D0
      rat = 1.0D0

CCC  NEXT bit only needed if fitting scaling function CCC
      do i=1,nbins       
        psip = psimin+dpsi*(i-1)            
       FY1 = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
       FY2 = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     > (1. + exp(xvalc(10) * psip))*(1.0-abs(psip)/psimax)**2.0* 
     > (1.0+abs(psip)/psimax)**2.0

       if(psip.GT.psimax) FY2 = 0.0
       FY2 = max(0.0,FY2)
       int1 = int1+fy1  
       int2 = int2+fy2
      enddo

      rat= 1.00/int2/dpsi

CCC   


c      write(6,*) int1*0.04,int2*0.04,rat

c      call f1f2in09(Z,A,q2,w2,xvalc,f1t,f2t,r)

      call gsmearing(Z,A,w2,q2,xvalc,f1i,f2i,fLi)
      if(fLi.LT.0.0) FLi = 0.0

      
c      write(6,*) "gsmearing:  ", w2,q2, f1,f2,fL

c      call smearing(Z,A,w2,q2,xvalc,f1,f2,fL)
      
c      write(6,*) "smearing:  ",w2,q2, f1,f2,fL

      r = fL/2.0D0/x/f1

c      f1 = f1
c      f2 = f2
c      fL = fL

c      if(w2.GT.1.5) write(6,2000) w2,q2,f1,f1t,f2,f2t


      f1qe = 0.0
      f2qe = 0.0
      call qenuc21off(Z,A,q2,w2,xvalc,f1qe,f2qe)
      
      f1qe = f1qe*rat  !!! renormalize
      f2qe = f2qe*rat  !!! renormalize
      fLqe =  kappa2*f2qe-2.0*x*f1qe 
      if(fLqe.LT.0.0) flqe = 0.0
      

      f1mec = 0.0
      fLmec = 0.0
      call MEC2021(Z,A,w2,q2,xvalc,f1mec)
     
      f2mec = 2.0*x*f1mec/kappa2
      


      if(type.EQ.1) then
        f1 = f1i + f1qe + f1mec
        f2 = f2i + f2qe + f2mec
c     fL = kappa2*f2-2.0*x*f1
        fL = fLi+fLqe
      elseif(type.EQ.2) then
        f1 = f1qe
        f2 = f2qe
        fL = fLqe
      elseif(type.EQ.3) then
        f1 = f1i
        f2 = f2i
        fL = fLi
      elseif(type.EQ.4) then         
        f1 = f1mec
        f2 = f2mec
        fL = fLmec
      elseif(type.EQ.5) then         
        f1 = f1qe+f1mec
        f2 = f2qe+f2mec
        fL = fLqe+fLmec            
      endif
      if(fL.LT.0.0) fL = 0.0   
       
       sigt = f1
       sigl = fL/2./x

c      write(6,*) "2  ",w2,q2,x,sigt,sigl

 2000 format(6f10.4)
          
      return

      end
      




