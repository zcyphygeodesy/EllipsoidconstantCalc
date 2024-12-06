      subroutine Normalconst(GRS,BLH,NFD)
      !NFD(2)-the normal geopotential, normal gravity
      implicit none
 	real*8::BLH(3),NFD(2),val(3),rln(3),pn(40),dp1(40),dp2(40)
      integer::i,j,n
	real*8::GRS(6),djn(80),pi,RAD,gm,ae,ww,rr,t,u,vdn(5),atr
!---------------------------------------------------------------------
      call normdjn(GRS,djn)
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
	gm=GRS(1);ae=GRS(2);ww=GRS(4);vdn=0.d0
      call BLH_RLAT(GRS,BLH,rln)
      rr=rln(1);t=dsin(rln(2)*RAD);u=dcos(rln(2)*RAD)
      call LegPn_dt2(pn,dp1,dp2,40,t)
      do n=1,20
        atr=dexp(dble(2*n)*dlog(ae/rr))*djn(2*n)
        vdn(1)=vdn(1)+atr*pn(2*n)
        vdn(2)=vdn(2)+atr*pn(2*n)*dble(2*n+1)
        vdn(3)=vdn(3)+atr*dp1(2*n)
      enddo
      val(1)=gm/rr*(1.d0-vdn(1))+(ww*rr*u)**2/2.d0
      val(2)=-gm/rr**2*(1.d0-vdn(2))+ww**2*rr*u**2 
      val(3)=-gm/rr**2*vdn(3)+ww**2*rr*u*t         
      val(2)=dsqrt(val(2)**2+val(3)**2)
      NFD(1)=val(1);NFD(2)=val(2)
      end
