!  EllipsoidconstCalc.f90 
!
!  FUNCTIONS:
!  EllipsoidconstCalc - Entry point of console application.
!
!****************************************************************************

      program EllipsoidconstCalc
      implicit none
      integer::i,n
	real*8::GRS(6), a, b, gm, j2, f, U0, w, R, E, M, c, e2, e21,pm,q0,c20
	real*8::fa, ga, gb,BLH(3),val(2),rslt(18)
!---------------------------------------------------------------------
      GRS=1.d0 !>0 required
      GRS(1)= 3.986004415d14; GRS(2)=6378136.3d0; GRS(4) = 7.292115d-5 !3 known constants
      GRS(3)=1.082636277388d-3 !given j2 or
      !GRS(5)=1.d0/298.25641153d0 !given 1/f
      call ELLIPSOIDPARA(GRS)!calculate the ellipsoid constants
	BLH(1) = 0.0; BLH(2) = 110.0; BLH(3) = 0.0 !at the ellipsoidal Equator
	call NORMALCONST(GRS,BLH,val);ga = val(2);
	BLH(1) = 90.0; BLH(2) = 110.0; BLH(3) = 0.0 !at the ellipsoidal north pole
	call NORMALCONST(GRS,BLH,val);gb = val(2);
	gm = GRS(1); a=GRS(2); j2=GRS(3); w=GRS(4); f=GRS(5); U0=GRS(6)
      c20=-j2/dsqrt(5.d0)
	b = (1.d0 - f)*a; R = (a*a*b)**(1.d0 / 3.d0); E = dsqrt(a*a - b*b);
	e2 = 2.d0*f - f * f; e21 = 2.d0*f + 3.d0*f*f;
	M = a * (1.d0 - e2); c = a * a / b;
	pm = w * w*b / gm * a*a;
	fa = b*gb / a / ga - 1.d0;
      write (*,'(a60,f12.9)')'Geocentric gravitational constant GM(e14m2/s2) of the Earth',gm*1.d-14
      write (*,'(a60,f9.6)')'Mean angular velocity w(e-5/s) of the Earth',w*1.d5
      write (*,'(a60,f11.2)')'Major semi axis a(m) of the Earth',a
      write (*,'(a60,f15.12)')'Dynamic form factor J2(e-3):',j2*1.d3
      write (*,'(a60,f16.12)')'Geopotential coefficient c20(e-3) from GM',c20*1.d3
      write (*,'(a60,f13.8)')'Reciprocal 1/f of ellipsoid flattening',1.d0/f
      write (*,'(a60,f14.4)')'Normal ellipsoid geopotential U0=WG',U0
      !-----------------------------------------------
      write (*,'(a60,f13.4)')'Minor semi axis of the Earth b(m)',b
      write (*,'(a60,f13.4)')'Radius of sphere of same volume R(m)',R
      write (*,'(a60,f12.4)')'Linear eccentricity E(m)',E
      write (*,'(a60,f21.18)')'Square of first eccentricity e2',e2
      write (*,'(a60,f21.18)')'Square of second eccentricity e21',e21
      write (*,'(a60,f13.4)')'Equatorial curvature radius M(m)',M
      write (*,'(a60,f13.4)')'Polar radius of curvature c(m)',c
      !-----------------------------------------------
      write (*,'(a60,f15.10)')'Gravity flattening reciprocal 1/fk',1.d0/fa
      write (*,'(a60,f13.4)')'Geodetic parameter m',m
      write (*,'(a60,f13.10)')'Normal gravity at equator ga(m/s2)',ga
      write (*,'(a60,f13.10)')'Normal gravity at pole gp(m/s2)',gb
      pause
      end
