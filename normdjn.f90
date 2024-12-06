      subroutine normdjn(GRS,djn)
!输入GRS－gm,ae,j2,omega,1/f
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     djn  array of even zonals of normal ellipsoid
!c     ex. djn(2)=j2=GRS(3)
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      real*8 GRS(6),djn(80)
	ae=GRS(2)
	ff=GRS(5)
	do 101 n2=1,80,1
	djn(n2)=0.d0
  101	continue
	djn(2)=GRS(3)
      esq=(2.d0-ff)*ff!e**2
!c     compute normal even zonals from 4 to 20
      do 200 n2=4,80,2
        n=n2/2
        djn(n2)=(-1.d0)**(n+1)*3.d0*esq**n/(n2+1.d0)/(n2+3.d0)
     >        *(1.d0-dble(n)+5.d0*dble(n)*djn(2)/esq)
  200 continue
      return
      end
