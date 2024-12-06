## Fortran codes for calculation of the geometric and physical constants of Earth ellipsoid
https://www.zcyphygeodesy.com/en/h-nd-116.html
## [Algorithm purpose]
    From four basic parameters of the Earth ellipsoid, calculate the main geometric and physical derived constants of the Earth ellipsoid. The fourth basic parameter can be selected from the second-degree zonal harmonic coefficient C̅₂₀ from global geopotential model, dynamic form factor J₂, reciprocal 1/f of the ellipsoid flattening and ellipsoid normal geopotential U₀.
    PAGravf4.5 suggests that the scale parameters (GM, a) of global geopotential model, second-degree zonal harmonic coefficient C̅₂₀ and the mean rotation angular velocity ω should be employed as the four basic parameters of the normal ellipsoid. Using such a normal ellipsoid as the reference datum, the second-degree zonal harmonic term of anomalous gravity field is always zero, which is beneficial to improve the performance of the gravity field approach.
## [Main program for entrance]
    EllipsoidconstCalc.f90
    Given the dynamic form factor J2 in the test program.
    Output all the geometric and physical constants of the Earth ellipsoid on the computer scree, which include:
    Geocentric gravitational constant GM(e14m2/s2) of the Earth
    Mean angular velocity w(e-5/s) of the Earth
    Major semi axis a(m) of the Earth
    Dynamic form factor J2(e-3)
    Geopotential coefficient c20(e-3) from GM
    Reciprocal 1/f of ellipsoid flattening
    Normal ellipsoid geopotential U0=WG
    Minor semi axis of the Earth b(m)
    Radius of sphere of same volume R(m)
    Linear eccentricity E(m)
    Square of first eccentricity e2
    Square of second eccentricity e21
    Equatorial curvature radius M(m)
    Polar radius of curvature c(m)
    Gravity flattening reciprocal 1/fk
    Geodetic parameter m
    Normal gravity at equator ga(m/s2)
    Normal gravity at pole gp(m/s2)
## (1) Calculation module for the basic parameters of the Earth ellipsoid
    ELLIPSOIDPARA(GRS)
    Input parameters: GRS(1), GRS(2), GRS(4) and one of GRS(3), GRS(5), GRS(6).
    GRS(1) - Geocentric gravitational constant GM. 
    GRS(2) - Major semi axis of the Earth. 
    GRS(3)- Dynamic form factor. 
    GRS(4)- Mean angular velocity omega (e-5/s) of the Earth.
    GRS(5)- ellipsoid flattening f.
    GRS(6)- Normal ellipsoid geopotential.
    Return parameters: GRS(1:6)
## (2) Calculation module for the normal gravity potential and normal gravity
    Normalconst(GRS,BLH,NFD)
    Input：BLH(3)-Latitude, longitude (degree decimal) and ellipsoidal height (m) of the calculation point
    Return parameters: NFD(1)- the normal gravity potential. NFD(2)- normal gravity.
## (3) Calculation module for normal geopotential coefficients
    normdjn(GRS, djn)
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Algorithm module for transforming ellipsoid geodetic coordinates into spherical coordinates
    BLH_RLAT(GRS,BLH,RLAT)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.

DOS executable test file and all input and output data.
