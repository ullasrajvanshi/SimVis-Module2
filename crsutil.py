"""
# ----------------------------------------------------------------------------
# crsutil.py    CRSUTIL Coordinate and Time Reference System Toolbox.
# Version 1.0 (1 November 2019).
# Created by: Hans van der Maarel and Ullas Rajvanshi
# Date:       13 Nov 2019
# Modified:   -
#
#   Copyright: Hans van der Marel, Ullas Rajvanshi, Delft University of Technology
#   Email:     H.vanderMarel@tudelft.nl ; U.Rajvanshi@student.tudelft.nl
#   Github:    www.github.com/ullasrajvanshi
# ----------------------------------------------------------------------------
# Functions created:
# Coordinate transformations (ECEF reference frame)
#   inqell      - Semi-major axis, flattening and GM for various ellipsoids
#   xyz2plh     - Cartesian Coordinates to Ellipsoidal coordinates
#   plh2xyz     - Ellipsoidal coordinates to Cartesian Coordinates
#   xyz2neu     - North, East, Up (dN, dE, dU) to Cartesian delta's (dX, dY, dZ)
#   neu2xyz     - Cartesian delta's (dX, dY, dZ) to North, East, Up (dN, dE, dU)
#   plh2neu     - Ellipsoidal (Lat,Lon,Hgt) to North,East,Up (dN, dE, dU)
#   xyz2zas     - Cartesian coordinates to Zenith angle, azimuth and distance
#   zas2xyz     - Zenith angle, azimuth and distance to cartesian coordinates
#
# UT1 to GMST, and ECI/ECEF, conversions
#
#   ut2gmst    - Compute Greenwich Mean Siderial Time from UT1
#   ecef2eci   - Convert position and velocity from ECEF to ECI reference frame
#   eci2ecef   - Convert position and velocity from ECI to ECEF reference frame
#
# Keplerian elements
#
#   vec2orb     - Convert inertial state vector into Keplerian elements
#   orb2vec     - Convert Keplerian elements into inertial state vector
#   kepler      - Compute mean anomaly from eccentric anomaly (Kepler's equation)
#   keplernu    - Compute mean anomaly from true anomaly (Kepler's equation)
#   keplerm     - Compute eccentric/true from mean anomaly solving Kepler's eqn
#
# ----------------------------------------------------------------------------

"""
# Importing the Libraries
import numpy as np

def inqell():
    ell = ['AIRY', 'BESSEL', 'CLARKE', 'INTERNATIONAL', 'HAYFORD', 'GRS80', 'WGS-84']

    par = np.array([
        [6377563.396, 299.324964, np.NaN],
        [6377397.155, 299.1528128, np.NaN],
        [6378249.145, 293.465, np.NaN],
        [6378388.0, 297.00, np.NaN],
        [6378388.0, 297.00, 3.986329e14],
        [6378137.0, 298.257222101, 3.986005e14],
        [6378137.0, 298.257223563, 3.986005e14]])

    i = 0

    for j in range(0, par.shape[0]):
        i = j
    if i == 0:
        i = par.shape[0]
    # module only taking WGS-84 at the moment by default
    # need to add a string comparison and based on that add the current parameter value
    a = par[i, 0]
    f = 1 / par[i, 1]
    GM = par[i, 2]

    return a, f, GM


def xyz2plh(xyz, ellipse='WGS-84', method=0):
    """
    #(c) Hans van der Marel, Delft University of Technology, 1995,2013
    :param xyz:Nx3 matrix XYZ with in the rows cartesian coordinates X, Y and Z
    :param ellipse: allows to specify the ellipsoid. ELLIPS is a text
    is a text string with the name of the ellipsoid or a vector with the
    semi-major axis a and flattening 1/f. Default for ellips is 'WGS-84'
    :param method:  uses the more conventional iterative method
    instead of Bowring's method (the default method). Bowring's method is
    faster, but can only be used on the surface of the Earth. The iterative
    method is slower and less precise on the surface of the earth, but should
    be used above 10-20 km of altitude.
    :return:  Nx3 matrix PLH with ellipsoidal coordinates
    Phi, Lambda and h. Phi and Lambda are in radians, h is in meters
    """

    a, f, GM = inqell()
    # excentricity e(squared) and semi - minor axis
    e2 = 2 * f - f ** 2
    b = (1 - f) * a
    [m, n] = xyz.shape

    if n == 3 and m == 3:
        xyz = xyz.transpose()

    r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2);

    if method == 1:
        # compute phi via iteration
        Np = xyz[:, 2]
        for i in range(0, 4):
            phi = np.arctan((xyz[:, 2] + e2 * Np) / r)
            N = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)
            Np = N * np.sin(phi)

    else:
        # compute phi using B.R.Bowring's equation (default method)
        u = np.arctan2(xyz[:, 2] * a, r * b)
        phi = np.arctan2(xyz[:, 2] + (e2 / (1 - e2) * b) * np.sin(u) ** 3, r - (e2 * a) * np.cos(u) ** 3)
        N = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)

    plh = np.array([phi, np.arctan2(xyz[:, 1], xyz[:, 0]), r / np.cos(phi) - N])
    plh = plh.transpose()

    return plh


def plh2xyz(plh, ellipse='WGS-84'):
    """

    :param plh:Nx3 matrix PLH with in the rows ellipsoidal coordinates Phi, Lambda and h
    :param ellipse: specify the ellipsoid. ELLIPS is a text is a text string with the name of the ellipsoid or a vector with the  semi-major axis a and flattening 1/f. Default for ellips is 'WGS-84'.
    :return:Nx3 matrix XYZ with cartesian coordinates X, Y and Z.
    """
    a, f, GM = inqell()
    # excentricity e(squared)
    e2 = 2 * f - f ** 2
    [m, n] = plh.shape

    if n == 3 and m == 3:
        plh = plh.transpose()

    N = a / np.sqrt(1 - e2 * np.sin(plh[:, 0]) ** 2)
    xyz = np.array([
        (N + plh[:, 2]) * np.cos(plh[:, 0]) * np.cos(plh[:, 1]),
        (N + plh[:, 2]) * np.cos(plh[:, 0]) * np.sin(plh[:, 1]),
        (N - e2 * N + plh[:, 2]) * np.sin(plh[:, 0])
    ])
    xyz = xyz.transpose()

    return xyz


# ----------------------------------------------------------------------------
#                            KEPLERIAN ELEMENTS
# ----------------------------------------------------------------------------

def vec2orb(svec, GM=3986004418e5):
    """
    VEC2ORB   Convert inertial state vector into Keplerian elements.
    ORB=VEC2ORB(SVEC) converts a 6-element inertial state vector SVEC
    with cartesian position and velocity (X, Y, Z, Xdot, Ydot, Zdot) into
    a vector ORB with 6 Keplerian elements, with
          ORB(:,1)    Semi-major axis (meters),
          ORB(:,2)    Eccentricity (unity),
          ORB(:,3)    Inclination (radians),
          ORB(:,4)    Right ascension of the ascending node (radians),
         ORB(:,5)    Argument of the pericenter (radians),
          ORB(:,6)    True anomaly (radians).
    This routine is fully vectorized, if SVEC is a matrix then ORB will be
    a matrix of the same size. One of the dimensions of the input matrix must
    be 6. Units are meter, meter/sec or radians.

    ORB=VEC2ORB(SVEC,GM) provides an optional gravitational parameter of the
    central body. Default for GM [meter**3/sec**2] is the IERS 1996 standard
    value for the Earth (GM=3986004418e5)

    [ORB,CORBTYPE]=VEC2ORB(SVEC) also outputs a character array CORBTYPE
    with the orbit type, with
        'ei'   elliptical inclined     (all Kepler elements defined)
        'ci'   circular inclined       (w =0, nu=arglat)
        'ee'   elliptical equatorial   (w=lonper, omega=0)
        'ce'   circular equatorial     (w=0, omega=0, nu=truelon)
    orbits are "circular" if the eccentricity < eps and "equatorial" if
    the inclination < eps, with eps=1e-8.
    """
    if type(svec == list):
        svec = np.array([svec])
    [m, n] = svec.shape

    if n != 6 and m == 6:
        svec = svec.transpose()

    # Inner products (rrdot = R.V , r=sqrt(R.R) , vsq = V.V )

    rrdot = svec[:, 0] * svec[:, 3] + svec[:, 1] * svec[:, 4] + svec[:, 2] * svec[:, 5]
    r = np.sqrt(svec[:, 0] * svec[:, 0] + svec[:, 1] * svec[:, 1] + svec[:, 2] * svec[:, 2])
    vsq = svec[:, 3] * svec[:, 3] + svec[:, 4] * svec[:, 4] + svec[:, 5] * svec[:, 5]

    # Angular momentum vector (H = R x V)

    hx = svec[:, 1] * svec[:, 5] - svec[:, 2] * svec[:, 4]
    hy = svec[:, 2] * svec[:, 3] - svec[:, 0] * svec[:, 5]
    hz = svec[:, 0] * svec[:, 4] - svec[:, 1] * svec[:, 3]

    hsini2 = hx * hx + hy * hy
    hsq = hsini2 + hz * hz
    h = np.sqrt(hsq)

    # Semi-major axis

    ainv = 2 / r - vsq / GM
    a = 1. / ainv

    # Eccentricity

    ome2 = hsq * ainv / GM
    ecc = np.sqrt(1.0 - ome2)
    ecc[ome2 > 1] = 0  # special handling of negative values

    # Inclination (0...pi)

    incl = np.arccos(hz / h)
    # Determine orbit type (for handling of special cases)
    small = 1e-8
    # idxecc = ecc < small
    # idxincl = np.logical_or((incl < small), (abs(incl - np.pi) < small))
    # Appendinng this in future version in python :
    # idx_ce = ( idxecc & idxincl );      % circular equatorial => w=0, omega=0, nu=truelon
    # idx_ci = ( idxecc & ~idxincl );     % circular inclined => w =0, nu=arglat
    # idx_ee = ( ~idxecc & idxincl );     % elliptical equatorial => w=lonper, omega=0
    # idx_ei = ( ~idxecc & ~idxincl );    % elliptical inclined
    #
    # orbtype(idx_ei)=0;
    # orbtype(idx_ee)=1;
    # orbtype(idx_ci)=2;
    # orbtype(idx_ce)=3;
    #
    # corbdef=['ei';'ee';'ci';'ce'];
    # corbtype=corbdef(orbtype+1,:);

    # Standard handling of elliptical inclined orbits...

    # The computations below do not do special hanling of circular or equatorial
    # orbits. This is possible because atan2(0,0)=0 is defined in Matlab, however,
    # the some of the angles will actually be undefined.

    # Longitude of ascending node (0...2*pi)
    omega = np.arctan2(hx, -hy)
    idx = (omega < 0).nonzero()
    omega[idx] = omega[idx] + 2 * np.pi

    # True anomaly (0...2*pi)

    resinf = a * ome2 * rrdot / h
    recosf = a * ome2 - r
    nu = np.arctan2(resinf, recosf)
    idx = (nu < 0).nonzero()
    nu[idx] = nu[idx] + 2 * np.pi

    # Argument of perigee (0...2*pi)

    suprod = -hz * (svec[:, 0] * hx + svec[:, 1] * hy) + svec[:, 2] * hsini2
    cuprod = h * (-svec[:, 0] * hy + svec[:, 1] * hx)
    w = np.arctan2(suprod * recosf - cuprod * resinf, cuprod * recosf + suprod * resinf)
    idx = (w < 0).nonzero()
    w[idx] = w[idx] + 2 * np.pi

    orb = np.array([a, ecc, incl, omega, w, nu])
    orb = orb.transpose()

    return orb


def orb2vec(orb, GM=3986004418e5):
    """
    ORB2VEC   Convert Keplerian elements into inertial state vector.
    SVEC=ORB2VEC(ORB) converts the vector ORB with 6 Keplerian elements
    into the 6-element inertial state vector SVEC with cartesian position and
    velocity (X, Y, Z, Xdot, Ydot, Zdot). The Keplerian elements are
          ORB(:,1)    Semi-major axis (meters),
          ORB(:,2)    Eccentricity (unity),
          ORB(:,3)    Inclination (radians),
          ORB(:,4)    Right ascension of the ascending node (radians),
          ORB(:,5)    Argument of the pericenter (radians),
          ORB(:,6)    True anomaly (radians).
    This routine is fully vectorized, if ORB is a matrix then SVEC
    will be a matrix of the same size. One of the dimensions of the input
    matrix must be 6. The units are meter, meter/sec or radians.

    SVEC=ORB2VEC(ORB,GM) provides an optional gravitational parameter of the
    central body. Default for GM [meter**3/sec**2] is the IERS 1996 standard
    value for the Earth (GM=3986004418e5)
    """
    if type(orb) == list:
        orb = np.array([orb])

    [m, n] = orb.shape

    if n != 6 and m == 6:
        orb = orb.transpose()

    # Compute position (rx,ry) and velocity (vx,vy) in orbital plane (perifocal system)

    ecc = orb[:, 1]  # Eccentricity
    cosnu = np.cos(orb[:, 5])  # Cosine and sine of true anomaly (nu)
    sinnu = np.sin(orb[:, 5])

    p = orb[:, 0] * (1.0 - ecc ** 2)  # Parameter of the ellipse p=a*(1-e^2)

    r = p / (1.0 + ecc * cosnu)  # Length of position vector

    rx = r * cosnu  # Position (rx,ry) in orbital plane
    ry = r * sinnu

    p[abs(p) < 0] = 0.0001  # Protect against division by zero
    tmp = np.sqrt(GM) / np.sqrt(p)

    vx = -tmp * sinnu  # Velocity (vx,vy) in orbital plane
    vy = tmp * (ecc + cosnu)

    # Convert into inertial frame (3-1-3 Euler rotations)

    cosincl = np.cos(orb[:, 2])  # Cosine and sine of inclination (incl)
    sinincl = np.sin(orb[:, 2])
    cosomega = np.cos(orb[:, 3])  # Cosine and sine of longitude of ascending node (omega)
    sinomega = np.sin(orb[:, 3])
    cosw = np.cos(orb[:, 4])  # Cosine and sine of argument of perigee (w)
    sinw = np.sin(orb[:, 4])

    rx0 = cosw * rx - sinw * ry  # Cosine and sine of argument of latitude u=w+nu
    ry0 = cosw * ry + sinw * rx

    vx0 = cosw * vx - sinw * vy
    vy0 = cosw * vy + sinw * vx

    svec = np.array([rx0 * cosomega - ry0 * cosincl * sinomega,
                     rx0 * sinomega + ry0 * cosincl * cosomega,
                     ry0 * sinincl,
                     vx0 * cosomega - vy0 * cosincl * sinomega,
                     vx0 * sinomega + vy0 * cosincl * cosomega,
                     vy0 * sinincl])

    svec = svec.transpose()

    return svec


def kepler(E, ecc):
    """
    KEPLER    Compute mean anomaly from eccentric anomaly (Kepler's equation).
    M=KEPLER(E,ECC) computes the mean anomaly M [rad] from the eccentric
    anomaly E [rad] and eccentricity ECC. This routine is fully vectorized,
    if E is a vector, M will be a vector with the same dimensions.

    This routine should only be used for elliptical orbits. Parabolic and
    hyperbolic orbits are not supported and give false results (this is
    nowhere checked for).
    """
    M = E - ecc * np.sin(E)

    return M


def keplernu(nu, ecc):
    """
    KEPLERNU  Compute mean anomaly from true anomaly (Kepler's equation).
    M=KEPLERNU(NU,ECC) computes the mean anomaly M [rad] from the true
    anomaly NU [rad] and eccentricity ECC. This routine is fully vectorized,
    if NU is a vector, M will be a vector with the same dimensions.

    [M,E]=KEPLERNU(NU,ECC) returns also the eccentric anomaly E [rad].

    This routine should only be used for elliptical orbits. Parabolic and
    hyperbolic orbits are not supported and give false results (this is
    nowhere checked for).
    """
    denom = 1.0 + ecc * np.cos(nu)
    sine = (np.sqrt(1.0 - ecc * ecc) * np.sin(nu)) / denom
    cose = (ecc + np.cos(nu)) / denom
    E = np.arctan2(sine, cose)
    # Compute mean anomaly
    M = E - ecc * np.sin(E)

    return M, E


def keplerm(M, ecc, TOL=1e-10):
    """
     KEPLERM  Compute eccentric/true from mean anomaly solving Kepler's eqn.
    E=KEPLERM(M,ECC) computes the eccentric anomaly E [rad] from the mean
    anomaly M [rad] and eccentricity ECC by solving Kepler's equation
    M=E-ECC*sin(E) iteratively using Newton's method. This routine is fully
    vectorized, if M is a vector, E will be a vector with the same dimensions.

    [E,NU]=KEPLERM(M,ECC) returns also the true anomaly NU [rad].

    E=KEPLERM(M,ECC,TOL) uses TOL as the cutoff criterion for the iterations,
    the default value is 1e-10.

    This routine should only be used for elliptical orbits. Parabolic and
    hyperbolic orbits are not supported and give false results (this is
    nowhere checked for).
    """
    E = M  # Use M for the first value of E
    [m, n] = E.shape
    f = np.ones((m, n))  # Newton's method for root finding
    while max(abs(f)) > TOL:
        f = M - E + ecc * np.sin(E)  # Kepler's Equation
        fdot = -1 + ecc * np.cos(E)  # Derivative of Kepler's equation
        E = E - f / fdot

    sinnu = -1 * np.sqrt(1 - ecc * ecc) * np.sin(E) / fdot
    cosnu = (ecc - np.cos(E)) / fdot

    nu = np.arctan2(sinnu, cosnu)  # True anomaly

    return E, nu
