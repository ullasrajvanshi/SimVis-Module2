# TLE (Satellite NORAD Two Line Elements) Toolbox
# Version 1.2-l (September 13, 2017) - Limited edition for CIE4604
#
# Get and read NORAD Two Line Elements in TLE structure array
#   tleget      - Retrieve NORAD Two Line Elements from www.celestrak.com
#   tleread     - Read NORAD Two Line Elements from file.
#
# Plotting of satellite positions, elevation, azimuth, visibility, rise/set, ...
#   tleplot1    - Plot satellite position and velocity from NORAD Two Line Elements.
#
# Find satellites and select dates
#   tlefind     - Find named satellites in the NORAD Two Line Elements.
#   tledatenum  - Compute Matlab datenumbers from a date range.
#
# Compute satellite positions and orbit propagation
#   tle2vec1    - Satellite position and velocity from NORAD Two Line Elements.
#   tle2orb     - Compute orbital elements from NORAD Two Line Elements
#
# Examples:
#   tle=tleread('resource-10-oct-2017.tle');
#   tlefind(tle,'SENTINEL');
#   tleplot1(tle,{'2017-10-10 0:00', 24*60 ,1},'SENTINEL-1A',[ 52 4.8  0 ]);
#
#
# (c) Hans van der Marel & Ullas Rajvanshi Delft University of Technology, 2012-2017.
# !/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011 - 2018

# Author(s):

#   Esben S. Nielsen <esn@dmi.dk>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Adding reading TLE files
# !/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011 - 2018

# Author(s):

#   Esben S. Nielsen <esn@dmi.dk>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import io
import logging
import datetime
import crsutil as crs

try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
import os
import glob
import numpy as np

TLE_URLS = ('http://celestrak.com/NORAD/elements/weather.txt',
            'http://celestrak.com/NORAD/elements/resource.txt',
            'https://www.celestrak.com/NORAD/elements/cubesat.txt',
            'http://celestrak.com/NORAD/elements/stations.txt',
            'https://www.celestrak.com/NORAD/elements/sarsat.txt',
            'https://www.celestrak.com/NORAD/elements/noaa.txt',
            'https://www.celestrak.com/NORAD/elements/amateur.txt')

LOGGER = logging.getLogger(__name__)
PKG_CONFIG_DIR = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'etc')


def read_platform_numbers(in_upper=False, num_as_int=False):
    """Read platform numbers from $PPP_CONFIG_DIR/platforms.txt if available."""
    out_dict = {}
    os.getenv('PPP_CONFIG_DIR', PKG_CONFIG_DIR)
    platform_file = None
    if 'PPP_CONFIG_DIR' in os.environ:
        platform_file = os.path.join(os.environ['PPP_CONFIG_DIR'], 'platforms.txt')
    if not platform_file or not os.path.isfile(platform_file):
        platform_file = os.path.join(PKG_CONFIG_DIR, 'platforms.txt')

    try:
        fid = open(platform_file, 'r')
    except IOError:
        # LOGGER.error("Platform file %s not found.", platform_file)
        return out_dict
    for row in fid:
        # skip comment lines
        if not row.startswith('#'):
            parts = row.split()
            if len(parts) < 2:
                continue
            if in_upper:
                parts[0] = parts[0].upper()
            if num_as_int:
                parts[1] = int(parts[1])
            out_dict[parts[0]] = parts[1]
    fid.close()

    return out_dict


SATELLITES = read_platform_numbers(in_upper=True, num_as_int=False)
'''
The platform numbers are given in a file $PPP_CONFIG/platforms.txt
in the following format:

.. literalinclude:: ../../etc/platforms.txt
  :language: text
  :lines: 4-
'''


def read(platform, tle_file=None, line1=None, line2=None):
    """Read TLE for `platform` from `tle_file`

    File is read from `line1` to `line2`, from the newest file provided in the
    TLES pattern, or from internet if none is provided.
    """
    return Tle(platform, tle_file=tle_file, line1=line1, line2=line2)


def fetch(destination):
    """Fetch TLE from internet and save it to `destination`."""
    with io.open(destination, mode="w", encoding="utf-8") as dest:
        for url in TLE_URLS:
            response = urlopen(url)
            dest.write(response.read().decode("utf-8"))


class ChecksumError(Exception):
    """ChecksumError."""
    pass


class Tle(object):
    """Class holding TLE objects."""

    def __init__(self, platform, tle_file=None, line1=None, line2=None):
        self._platform = platform.strip().upper()
        self._tle_file = tle_file
        self._line1 = line1
        self._line2 = line2

        self.satnumber = None
        self.classification = None
        self.id_launch_year = None
        self.id_launch_number = None
        self.id_launch_piece = None
        self.epoch_year = None
        self.epoch_day = None
        self.epoch = None
        self.mean_motion_derivative = None
        self.mean_motion_sec_derivative = None
        self.bstar = None
        self.ephemeris_type = None
        self.element_number = None
        self.inclination = None
        self.right_ascension = None
        self.excentricity = None
        self.arg_perigee = None
        self.mean_anomaly = None
        self.mean_motion = None
        self.orbit = None

        self._read_tle()
        self._checksum()
        self._parse_tle()

    @property
    def line1(self):
        """Return first TLE line."""
        return self._line1

    @property
    def line2(self):
        """Return second TLE line."""
        return self._line2

    @property
    def platform(self):
        """Return satellite platform name."""
        return self._platform

    def _checksum(self):
        """Performs the checksum for the current TLE."""
        for line in [self._line1, self._line2]:
            check = 0
            for char in line[:-1]:
                if char.isdigit():
                    check += int(char)
                if char == "-":
                    check += 1

            if (check % 10) != int(line[-1]):
                raise ChecksumError(self._platform + " " + line)

    def _read_tle(self):
        """Read TLE data."""
        if self._line1 is not None and self._line2 is not None:
            tle = self._line1.strip() + "\n" + self._line2.strip()
        else:
            def _open(filename):
                return io.open(filename, 'rb')

            if self._tle_file:
                urls = (self._tle_file,)
                open_func = _open
            elif "TLES" in os.environ:
                # TODO: get the TLE file closest in time to the actual satellite
                # overpass, NOT the latest!
                urls = (max(glob.glob(os.environ["TLES"]),
                            key=os.path.getctime),)
                LOGGER.debug("Reading TLE from %s", urls[0])
                open_func = _open
            else:
                LOGGER.debug("Fetch TLE from the internet.")
                urls = TLE_URLS
                open_func = urlopen

            tle = ""
            designator = "1 " + SATELLITES.get(self._platform, '')
            for url in urls:
                fid = open_func(url)
                for l_0 in fid:
                    l_0 = l_0.decode('utf-8')
                    if l_0.strip() == self._platform:
                        l_1 = next(fid).decode('utf-8')
                        l_2 = next(fid).decode('utf-8')
                        tle = l_1.strip() + "\n" + l_2.strip()
                        break
                    if (self._platform in SATELLITES and
                            l_0.strip().startswith(designator)):
                        l_1 = l_0
                        l_2 = next(fid).decode('utf-8')
                        tle = l_1.strip() + "\n" + l_2.strip()
                        LOGGER.debug("Found platform %s, ID: %s",
                                     self._platform,
                                     SATELLITES[self._platform])
                        break
                fid.close()
                if tle:
                    break

            if not tle:
                raise KeyError("Found no TLE entry for '%s'" % self._platform)

        self._line1, self._line2 = tle.split('\n')

    def _parse_tle(self):
        """Parse values from TLE data."""

        def _read_tle_decimal(rep):
            """Convert *rep* to decimal value."""
            if rep[0] in ["-", " ", "+"]:
                digits = rep[1:-2].strip()
                val = rep[0] + "." + digits + "e" + rep[-2:]
            else:
                digits = rep[:-2].strip()
                val = "." + digits + "e" + rep[-2:]

            return float(val)

        self.satnumber = self._line1[2:7]
        self.classification = self._line1[7]
        self.id_launch_year = self._line1[9:11]
        self.id_launch_number = self._line1[11:14]
        self.id_launch_piece = self._line1[14:17]
        self.epoch_year = self._line1[18:20]
        self.epoch_day = float(self._line1[20:32])
        self.epoch = \
            np.datetime64(datetime.datetime.strptime(self.epoch_year, "%y") +
                          datetime.timedelta(days=self.epoch_day - 1), 'us')
        self.mean_motion_derivative = float(self._line1[33:43])
        self.mean_motion_sec_derivative = _read_tle_decimal(self._line1[44:52])
        self.bstar = _read_tle_decimal(self._line1[53:61])
        try:
            self.ephemeris_type = int(self._line1[62])
        except ValueError:
            self.ephemeris_type = 0
        self.element_number = int(self._line1[64:68])

        self.inclination = np.deg2rad(float(self._line2[8:16]))
        self.right_ascension = np.deg2rad(float(self._line2[17:25]))
        self.excentricity = int(self._line2[26:33]) * 10 ** -7
        self.arg_perigee = np.deg2rad(float(self._line2[34:42]))
        self.mean_anomaly = np.deg2rad(float(self._line2[43:51]))
        self.mean_motion = 2 * np.pi * (float(self._line2[52:63]))
        self.orbit = int(self._line2[63:68])
        mu = 398600.5  # Earth gravitational parameter(WGS84)[km ^ 3 / s ^ 2]
        self.a0 = (mu / (self.mean_motion / (24 * 3600)) ** 2) ** (1 / 3) * 1000  # Semi-major axis [km]

    def __str__(self):
        import pprint
        import sys
        if sys.version_info < (3, 0):
            from StringIO import StringIO
        else:
            from io import StringIO
        s_var = StringIO()
        d_var = dict(([(k, v) for k, v in
                       list(self.__dict__.items()) if k[0] != '_']))
        pprint.pprint(d_var, s_var)
        return s_var.getvalue()[:-1]


def tledatenum(y1, m1, d1, y2, m2, d2, interval=1):
    """"
    .. write some
    """

    start = datetime.date.toordinal(datetime.date(y1, m1, d1)) + 366
    end = datetime.date.toordinal(datetime.date(y2, m2, d2)) + 366
    # program
    step = interval / (24 * 60)
    t = np.arange(start, end, step)
    t = np.array([t])
    t = t.transpose()

    return t


def tle2orb(tle, t):
    """
    Using J2 propogation
    """
    J2 = 0.00108262998905  # J2
    Re = 6378136  # [m]   radius of the Earth
    mu = 3986004418e5

    # Initialize orbit propagation

    nepoch = t.shape[0]
    epoch_tle = tle.epoch_day
    doy = np.floor(epoch_tle)
    # t0 = datetime.date.toordinal(datetime.date(2000 + int(tle.epoch_year), doy)) + 366
    t0 = datetime.date.toordinal(
        datetime.date(2000 + int(tle.epoch_year), 1, 1) + datetime.timedelta(doy)) + 365 + epoch_tle - doy

    # Orbit propogation

    # Compute rate of change of orbital elements
    # draan/dt = s.*cos(inclination)
    # dargp/dt = -0.5*s.*(5*cos(inclination-1).^2)
    # dM/dt = -0.5*s.*sqrt(1-e.^2).*(3*cos(inclination).^2 -1)

    # with s=-J2*3/2*sqrt(mu/a^3)*(Re/p)^2

    # dM/dt is not needed for two line element propagation, but computed nevertheless.

    p = tle.a0 * (1 - tle.excentricity ** 2)
    s = -J2 * 3 / 2 * np.sqrt(mu / tle.a0 ** 3) * (Re / p) ** 2  # identical to s = -J2 * 3 / 2 * n0 * (Re / p) ^ 2;
    odot = s * np.cos(tle.inclination) * 86400
    wdot = -0.5 * s * (5 * np.cos(tle.inclination) ** 2 - 1) * 86400
    mdot = -0.5 * s * np.sqrt(1 - tle.excentricity ** 2) * (3 * np.cos(tle.inclination) ** 2 - 1) * 86400
    raan = tle.right_ascension + odot * (t - t0)  # Right ascension of ascending node[rad]
    argp = tle.arg_perigee + wdot * (t - t0)  # Argument of periapsis[rad]
    m = tle.mean_anomaly + tle.mean_motion * (t - t0)  # Mean anomaly[rad]
    E, nu = crs.keplerm(m + argp, tle.excentricity)  # Solve Kepler's equation for E + argp
    nu = nu - argp  # True anomaly  [rad]

    orb0 = np.tile([tle.a0, tle.excentricity, tle.inclination], [nepoch, 1])
    orb = np.concatenate((orb0, raan, argp, nu), axis=1)

    return orb


def tle2vec(tle, t):
    """
    Using J2 propogation
    """
    orb = tle2orb(tle, t)
    vec = crs.orb2vec(orb)

    xsat = np.array([vec[:, 0], vec[:, 1], vec[:, 2]])
    xsat = xsat.transpose()
    vsat = np.array([vec[:, 3], vec[:, 4], vec[:, 5]])
    vsat = vsat.transpose()

    return xsat, vsat
