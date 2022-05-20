# -*- coding: utf-8 -*-
"""
*geoproj* is a package with functions transform between geographic
coordinates and apply cartographic projections

Copyright (c) 2019 Centro Sismológico Nacional, Universidad de Chile

This file is part of "geoproj".

"geoproj" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"geoproj" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "geoproj". If not, see <https://www.gnu.org/licenses/>.
"""
from typing import Optional, Tuple, Sequence
import numpy as np
import math as mt
import geoproj.const as c
from geoproj.aux import trig, floatOrArray, floatOrInt


# ---------------------------------------------------------------------------
# Transverse Mercator projection

class TransverseMercator:
    """
    This class implements a Transverse Mercator projection.

    :ivar lambda0: longitude of the central meridian in radians.
    :ivar phi0: latitude of the center of the projection in radians.
    :ivar m0: auxiliary parameter.
    :ivar x0: false zero
    :ivar y0: false zero
    """
    def __init__(self, lon: float, lat: float, x0: floatOrInt = 0.0,
                 y0: floatOrInt = 0.0, rad: bool = False):
        """Initializes the object

        :param lon: longitude of the central meridian.
        :param lat: latitude of the center of the projection.
        :param x0: east coordinate false zero.
        :param y0: north coordinate false zero.
        :param rad: indicates if geographical coordinates are in radians or
            degrees.
        """
        if rad:
            self.lambda0 = lon
            self.phi0 = lat
            self.lon_ref = c.rad2deg*lon
            self.lat_ref = c.rad2deg*lat
        else:
            self.lon_ref = lon
            self.lat_ref = lat
            self.lambda0 = c.deg2rad*lon
            self.phi0 = c.deg2rad*lat
        self.m0 = _m_aux(self.phi0)
        self.x0 = x0
        self.y0 = y0

    def reset(self, lon: float, lat: float , x0: floatOrInt = 0.0,
              y0: floatOrInt = 0.0, rad: bool = False) -> None:
        """Resets the center of the TM projection.

        See :py:meth:`geo.TransverseMercator.__init__`."""
        if rad:
            self.lambda0 = lon
            self.phi0 = lat
        else:
            self.lambda0 = c.deg2rad*lon
            self.phi0 = c.deg2rad*lat
        self.m0 = _m_aux(self.phi0)
        self.x0 = x0
        self.y0 = y0

    def set_origin(self, lon_ref: float, lat_ref: float,
                   rad: bool = False) -> None:
        """Set the false zero.

        :param lon_ref: longitude of reference point.
        :param lat_ref: latitude of reference point
        :param rad: indicates if geographical coordinates are in radians or
            degrees.
        """
        e1, n1 = self(lon_ref, lat_ref, rad=rad)
        self.x0, self.y0 = -e1, -n1

    def __call__(self, lon: floatOrArray, lat: floatOrArray,
                 rad: bool = False) -> Tuple:
        """Apply transverse Mercator projection on geographical coordinates.

        :param lon: longitude.
        :param lat: latitude.
        :param rad: indicates if geographical coordinates are in radians or
            degrees.
        :return: east and north coordinates in metres.
        """
        if rad:
            lon0 = lon
            lat0 = lat
        else:
            lon0, lat0 = c.deg2rad*lon, c.deg2rad*lat
        cos_lat = np.cos(lat0)
        sin_lat = np.sin(lat0)
        aux0 = np.sqrt(1-c.e2*sin_lat*sin_lat)
        rho = c.a*(1-c.e2)/(aux0**3)
        nu = c.a/aux0
        psi = nu/rho
        t2 = (sin_lat/cos_lat)**2  # t2 = tan(lat)**2
        omega = lon0 - self.lambda0

        aux = (omega*omega)*cos_lat*cos_lat
        # aux0 = nu*omega**2*sin_lat*cos_lat
        aux0 = (nu*omega*omega)*sin_lat*cos_lat
        aux1 = aux0/2
        aux0 *= aux  # aux0 = nu*omega**4*sin_lat*cos_lat**3
        aux2 = aux0*(psi*(4*psi+1)-t2)/24
        aux0 *= aux  # aux0 = nu*omega**6*sin_lat*cos_lat**5
        aux3 = aux0*(
            psi*(psi*(psi*(8*psi*(11-24*t2) - 28*(1-6*t2)) +
                      1-32*t2) - 2*t2) + t2*t2)/720
        aux0 *= aux  # aux0 = nu*omega**8*sin_lat*cos_lat**7
        aux4 = aux0*(1385 + t2*(-3111+t2*(543-t2)))/40320

        # north coordinate
        y = self.y0 + c.ko*(
                _m_aux(lat0) - self.m0 + aux1 + aux2 + aux3 + aux4)

        aux0 = aux  # aux0 = omega**2*cos_lat**2
        aux1 = aux0*(psi - t2)/6
        aux0 *= aux  # aux0 = omega**4*cos_lat**4
        aux2 = aux0*(psi*(psi*(4*psi*(1-6*t2) + 1+8*t2) - 2*t2) + t2*t2)/120
        aux0 *= aux  # aux0 = omega**6*cos_lat**6
        aux3 = aux0*(61 + t2*(-479 + t2*(179-t2)))/5040

        # east coordinate
        x = self.x0 + c.ko*nu*omega*cos_lat*(1 + aux1 + aux2 + aux3)

        return x, y

    def inv(self, x: floatOrArray, y: floatOrArray,
            rad: bool = False) -> Tuple:
        """Inverse transverse Mercator projection.

        :param x: east coordinate in metres.
        :param y: north coordinate in metres.
        :param rad: indicates if geographical coordinates are given in radians
            or degrees.
        :return: longitude and latitude in degrees or radians.
        """
        e_aux = (x - self.x0)/c.ko
        m1 = self.m0 + (y - self.y0)/c.ko

        sigma = m1/c.g_aux
        lat1 = (sigma + c.no*(3/2.0 - c.no2*27/32.)*np.sin(2*sigma) +
                c.no2*(21/16.0 - c.no2*55/32.)*np.sin(4*sigma) +
                c.no2*(c.no*(151/96.)*np.sin(6*sigma) +
                       c.no2*(1097/512.)*np.sin(8*sigma)))
        cos_lat = np.cos(lat1)
        sin_lat = np.sin(lat1)
        aux0 = np.sqrt(1-c.e2*sin_lat*sin_lat)
        rho1 = c.a*(1-c.e2)/(aux0**3)
        nu1 = c.a/aux0
        psi1 = nu1/rho1
        t1 = sin_lat/cos_lat  # t1= tan(lat1)
        t2 = t1*t1  # t2 = tan(lat1)**2

        x1 = e_aux/nu1
        x2 = x1*x1

        aux1 = (psi1*(-4*psi1+9*(1-t2))+12*t2)/24.
        aux2 = -(360*t2*t2 + psi1*(180*t2*(5-3*t2) +
                                   psi1*(15*(15-98*t2+15*t2*t2) +
                                         psi1*(-12*(21-71*t2) +
                                               psi1*8*(11-24*t2)))))/720.
        aux3 = (1385 + t2*(3633 + t2*(4095 + t2*1575)))/40320.
        lat = lat1 + (x1*t1*e_aux/rho1)*(
            -0.5 + x2*(aux1 + x2*(aux2 + x2*aux3)))

        aux1 = -(psi1 + 2*t2)/6.
        aux2 = (24*t2*t2 + psi1*(
            72*t2 + psi1*(9-68*t2 - psi1*4*(1-6*t2))))/120.
        aux3 = -(61+t2*(662+t2*(1320+t2*720)))/5040.
        lon = self.lambda0 + (x1/cos_lat)*(
            1.0 + x2*(aux1 + x2*(aux2 + x2*aux3)))

        if rad:
            return lon, lat
        else:
            return c.rad2deg*lon, c.rad2deg*lat

    def east_north(self, lon: floatOrArray, lat: floatOrArray) -> Tuple:
        """East and north direction(s) in (x, y) coordinates

        :param lon: longitude.
        :param lat: latitude.
        :return: east direction(s), north direction(s) as returned by
        """
        # north direction
        x, y = self(lon, lat + c.delta_half)
        aux1, aux2 = self(lon, lat - c.delta_half)
        x -= aux1
        y -= aux2
        aux1 = np.sqrt(x*x + y*y)
        x /= aux1
        y /= aux1
        return (y, -x), (x, y)

    def east(self, lon: floatOrArray, lat: floatOrArray) -> Tuple:
        """East direction(s) in (x, y) coordinates

        :param lon: longitude.
        :param lat: latitude.
        :return: east direction(s): "x" and "y" components
        """
        x, y = self(lon + c.delta_half, lat)
        aux1, aux2 = self(lon - c.delta_half, lat)
        x -= aux1
        y -= aux2
        aux1 = np.sqrt(x*x + y*y)
        x /= aux1
        y /= aux1
        return x, y

    def north(self, lon: floatOrArray, lat: floatOrArray) -> Tuple:
        """North direction(s) in (x, y) coordinates

        :param lon: longitude.
        :param lat: latitude.
        :return: north direction(s): "x" and "y" components
        """
        x, y = self(lon, lat + c.delta_half)
        aux1, aux2 = self(lon, lat - c.delta_half)
        x -= aux1
        y -= aux2
        aux1 = np.sqrt(x*x + y*y)
        x /= aux1
        y /= aux1
        return x, y


def _m_aux(lat: float) -> float:
    """Auxiliary function for Transverse Mercator projection.

    :param lat: latitude in radians.
    """
    return c.a*(c.A0*lat - c.A2*np.sin(2*lat) +
                c.A4*np.sin(4*lat) - c.A6*np.sin(6*lat))
# -----------------------------------------------------------------------------


def diff_easting_northing(lat_ref: float, h_ref: Optional[float] = None,
                          deg: bool = True) -> Tuple:
    r"""
    Coefficients to compute easting and northing equivalent to small changes
    in latitude and longitude. :math:`\kappa_e, \kappa_n`.

    Using the definition given in :py:func:`coords.geo2ecef` after some algebra
    we obtain

    .. math ::
        \frac{\partial \vec{r}_{ecef}}{\partial\lambda} =
        \kappa_e \, \hat{e}

        \frac{\partial \vec{r}_{ecef}}{\partial\varphi} =
        \kappa_n \, \hat{n}

    where

    .. math ::
        \kappa_e = \left\{\frac{a}{\sqrt{1-e^2 \sin^2 \varphi}} + h\right\}
        \cos\varphi

        \kappa_n = \frac{a(1-e^2)}{(1-e^2 \sin^2 \varphi)^{\frac{3}{2}}} + h

    and :math:`\hat{e},\,\hat{n}` are the east and north unit vectors. We
    conclude that for small variations in longitude and latitude, the easting
    and northing are given by

    .. math ::
        \Delta E = \kappa_e \, \Delta\lambda

        \Delta N = \kappa_n \, \Delta\varphi


    :param lat_ref: reference latitude :math:`\lambda` in degrees.
    :param h_ref: reference height :math:`h` in metres.
    :param deg: coefficients in metres/degree or metres/radian?
    :return: :math:`\kappa_e,\,\kappa_n` in metres.
    """
    h_ref = 0.0 if h_ref is None else h_ref
    cos_phi, sen_phi = trig(lat_ref*c.deg2rad)
    aux1 = 1.0/np.sqrt(1.0-c.e2*sen_phi*sen_phi)
    if deg:
        return (c.deg2rad*((c.a*aux1 + h_ref)*cos_phi),
                c.deg2rad*(c.a*(1-c.e2)*aux1*aux1*aux1 + h_ref))
    else:
        return ((c.a*aux1 + h_ref)*cos_phi,
                c.a*(1-c.e2)*aux1*aux1*aux1 + h_ref)


def easting_northing(d_lon: floatOrArray, d_lat: floatOrArray, lat_ref: float,
                     h_ref: Optional[float] = None) -> Tuple:
    r"""
    Variations in easting and northing equivalent to small changes
    in latitude and longitude.

    .. math ::
        \Delta E = (\frac{a}{\sqrt{1-e^2 \sin^2 \lambda}} + h)\cos\lambda
            \, \Delta\varphi[rad]

        \Delta N = (\frac{a(1-e^2)}{(1-e^2 \sin^2 \lambda)^{\frac{3}{2}}} + h)
            \, \Delta\lambda[rad]


    :param lat_ref: latitude :math:`\lambda` in degrees.
    :param h_ref: height :math:`h` in metres.
    :param d_lon: variation in longitude :math:`\Delta\varphi` in degrees.
    :param d_lat: variation in latitude :math:`\Delta\lambda` in degrees.
    :return: :math:`\Delta E, \Delta N` in metres.
    """
    kappa_e, kappa_n = diff_easting_northing(lat_ref, h_ref=h_ref, deg=True)
    return kappa_e*d_lon, kappa_n*d_lat


def delta_lon_lat(d_east: floatOrArray, d_north: floatOrArray,
                  lat_ref: float, h_ref: Optional[float] = None) -> Tuple:
    r"""
    Latitude and longitude variations due to small changes in
    easting and northing.

    .. math ::
        \Delta\varphi[rad] =
            \frac{\Delta E}{(\frac{a}{\sqrt{1-e^2 \sin^2 \lambda}} + h)
            \cos\lambda}

        \Delta\lambda[rad] =
            \frac{\Delta N}{\frac{a(1-e^2)}{(1-e^2 \sin^2 \lambda)
            ^{\frac{3}{2}}} + h}

    :param d_east: variation in easting :math:`\Delta E` in metres.
    :param d_north: variation in northing :math:`\Delta N` in metres.
    :param lat_ref: latitude :math:`\lambda` in degrees.
    :param h_ref: height :math:`h` in metres.
    :return: :math:`\Delta \varphi, \Delta \lambda` in degrees.
    """
    h_ref = 0.0 if h_ref is None else h_ref
    cos_psi, sen_psi = trig(lat_ref*c.deg2rad)
    aux1 = 1.0/np.sqrt(1.0-c.e2*sen_psi*sen_psi)
    return (c.rad2deg*d_east/((c.a*aux1 + h_ref)*cos_psi),
            c.rad2deg*d_north/(c.a*(1-c.e2)*aux1*aux1*aux1 + h_ref))


def distance_approx(lon_lat_1: Sequence, lon_lat_2: Sequence,
                    rough: bool = True) -> float:
    r"""Approximate distance between two near points on the earth surface.

    .. math::
        d = a\sqrt{(\lambda_1-\lambda_2)^2 +
        \cos^2(\tfrac{1}{2}(\lambda_1+\lambda_2))(\varphi_1-\varphi_2)^2}

        \text{where all angles are in radians}

    :param lon_lat_1: longitude and latitude :math:`\varphi_1, \lambda_1` of
        the first point in degrees.
    :param lon_lat_2: longitude and latitude :math:`\varphi_1, \lambda_1` of
        the second point in degrees.
    :param rough: rough approximation?
    :return: distance :math:`d` in metres.

    .. warning:: The points must be close, less than 100 km.
    """
    if rough:
        phi_1 = c.deg2rad*lon_lat_1[1]
        phi_2 = c.deg2rad*lon_lat_2[1]
        cos_lat = mt.cos(0.5*(phi_1 + phi_2))
        return c.a*mt.sqrt((cos_lat*c.deg2rad*(lon_lat_2[0]-lon_lat_1[0]))**2
                           + (phi_2-phi_1)**2)
    de, dn = easting_northing(lon_lat_2[0] - lon_lat_1[0],
                              lon_lat_2[1] - lon_lat_1[1],
                              0.5*(lon_lat_1[1] + lon_lat_2[1]))
    return mt.sqrt(de*de + dn*dn)


# ---------------------------------------------------------------------------
# "Web Mercator" projection, a.k.a. epsg:3857

def geo2web_merc(lon: floatOrArray, lat:floatOrArray,
                 rad: bool = False) -> Tuple:
    r"""From longitude-latitude to web Mercator projection (epsg:3857).

    Given the longitude :math:`\lambda` and latitude :math:`\varphi`,
    the easting :math:`x` and northing :math:`y` are

    .. math::
        x = a\,\lambda

        y = a\,\log(\tan(\frac{\pi}{4} + \frac{\varphi}{2}))

    :param lon: longitudes(s) in degrees or radians
    :param lat: latitude(s) in degrees or radians
    :param rad: angles in radians?
    :return: easting(s), northing(s)

    .. note:: Both pair of coordinates are defined in the WGS84 ellipsoid
    """
    if rad:
        ln = lon
        lt = lat
    else:
        ln = c.deg2rad*lon
        lt = c.deg2rad*lat
    return c.a*ln, c.a*np.log(np.tan(0.5*lt + 0.785398163397448270))


def web_merc2geo(x: floatOrArray, y: floatOrArray, rad: bool = False) -> Tuple:
    r"""From web Mercator projection (epsg:3857) to longitude-latitude.

    Given the easting :math:`x` and northing :math:`y`,
    the longitude :math:`\lambda` and latitude :math:`\varphi` are

    .. math::
        \lambda = \frac{x}{a}

        \varphi = 2\,\arctan\left(e^{\frac{y}{a}}\right) - \frac{\pi}{2}

    :param x: easting(s) in metres
    :param y: northing(s) in metres
    :param rad: return angles in radians?
    :return: longitude(s), latitude(s)

    .. note:: Both pair of coordinates are defined in the WGS84 ellipsoid
    """
    lon = x/c.a
    lat = 2*np.arctan(np.exp(y/c.a)) - 1.5707963267948966
    if not rad:
        lon *= c.rad2deg
        lat *= c.rad2deg
    return lon, lat
