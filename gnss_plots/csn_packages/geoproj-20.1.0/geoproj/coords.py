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
from typing import Union, Optional, Tuple, Sequence
import numpy as np
import math as mt
from geoproj.proj import diff_easting_northing
from geoproj.aux import trig, floatOrArray
import geoproj.const as c
TupleOrArray = Union[Tuple, np.ndarray]
NumberOrArray = Union[float, np.ndarray]


def geo2ecef(lon: NumberOrArray, lat: NumberOrArray,
             h: Optional[NumberOrArray] = None, deg: bool = True) -> Tuple:
    r"""
    From geodetic coordinates to Earth-centered earth-fixed Coordinates (ECEF).

    .. math::
        \vec{r}_{ecef}(\varphi, \lambda, h) =
        \begin{bmatrix} x \\ y \\ z
        \end{bmatrix} =
        \begin{bmatrix}
        (N + h)\cos\varphi\cos\lambda \\
        (N + h)\cos\varphi\sin\lambda \\
        (N(1 - e^2) + h)\sin\varphi
        \end{bmatrix}

        \text{where } N(\varphi) \equiv \frac{a}{\sqrt{1 - e^2\sin^2\varphi}}

        \text{and $a$, $e$ are the semi-major axis and eccentricity of the
            WGS-84 reference ellipsoid.}

    :param lon: longitude(s) :math:`\lambda`.
    :param lat: latitude(s) :math:`\varphi`.
    :param h: ellipsoidal height(s) in metres :math:`h`. If None, zero
        is assumed.
    :param deg: degrees or radians?
    :return: ECEF coordinates :math:`\bm{x},\, \bm{y},\, \bm{z}`.
    """
    cos_lon, sin_lon, cos_lat, sin_lat = _trig_lon_lat(lon, lat, deg)
    n_aux = c.a/np.sqrt(1.0-c.e2*sin_lat*sin_lat)
    if h is None:
        h1 = (0.0 if isinstance(lon, float) else np.zeros(lon.size))
    else:
        h1 = h
    x = (n_aux + h1)*cos_lat*cos_lon
    y = (n_aux + h1)*cos_lat*sin_lon
    z = (n_aux*(1.0-c.e2) + h1)*sin_lat
    return x, y, z


def ecef2geo(x: NumberOrArray, y: NumberOrArray, z: NumberOrArray) -> Tuple:
    r"""
    From Earth-centered earth-fixed coordinates (ECEF) to geodetic coordinates.

    This is the inverse function of :py:func:`coords.geo2ecef`.

    :param x: x ECEF coordinate(s).
    :param y: y ECEF coordinate(s).
    :param z: z ECEF coordinate(s).
    :return: longitude(s) in degrees, latitude(s) in degrees, ellipsoidal
        height(s) in metres.
    """
    if isinstance(x, float):
        aux = np.sqrt(x*x + y*y)
        alpha = aux/c.a
        beta = z/aux
        # mu = tan(phi)
        mu = beta
        aux1 = beta - mu*(1-c.e2/(alpha*np.sqrt(1+mu*mu*(1.0-c.e2))))
        while np.max(np.abs(aux1)) > 1.0e-12:
            mu2 = mu*mu
            aux2 = np.sqrt(1+mu2*(1.0-c.e2))
            aux1 = beta - mu*(1-c.e2/(alpha*aux2))
            mu += aux1/(1-c.e2/(alpha*aux2*aux2*aux2))
        mu2 = mu*mu
        lon, lat = c.rad2deg*np.arctan2(y, x) % 360., c.rad2deg*np.arctan(mu)
        return lon, lat, np.sqrt(1+mu2)*(aux - c.a/np.sqrt(1+mu2*(1-c.e2)))
    else:
        aux = tuple(ecef2geo(x1, x2, x3) for x1, x2, x3 in zip(x, y, z))
        return (np.array(tuple(point[0] for point in aux)),
                np.array(tuple(point[1] for point in aux)),
                np.array(tuple(point[2] for point in aux)))


def ecef_aux1(cos_lat: float, sin_lat: float, h: float) -> Tuple:
    r"""Computes :math:`m_\varphi`, :math:`m_\lambda` and :math:`N(\lambda)`.

    .. math::
        m_\varphi &\equiv
        |\frac{\partial \vec{r}_{ecef}(\varphi, \lambda, h)}
        {\partial \varphi}|, \\ m_\lambda &\equiv
        |\frac{\partial \vec{r}_{ecef}(\varphi, \lambda, h)}
        {\partial \lambda}|, \\
        N(\lambda) &\equiv \frac{a}{\sqrt{1 - e^2\sin^2\lambda}}

    :param cos_lat: cosine of latitude :math:`\cos\varphi`.
    :param sin_lat: sine of latitude :math:`\sin\varphi`.
    :param h: altitude :math:`h` in metres.
    :return: :math:`m_\varphi`, :math:`m_\lambda`, :math:`N(\lambda)`.
    """
    aux = 1.0/(1.0 - c.e2*sin_lat*sin_lat)
    n_aux = c.a*mt.sqrt(aux)
    return (n_aux + h)*cos_lat, n_aux*aux*(1.0 - c.e2) + h, n_aux


def ecef_aux2(cos_lon: float, sin_lon: float, cos_lat: float, sin_lat: float,
              h: float) -> Tuple:
    r"""
    Computes the position vector in ECEF
    :math:`\vec{r}_{ecef}(\varphi, \lambda, h)`,
    :math:`m_\varphi` and :math:`m_\lambda`, from the trigonometric functions
    of geographical coordinates. See :py:func:`coords.ecef_aux1` and
    :py:func:`coords.geo2ecef` functions for more details.

    :param cos_lon: cosine of longitude :math:`\cos\lambda`.
    :param sin_lon: sine of longitude :math:`\sin\lambda`.
    :param cos_lat: cosine of latitude :math:`\cos\varphi`.
    :param sin_lat: sine of latitude.
    :param h: altitude :math:`h` in metres.
    :return: position vector :math:`\vec{r}_{ecef}(\varphi, \lambda, h)`,
        :math:`m_\varphi`, :math:`m_\lambda`.
    """
    m_lon, m_lat, n_aux = ecef_aux1(cos_lat, sin_lat, h)
    return (np.array([(n_aux + h)*cos_lat*cos_lon,
                      (n_aux + h)*cos_lat*sin_lon,
                      (n_aux*(1.0 - c.e2) + h)*sin_lat]),
            m_lon,
            m_lat)


def ecef2enu_matrix(lon: float, lat: float, deg: bool = True) -> np.ndarray:
    r"""
    Rotation matrix from ECEF to local coordinates (east, north, up).

    Rotation matrix from Earth-centered earth-fixed Coordinates to
    (east, north, up). As input, uses longitude and latitude of the point
    where the transformation wants to be applied.

    .. math::
        R_{ecef2enu}(\lambda, \varphi) = \begin{bmatrix}
        -\sin\lambda & \cos\lambda & 0 \\
        -\sin\varphi\cos\lambda & -\sin\varphi\sin\lambda & \cos\varphi \\
        \cos\varphi\cos\lambda & \cos\varphi\sin\lambda & \sin\varphi
        \end{bmatrix}

        \vec{r}_{enu} = R_{ecef2enu}(\lambda, \varphi) \vec{r}_{ecef}

    :param lon: longitude :math:`\lambda`.
    :param lat: latitude :math:`\varphi`
    :param deg: degrees or radians?
    :return: rotation matrix :math:`R_{ecef2enu}(\varphi, \lambda)`.
    """
    return ecef2enu_matrix_trig(*_trig_lon_lat(lon, lat, deg))


def ecef2enu_matrix_trig(cos_lon: float, sin_lon: float,
                         cos_lat: float, sin_lat: float) -> np.ndarray:
    """See :py:func:`geoproj.coords.ecef2enu_matrix`"""
    e_, n_, v_ = _enu_trig(cos_lon, sin_lon, cos_lat, sin_lat, as_arrays=False)
    return np.append((e_, n_), (v_,), axis=0)


def east_north_up(lon: float, lat: float, deg: bool = True,
                  as_arrays: bool = True) -> Tuple[TupleOrArray]:
    r"""East, north and vertical directions in ECEF coordinates.

    :param lon: longitude :math:`\lambda`.
    :param lat: latitude :math:`\varphi`.
    :param deg: degrees or radians?
    :param as_arrays: return as numpy arrays or tuples?
    :return: unit east :math:`\hat{e}`, north :math:`\hat{n}` and up
        :math:`\hat{z}` directions. See :py:func:`coords._east_trig`,
        :py:func:`coords._north_trig` and :py:func:`coords._up_trig` functions.
    """
    return _enu_trig(*_trig_lon_lat(lon, lat, deg), as_arrays=as_arrays)


def east(lon: float, deg: bool = True, as_array: bool = True) -> TupleOrArray:
    r"""East direction in ECEF coordinates.

    .. math::
        \hat{e} = \begin{bmatrix} -\sin\lambda \\ \cos\lambda \\ 0
        \end{bmatrix}

    :param lon: longitude :math:`\cos\lambda`.
    :param deg: degrees or radians?
    :param as_array: return as numpy array or tuple?
    :return: unit east direction :math:`\hat{e}`.
    """
    return _east_trig(*trig(c.deg2rad*lon if deg else lon), as_array=as_array)


def north(lon: float, lat: float, deg: bool = True,
          as_array: bool = True) -> TupleOrArray:
    r"""North direction in ECEF coordinates.

    .. math::
        \hat{n} = \begin{bmatrix}
        -\sin\varphi\cos\lambda \\ -\sin\varphi\sin\lambda \\ \cos\varphi
        \end{bmatrix}

    :param lon: longitude :math:`\cos\lambda`.
    :param lat: latitude :math:`\cos\varphi`.
    :param deg: degrees or radians?
    :param as_array: return as numpy array or tuple?
    :return: unit north direction :math:`\hat{n}`.
    """
    return _north_trig(*_trig_lon_lat(lon, lat, deg), as_array=as_array)


def up(lon: float, lat: float, deg: bool = True,
       as_array: bool = True) -> TupleOrArray:
    r"""Upwards vertical direction in ECEF coordinates.

    .. math::
        \hat{z} = \begin{bmatrix}
        \cos\varphi\cos\lambda \\ \cos\varphi\sin\lambda \\ \sin\varphi
        \end{bmatrix}

    :param lon: longitude :math:`\cos\lambda`.
    :param lat: latitude :math:`\cos\varphi`.
    :param deg: degrees or radians?
    :param as_array: return as numpy array or tuple?
    :return: unit up direction :math:`\hat{z}`.
    """
    return _up_trig(*_trig_lon_lat(lon, lat, deg), as_array=as_array)


def _enu_trig(cos_lon: float, sin_lon: float, cos_lat: float, sin_lat: float,
              as_arrays: bool = True) -> Tuple:
    r"""East, north and up directions in ECEF coordinates.

    :param cos_lon: cosine of longitude :math:`\cos\lambda`.
    :param sin_lon: sine of longitude :math:`\sin\lambda`.
    :param cos_lat: cosine of latitude :math:`\cos\varphi`.
    :param sin_lat: sine of latitude :math:`\sin\varphi`.
    :param as_arrays: return as numpy arrays or tuples?
    :return: unit east :math:`\hat{e}`, north :math:`\hat{n}` and up
        :math:`\hat{z}` directions. See :py:func:`geo.east_tr`,
        :py:func:`geo.north_tr` and :py:func:`geo.up_tr` functions.
    """
    return (_east_trig(cos_lon, sin_lon,
                       as_array=as_arrays),
            _north_trig(cos_lon, sin_lon, cos_lat, sin_lat,
                        as_array=as_arrays),
            _up_trig(cos_lon, sin_lon, cos_lat, sin_lat,
                     as_array=as_arrays))


def _east_trig(cos_lon: float, sin_lon: float,
               as_array: bool = True) -> TupleOrArray:
    """See :py:func:`geoproj.coords.east`"""
    aux = (-sin_lon, cos_lon, 0.0)
    return np.array(aux) if as_array else aux


def _north_trig(cos_lon: float, sin_lon: float, cos_lat: float, sin_lat: float,
                as_array: bool = True) -> TupleOrArray:
    """See :py:func:`geoproj.coords.north`"""
    aux = (-sin_lat*cos_lon, -sin_lat*sin_lon, cos_lat)
    return np.array(aux) if as_array else aux


def _up_trig(cos_lon: float, sin_lon: float, cos_lat: float, sin_lat: float,
             as_array: bool = True) -> TupleOrArray:
    """See :py:func:`geoproj.coords.up`"""
    aux = (cos_lat*cos_lon, cos_lat*sin_lon, sin_lat)
    return np.array(aux) if as_array else aux


def _trig_lon_lat(lon: floatOrArray, lat: floatOrArray,
                  deg: bool = True) -> Tuple:
    cos_lon, sin_lon = trig(c.deg2rad*lon if deg else lon)
    cos_lat, sin_lat = trig(c.deg2rad*lat if deg else lat)
    return cos_lon, sin_lon, cos_lat, sin_lat


def wgs84_radius_prime_vertical(lat: float, deg: bool = True) -> float:
    r"""Computes the radius of curvature in the prime vertical

    .. math::
        \text{where } N(\varphi) \equiv \frac{a}{\sqrt{1 - e^2\sin^2\varphi}}

    :param lat: latitude :math:`\varphi`
    :param deg: degrees or radians?
    :return: radius of curvature :math:`N`
    """
    sin_lat = mt.sin(c.deg2rad*lat if deg else lat)
    return c.a/np.sqrt(1.0 - c.e2*sin_lat*sin_lat)


def wgs84_radius_meridian(lat: float, deg: bool = True) -> float:
    r"""Computes the radius of curvature in the meridian for the WGS-84
    reference ellipsoid

    .. math::
        \text{where } M(\varphi) \equiv
            \frac{a(1 - e^2)}{(\sqrt{1 - e^2\sin^2\varphi})^{\frac{3}{2}}}

    :param lat: latitude :math:`\varphi`
    :param deg: degrees or radians?
    :return: radius of curvature :math:`M`
    """
    sin_lat = mt.sin(c.deg2rad*lat if deg else lat)
    return c.a*(1. - c.e2)/mt.pow(1.0 - c.e2*sin_lat*sin_lat, 1.5)


def wgs84_radii(lat: float, deg: bool = True) -> Tuple:
    r"""Computes the radii of curvature in the prime vertical and in the
    meridian for the WGS-84 reference ellipsoid

    .. math::
        \text{where } N(\varphi) \equiv \frac{a}{\sqrt{1 - e^2\sin^2\varphi}}

        \text{where } M(\varphi) \equiv
            \frac{a(1 - e^2)}{(\sqrt{1 - e^2\sin^2\varphi})^{\frac{3}{2}}}

    :param lat: latitude :math:`\varphi`
    :param deg: degrees or radians?
    :return: radii of curvature :math:`N`, :math:`M`
    """
    sin_lat = mt.sin(c.deg2rad*lat if deg else lat)
    aux = 1.0 - c.e2*sin_lat*sin_lat
    n = c.a/np.sqrt(aux)
    return n, n*(1. - c.e2)/aux


def wgs84_gaussian_radius(lat: float, deg: bool = True) -> float:
    r"""Computes the gaussian radius of curvature for the WGS-84
    reference ellipsoid

    .. math::
        \text{where } M(\varphi) \equiv
            \frac{a(1 - e^2)}{(\sqrt{1 - e^2\sin^2\varphi})^{\frac{3}{2}}}

    :param lat: latitude :math:`\varphi`
    :param deg: degrees or radians?
    :return: radius of curvature :math:`M`
    """
    sin_lat = mt.sin(c.deg2rad*lat if deg else lat)
    return c.b/(1.0 - c.e2*sin_lat*sin_lat)


# -----------------------------------------------------------------------------
# geographic distance


class GeoDistance:
    """
    Computation of geographical distances in the WGS84 reference ellipsoid
    and through a spherical approximation
    """
    def __init__(self, lon_ref: float, lat_ref: float, deg: bool = True):
        """Initializes instance

        :param lon_ref: longitude of the reference point
        :param lat_ref: latitude of the reference point
        :param deg: degrees or radians?
        """
        self._lon_ref_rad = None
        self._lat_ref_rad = None
        self.cos_lon = None
        self.sin_lon = None
        self.cos_lat = None
        self.sin_lat = None
        self.r_gaussian = None
        self.set_reference_point(lon_ref, lat_ref, deg=deg)

    def great_circle(self, lon: float, lat: float, deg: bool = True) -> float:
        r"""Great circle distance between two points on a spherical earth

        :param lon: longitude of the second point :math:`\lambda_2`
        :param lat: latitude of the second point :math:`\varphi_2`
        :param deg: degrees or radians?
        :return: distance in metres
        """
        lat2_rad = (c.deg2rad * lat if deg else lat)
        r_gaussian = wgs84_gaussian_radius(0.5*(self._lat_ref_rad+lat2_rad),
                                           deg=False)
        return r_gaussian*self.angular(lon, lat, deg=deg, deg_out=False)

    def angular(self, lon: float, lat: float, deg: bool = True,
                deg_out: bool = True) -> float:
        r"""Great circle angular distance between two points on a spherical
        earth

        :param lon: longitude of the second point :math:`\lambda_2`
        :param lat: latitude of the second point :math:`\varphi_2`
        :param deg: degrees or radians?
        :param deg_out: distance degrees or radians?
        :return: angular distance
        """
        aux = mt.acos(self._dot(lon, lat, deg))
        if deg_out:
            return c.rad2deg*aux
        return aux

    def straight(self, lon: float, lat: float, height: float,
                 deg: bool = True) -> float:
        r"""Straight line distance between two points (approximate)

        :param lon: longitude of the second point :math:`\lambda_2`
        :param lat: latitude of the second point :math:`\varphi_2`
        :param height: height of the second point in metres
        :param deg: degrees or radians?
        :return: distance
        """
        cos_theta = self._dot(lon, lat, deg)
        aux = self.r_gaussian + height
        return mt.sqrt(self.r_gaussian*(self.r_gaussian - 2*aux*cos_theta) +
                       aux*aux)

    def vincenty(self, lon: float, lat: float, deg: bool = True) -> float:
        r"""Vincenty's formula for the distance between two points on the
        reference ellipsoid.

        :param lon: longitude of the second point :math:`\lambda_2`
        :param lat: latitude of the second point :math:`\varphi_2`
        :param deg: degrees or radians?
        :return: distance in metres
        """
        # todo implement !
        pass

    def set_reference_point(self, lon: float, lat: float, deg: bool = True):
        self._lon_ref_rad = (c.deg2rad*lon if deg else lon)
        self._lat_ref_rad = (c.deg2rad*lat if deg else lat)
        self.cos_lon, self.sin_lon = trig(self._lon_ref_rad)
        self.cos_lat, self.sin_lat = trig(self._lat_ref_rad)
        self.r_gaussian = wgs84_gaussian_radius(self._lat_ref_rad, deg=False)

    def _dot(self, lon, lat, deg):
        lon2_rad = (c.deg2rad*lon if deg else lon)
        cos_dlon = mt.cos(lon2_rad - self._lon_ref_rad)
        cos_lat2, sin_lat2 = trig(c.deg2rad*lat if deg else lat)
        return self.cos_lat*cos_lat2*cos_dlon + self.sin_lat*sin_lat2


def distance_wgs84(lon1: float, lat1: float, lon2: float, lat2: float,
                   deg: bool = True, method: str = 'great circle') -> float:
    """Vincenty's formula for the distance between two points on the reference
    ellipsoid.

    :param deg: degrees or radians?
    :param lon1: longitude of the first point :math:`\lambda_1`
    :param lat1: latitude of the first point :math:`\varphi_1`
    :param lon2: longitude of the second point :math:`\lambda_2`
    :param lat2: latitude of the second point :math:`\varphi_2`
    :param method: method to use: 'great circle', 'vincenty', ...
    :return: distance in metres
    """
    geo_dist = GeoDistance(lon1, lat1, deg=deg)
    if method == 'great circle':
        return geo_dist.great_circle(lon2, lat2, deg=deg)
    else:
        return geo_dist.vincenty(lon2, lat2, deg=deg)


# -----------------------------------------------------------------------------
# Topocentric coordinates

class Topocentric:
    """Useful to transform to topocentric coordinates from geographic or ECEF
    coordinates, given a fixed reference point.

    """
    def __init__(self, coords_ref: Sequence, c_type: str = 'geo'):
        """Initializes the instance

            :param coords_ref: coordinates of the reference point.
            :param c_type: type of coordinates for the reference point.
                Can be 'geo' (default) or 'xyz'.
            """
        if c_type == 'geo':
            self.geo_ref = coords_ref
            self.xyz_ref = geo2ecef(*coords_ref)
        else:  # ECEF reference
            self.geo_ref = ecef2geo(*coords_ref)
            self.xyz_ref = coords_ref

        self.kappa_e, self.kappa_n = diff_easting_northing(
            self.geo_ref[1], h_ref=self.geo_ref[2], deg=True)

        self._c_ln, self._s_ln = trig(self.geo_ref[0]*c.deg2rad)
        self._c_lt, self._s_lt = trig(self.geo_ref[1]*c.deg2rad)
        self._c_lt_c_ln = self._c_lt*self._c_ln
        self._c_lt_s_ln = self._c_lt*self._s_ln
        self._s_lt_c_ln = self._s_lt*self._c_ln
        self._s_lt_s_ln = self._s_lt*self._s_ln
        self.rot_ecef2topo = np.array(
            ((self._s_ln, self._c_ln, 0),
             (-self._s_lt_c_ln, -self._s_lt_s_ln, self._c_lt),
             (self._c_lt_c_ln, self._c_lt_s_ln, self._s_lt)))

    def __call__(self, coords: Sequence, c_type: str) -> Tuple:
        """Transforms to topocentric coordinates

        :param coords: coordinates. Three numbers or three 1d-arrays.
        :param c_type: type of coordinates to be used as input.
            'ecef': cartesian (ECEF),
            'geo': geographic -longitude, latitude and ellipsoidal height-.
            'delta_ecef': variation of ECEF coordinates,
            'delta_geo': variation of geographic coordinates.
        :return: east, north, up
        """
        if c_type == 'delta_geo':
            return self.delta_geo_to_topo(*coords)
        elif c_type == 'delta_ecef':
            return self.delta_ecef_to_topo(*coords)
        elif c_type == 'geo':
            return self.geo_to_topo(*coords)
        elif c_type == 'ecef':
            return self.ecef_to_topo(*coords)
        else:
            print(
                '   geo.Topocentric:  \"c_type\" must be \'geo\', \'ecef\', '
                '\'delta_geo\', \'delta_ecef\'.')
            return mt.nan, mt.nan, mt.nan

    def geo_to_topo(self, lon: floatOrArray, lat: floatOrArray,
                    h: floatOrArray) -> Tuple[floatOrArray]:
        return self.delta_geo_to_topo(
            lon - self.geo_ref[0], lat - self.geo_ref[1], h - self.geo_ref[2])

    def ecef_to_topo(self, x: floatOrArray, y: floatOrArray,
                     z: floatOrArray) -> Tuple[floatOrArray]:
        return self.delta_ecef_to_topo(
            x - self.xyz_ref[0], y - self.xyz_ref[1], z - self.xyz_ref[2])

    def delta_geo_to_topo(self, d_lon: floatOrArray, d_lat: floatOrArray,
                          d_h: floatOrArray) -> Tuple:
        return self.kappa_e*d_lon, self.kappa_n*d_lat, d_h

    def delta_ecef_to_topo(self, d_x: floatOrArray, d_y: floatOrArray,
                           d_z: floatOrArray) -> Tuple:
        return (-d_x*self._s_ln + d_y*self._c_ln,
                -d_x*self._s_lt_c_ln - d_y*self._s_lt_s_ln + d_z*self._c_lt,
                d_x*self._c_lt_c_ln + d_y*self._c_lt_s_ln + d_z*self._s_lt)

    def topo_to_geo(self, e: floatOrArray, n: floatOrArray, u: floatOrArray):
        return [a + b for a, b in zip(self.geo_ref,
                                      self.topo_to_delta_geo(e, n, u))]

    def topo_to_ecef(self, e: floatOrArray, n: floatOrArray,
                     u: floatOrArray) -> Tuple:
        return tuple(a + b for a, b in zip(self.xyz_ref,
                                           self.topo_to_delta_ecef(e, n, u)))

    def topo_to_delta_geo(self, e: floatOrArray, n: floatOrArray,
                          u: floatOrArray) -> Tuple:
        return e/self.kappa_e, n/self.kappa_n, u

    def topo_to_delta_ecef(self, e: floatOrArray, n: floatOrArray,
                           u: floatOrArray) -> Tuple:
        return (-e*self._s_ln - n*self._s_lt_c_ln + u*self._c_lt_c_ln,
                e*self._c_ln - n*self._s_lt_s_ln + u*self._c_lt_s_ln,
                n*self._c_lt + u*self._s_lt)

    def rotate_matrix(self, matrix: np.ndarray,
                      ecef2topo: bool = True) -> np.ndarray:
        if ecef2topo:
            return self.rot_ecef2topo.T.dot(matrix).dot(self.rot_ecef2topo)
        return self.rot_ecef2topo.dot(matrix).dot(self.rot_ecef2topo.T)
# -----------------------------------------------------------------------------


def lon_lat_displ_xyz(lon: float, lat: float, h: float,
                      disp: np.ndarray, get_pos: bool,
                      deg: bool = True) -> Tuple:
    r"""
    Latitude and longitude of a displaced point due to a spatial displacement.

    .. math::
        \Delta \vec{x} = m_\varphi \Delta \varphi \hat{e} +
            m_\lambda \Delta \lambda \hat{n} + \Delta h \hat{z}

        \text{where} \quad m_\varphi \equiv
        |\frac{\partial \vec{r}_{ecef}(\varphi, \lambda, h)}
        {\partial \varphi}|,\, m_\lambda \equiv
        |\frac{\partial \vec{r}_{ecef}(\varphi, \lambda, h)}{\partial \lambda}|
        \\
        \text{and } \{\hat{e}, \hat{n}, \hat{z}\}
        \text{ are the east, north and up unit vectors.}

    :param lon: longitude :math:`\varphi` of the point.
    :param lat: latitude :math:`\lambda` of the point.
    :param h: height :math:`h` in metres.
    :param disp: displacement vector :math:`\Delta \vec{x}` in metres and
        in ECEF (Geocentric Cartesian Coordinates).
    :param get_pos: controls if position vector is returned.
    :param deg: degrees or radians?
    :return: :math:`\Delta \varphi, \Delta \lambda, \Delta h,
        (\vec{r}(\varphi, \lambda, h))`.
    """
    cos_lon, sin_lon = trig(c.deg2rad*lon if deg else lon)
    cos_lat, sin_lat = trig(c.deg2rad*lat if deg else lat)
    # displacement in local coordinates
    d_rho_ = ecef2enu_matrix_trig(cos_lon, sin_lon, cos_lat, sin_lat).dot(disp)
    r_, m_phi, m_psi = ecef_aux2(cos_lon, sin_lon, cos_lat, sin_lat, h)
    lon_new = lon + d_rho_[0]/m_phi
    lat_new = lat + d_rho_[1]/m_psi
    if get_pos:
        return lon_new, lat_new, d_rho_[2], r_
    else:
        return lon_new, lat_new, d_rho_[1]/m_psi, d_rho_[2]


def reduced_latitude(lat: float, deg: bool = True) -> float:
    r"""Reduced latitude at a point (latitude on the auxiliary sphere)

    .. math::
        U = \text{arctan}((1 - f)\tan\lambda)

    :param lat: latitude :math:`\lambda` of the point.
    :param deg: degrees or radians?
    :return: reduced latitude :math:`U`
    """
    return _reduced_latitude_aux(*trig(c.deg2rad*lat if deg else lat))


def _reduced_latitude_aux(cos_lat, sin_lat):
    return mt.atan((1 - c.f)*sin_lat/cos_lat)
