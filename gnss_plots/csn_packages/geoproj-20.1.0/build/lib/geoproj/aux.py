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
from typing import Union, Tuple
import numpy as np
import math as mt
deg2rad = mt.pi/180.
floatOrArray = Union[float, np.ndarray]
floatOrInt = Union[float, int]


def trig(alpha: floatOrArray) -> Tuple:
    r"""Cosine and sine of an angle.

    :param alpha: angle :math:`\alpha` in radians.
    :return: :math:`\cos\alpha, \, \sin\alpha`
    """
    if isinstance(alpha, float):
        return mt.cos(alpha), mt.sin(alpha)
    else:
        return np.cos(alpha), np.sin(alpha)


def dms2deg(deg: int, mnt: int, sec: Union[float, int]) -> float:
    """(degrees, minutes, seconds) to floating point degrees.

    :param deg: degrees.
    :param mnt: minutes.
    :param sec: seconds.
    :return: degrees.
    """
    s = 1 if deg > 0 else -1
    return deg + s*(mnt + sec/60.)/60.


def deg2dms(d: Union[float, int]) -> Tuple:
    """Floating point degrees to (degrees, minutes, seconds).

    :param d: degrees.
    :return: degrees, minutes, seconds.
    """
    s = 1 if d > 0 else -1
    deg = int(d*s)
    aux = (d*s - deg)*60.
    mnt = int(aux)
    return s*deg, mnt, (aux - mnt)*60.


def dms2deg_arr(deg: np.ndarray, mnt: np.ndarray,
                sec: np.ndarray) -> np.ndarray:
    """Arrays: (degrees, minutes, seconds) to floating point degrees.

    :param deg: degrees.
    :param mnt: minutes.
    :param sec: seconds.
    :return: degrees.
    """
    s = 2*np.array(deg > 0, dtype=int) - 1
    return deg + s*(mnt + sec/60.)/60.


def deg2dms_arr(d: np.ndarray) -> Tuple:
    """Arrays: floating point degrees to (degrees, minutes, seconds).

    :param d: degrees.
    :return: degrees, minutes, seconds.
    """
    s = 2*np.array(d > 0, dtype=int) - 1
    deg = np.array(d*s, dtype=int)
    aux = (d*s - deg)*60.
    mnt = np.array(aux, dtype=int)
    return s*deg, mnt, (aux - mnt)*60.


def are_close(x, y, epsilon=1.e-8):
    return abs(x-y) < epsilon


def are_close_vec(x, y, epsilon=1.e-8):
    aux = x - y
    return aux.dot(aux) < epsilon*epsilon


def are_all_close(x, y, epsilon=1.e-8):
    return np.abs(x - y).max() < epsilon


def are_close_lonlat(lonlat1, lonlat2, epsilon=1.e-8):
    cos_lat = mt.cos(deg2rad*lonlat1[1])
    return (are_close(lonlat1[0]*cos_lat, lonlat2[0]*cos_lat,
                      epsilon=epsilon) and
            are_close(lonlat1[1], lonlat2[1], epsilon=epsilon))


def are_all_close_lonlat(lonlat1, lonlat2, epsilon=1.e-8):
    cos_lat = np.cos(deg2rad*lonlat1[1])
    return (are_all_close(lonlat1[0]*cos_lat, lonlat2[0]*cos_lat,
                          epsilon=epsilon) and
            are_all_close(lonlat1[1], lonlat2[1], epsilon=epsilon))
