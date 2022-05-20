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

"""Reproduces the API of networktools.geo, fixing some bugs in functions
that transform coordinates and removing unnecessary calls to sinus and cosinus
functions."""
from typing import Tuple
import math
import numpy as np
import geoproj.const as c
from geoproj.coords import (ecef2geo, geo2ecef,
    ecef2enu_matrix, ecef2enu_matrix_trig, east_north_up, _enu_trig,
    _east_trig, _north_trig, _up_trig)
from geoproj.aux import floatOrArray, floatOrInt


ecef2llh = ecef2geo
ecef2enu_rot = ecef2enu_matrix
ecef2enu_rot_tr = ecef2enu_matrix_trig
enu_tr = _enu_trig
east_tr = _east_trig
north_tr = _north_trig
up_tr = _up_trig


def excentricity(a0: floatOrInt, b0: floatOrInt):
    assert a0 >= b0, "a must be a superior value over b"
    aux = b0/a0
    return math.sqrt(1. - aux*aux)


def rad2deg(rlat: floatOrArray, rlon: floatOrArray) -> Tuple:
    return c.rad2deg*rlat, c.rad2deg*rlon


def deg2rad(dlat: floatOrArray, dlon: floatOrArray) -> Tuple:
    return c.deg2rad*dlat, c.deg2rad*dlon


def radius(latitude: float) -> Tuple:  # latitude in radians
    sin_lat = math.sin(latitude)
    aux = 1 - c.e2*sin_lat*sin_lat
    n = c.a/math.sqrt(aux)
    m = c.a*(1-c.e2)/math.pow(aux, 1.5)
    r = math.sqrt(m*n)
    return r, n, m


def llh2ecef(lat, lon, elev):
    return geo2ecef(lon, lat, h=elev, deg=False)


def ecef2neu(P0, dx, dy, dz):
    enu_vectors = east_north_up(P0[1], P0[0], deg=False, as_arrays=False)
    d_pos = (dx, dy, dz)
    enu_coords = tuple(sum(v*d for v, d in zip(vec, d_pos))
                       for vec in enu_vectors)
    return dict(N=enu_coords[1], E=enu_coords[0], U=enu_coords[2])


def get_from_ecef(ECEF_DATA):
    x = None
    y = None
    z = None
    try:
        #rprint("PRE ECEF: %s" %ECEF_DATA, flush =True)
        ecef = ECEF_DATA['ECEF']
        #gprint("ECEF: %s" %ecef)
        if 'X' in ecef:
            x = ecef['X']
            y = ecef['Y']
            z = ecef['Z']
        elif 'X_POS' in ecef:
            x = ecef['X_POS']
            y = ecef['Y_POS']
            z = ecef['Z_POS']
    except Exception as ex:
        print("get from ecef execption %s" % ex)
        raise ex
    return [x, y, z]


def get_vcv_matrix(VCV):
    return matrix_VCV(VCV['VCV_XX'], VCV['VCV_YY'], VCV['VCV_ZZ'],
                      VCV['VCV_XY'], VCV['VCV_XZ'], VCV['VCV_YZ'])


def rotate_vcv(R, VCV):
    return R.dot(VCV.dot(R.T))


def vcv2dict(C):
    return {'EE': C[0][0],
            'EN': C[0][1],
            'EU': C[0][2],
            'NN': C[1][1],
            'NU': C[1][2],
            'UU': C[2][2], }


def all_in_one_vcv(R, POSITION_VCV):
    vcv = get_vcv_matrix(POSITION_VCV)
    C = rotate_vcv(R, vcv)
    vcv_dict = vcv2dict(C)
    return vcv_dict


def matrix_VCV(Vxx, Vyy, Vzz, Vxy, Vxz, Vyz):
    return np.array([[Vxx, Vxy, Vxz],
                     [Vxy, Vyy, Vyz],
                     [Vxz, Vyz, Vzz]])
