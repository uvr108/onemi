# -*- coding: utf-8 -*-
"""
*gnss_redb* is a package that allows to get data from a RethinkDB databases
with GNSS data. Designed for internal use in the CSN.

Copyright (c) 2020 Centro Sismológico Nacional, Universidad de Chile

This file is part of "gnss_redb".

"gnss_redb" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"gnss_redb" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "gnss_redb". If not, see <https://www.gnu.org/licenses/>.
"""
from math import sqrt
import numpy as np
from geoproj.coords import ecef2geo, geo2ecef, east_north_up
from gnss_redb.path import ref_coords_file
__docformat__ = 'reStructuredText en'


def read_ref_coords():
    dict_coords = dict()
    with open(ref_coords_file, 'r') as reader:
        line = reader.readline().rstrip()
        while line != '':
            aux = line.split()
            h = 0.0 if (aux[3] == 'NaN') else float(aux[3])
            dict_coords[aux[0]] = (float(aux[1]), float(aux[2]), h)
            line = reader.readline()
    return dict_coords


class StationsCoordsConverter:

    def __init__(self):
        self._coords_conv = dict()

    def stations(self):
        return tuple(self._coords_conv.keys())

    def code_coords_generator(self):
        return ((code, conv.ref_llh)
                for code, conv in self._coords_conv.items())

    def add_station(self, code, ref_coords, is_llh=True, compute_coeff=False):
        coords_conv = self._coords_conv.get(code)
        station_added = (coords_conv is not None)
        if station_added:
            if coords_conv.ref_coords_not_set:
                coords_conv.set_ref_coords(ref_coords, is_llh=is_llh,
                                           compute_coeff=compute_coeff)
            elif compute_coeff and not coords_conv.coefficients_are_set:
                coords_conv.set_coefficients()
        else:
            self._coords_conv[code] = CoordsConverter(
                ref_coords, is_llh, compute_coeff)

    def ecef_2_enu(self, station, xyz, cov_xyz):
        try:
            cc = self._coords_conv[station]
            return cc.ecef2enu(*xyz), cc.cov_ecef_to_std_enu(cov_xyz)
        except KeyError:
            return None

    def get_ref_llh(self, station):
        return self._coords_conv[station].ref_llh

    def get_ref_ecef(self, station):
        return self._coords_conv[station].ref_ecef

    def ref_coords_not_set(self, station):
        """Returns station_added, ref_coords_set, coeffs_set"""
        cc = self._coords_conv.get(station)
        if cc is None:
            return True
        return cc.ref_coords_not_set


class CoordsConverter:

    def __init__(self, ref_coords, is_llh=True, compute_coeff=False):
        self._not_set = True
        self._coeffs_are_set = False
        self._enu_tuples = None
        self._var_e_coeff = None
        self._var_n_coeff = None
        self._var_u_coeff = None
        self._ref_llh = self._ref_xyz = None
        self.set_ref_coords(ref_coords, is_llh=is_llh,
                            compute_coeff=compute_coeff)

    def set_ref_coords(self, ref_coords, is_llh=True, compute_coeff=False):
        if self._not_set:
            self._not_set = (ref_coords is None)
            if self._not_set:
                self._ref_llh = self._ref_xyz = np.full(3, np.nan)
            else:
                self._not_set = False
                if is_llh:
                    self._ref_llh = ref_coords
                    self._ref_xyz = geo2ecef(*ref_coords)
                else:
                    self._ref_xyz = ref_coords
                    self._ref_llh = ecef2geo(*ref_coords)
                if compute_coeff:
                    self.set_coefficients()

    def set_coefficients(self):
            self._coeffs_are_set = True
            self._enu_tuples = east_north_up(self._ref_llh[0],
                                             self._ref_llh[1], as_arrays=False)
            self._var_e_coeff = self.std_coeff(self._enu_tuples[0])
            self._var_n_coeff = self.std_coeff(self._enu_tuples[1])
            self._var_u_coeff = self.std_coeff(self._enu_tuples[2])

    def ecef2enu(self, *xyz):
        if self._not_set:
            return self._ref_xyz
        d_xyz = tuple(xyz[k] - self._ref_xyz[k] for k in range(3))
        return tuple(sum(a*b for a, b in zip(self._enu_tuples[k], d_xyz))
                     for k in range(3))

    def cov_ecef_to_std_enu(self, cov_xyz):
        if self._not_set:
            return self._ref_xyz
        return (sqrt(self.dot(self._var_e_coeff, cov_xyz)),
                sqrt(self.dot(self._var_n_coeff, cov_xyz)),
                sqrt(self.dot(self._var_u_coeff, cov_xyz)))

    @property
    def ref_coords_not_set(self):
        return self._not_set

    @property
    def coefficients_are_set(self):
        return self._coeffs_are_set

    @property
    def ref_llh(self):
        return self._ref_llh

    @property
    def ref_ecef(self):
        return self._ref_xyz

    @staticmethod
    def std_coeff(vec):
        return (vec[0]*vec[0], 2*vec[0]*vec[1], 2*vec[0]*vec[2],
                vec[1]*vec[1], 2*vec[1]*vec[2], vec[2]*vec[2])

    @staticmethod
    def dot(vec_x, vec_y):
        return sum(a*b for a, b in zip(vec_x, vec_y))


def id_func(*x):
    return x

