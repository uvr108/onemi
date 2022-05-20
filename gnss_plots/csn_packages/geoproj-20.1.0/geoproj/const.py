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
"""Geographical constants and auxiliary quantities."""
# angle transformation: degrees to radians and radians to degrees
deg2rad = 0.017453292519943295  # Pi/180
rad2deg = 57.29577951308232  # 180/Pi
#  WGS84 reference ellipsoid constants
a = 6378137  # [m]
e2 = 6.69437999014e-3
f = 0.003352810664746797  # 1.0-sqrt(1.0-e2)
b = a*(1-f)
# Transverse Mercator projection constants
e4 = e2*e2
e6 = e4*e2
A0 = 1.0 - e2/4 - 3*e4/64 - 5*e6/256
A2 = (3/8.0)*(e2 + e4/4 + 15*e6/128)
A4 = (15/256.)*(e4 + 3*e6/4)
A6 = 35*e6/3072
ko = 1.0
no = (a-b)/(a+b)
no2 = no*no
g_aux = a*(1-no)*(1-no2)*(1+no2*(9/4.0 + (225/64.)*no2))
# auxiliary constants
l_km = 1000.
delta = 1.e-8
delta_half = 0.5*delta
