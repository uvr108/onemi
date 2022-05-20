import numpy as np
from geoproj.coords import geo2ecef, ecef2enu_matrix
from gnss_redb.coords import CoordsConverter

n = 100
lon_array = np.random.uniform(-180., 180., n)
lat_array = np.random.uniform(-90, 90, n)
h_array = np.random.uniform(-1e5, 1e5, n)

for k in range(n):
    lon = lon_array[k]
    lat = lat_array[k]
    h = h_array[k]
    x, y, z = geo2ecef(lon, lat, h)
    rot = ecef2enu_matrix(lon, lat)
    d_xyz = np.random.uniform(-1000, 1000, 3)
    enu_0 = rot.dot(d_xyz)
    sta = CoordsConverter((lon, lat, h), compute_coeff=True)
    enu = sta.ecef2enu(x + d_xyz[0], y + d_xyz[1], z + d_xyz[2])
    diff = max(a - b for a, b in zip(enu, enu_0))
    assert abs(diff) < 1.e-6
