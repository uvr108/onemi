from math import sqrt
import numpy as np
from scipy.stats import ortho_group
from geoproj.coords import geo2ecef, ecef2enu_matrix
from gnss_redb.coords import read_ref_coords, StationsCoordsConverter

dict_coords = read_ref_coords()
stations_cc = StationsCoordsConverter()


for station, ref_coords in dict_coords.items():
    stations_cc.add_station(station, ref_coords, compute_coeff=True)

for station, ref_coords in dict_coords.items():
    xyz = geo2ecef(*ref_coords)
    assert max(abs(a-b) for a, b in
               zip(xyz, stations_cc.get_ref_ecef(station))) < 0.0001
    assert max(abs(a-b) for a, b in
               zip(ref_coords, stations_cc.get_ref_llh(station))) < 1.e-9

    rot = ecef2enu_matrix(ref_coords[0], ref_coords[1])
    d_xyz = np.random.normal(0., 1000., (100, 3))
    xyz_vec = np.array(xyz)
    for k in range(d_xyz.shape[0]):
        rot_aux = ortho_group.rvs(3)
        cov_xyz = rot_aux.dot(np.diag(np.random.uniform(0., 0.1, 3))).dot(rot_aux.T)
        enu_aux = rot.dot(d_xyz[k])
        xyz_aux = xyz_vec + d_xyz[k]
        enu, std_enu = stations_cc.ecef_2_enu(
            station, xyz_aux, (cov_xyz[0, 0], cov_xyz[0, 1], cov_xyz[0, 2],
                               cov_xyz[1, 1], cov_xyz[1, 2], cov_xyz[2, 2]))
        assert max(abs(a-b) for a, b in zip(enu_aux, enu)) < 1.e-8
        cov_mtx_enu = rot.dot(cov_xyz).dot(rot.T)
        for i in range(3):
            assert abs(sqrt(cov_mtx_enu[i, i]) - std_enu[i]) < 1.e-8
