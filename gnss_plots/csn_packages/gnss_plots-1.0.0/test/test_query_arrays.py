from gnss_plots import GnssPlot

from datetime import datetime, timezone
from gnss_redb.db import ServerDB
from gnss_redb.parser import GsofEcefParser, GsofEnuParser, PppEnuParser


addresses = {
    'Atlas': ['10.54.217.15', 28015, [('collector', [GsofEcefParser])]],
    # 'Bellaco': ['10.54.218.66', 28045],
    'Test': ['10.54.218.39', 28015, [('enu_data', [GsofEnuParser]),
                                     ('amqp_collector', [PppEnuParser])]]
}
host, port, dbs_info = addresses['Atlas']
db_name = dbs_info[0][0]  # "collector"
parsers = dbs_info[0][1]

g_plot = GnssPlot(host, port, db_name=db_name, parsers=parsers)
g_plot.set_current_db(db_name=db_name)
tables = list(g_plot.get_tables_set())
table = tables[0]
t_array, enu, std_enu = g_plot.time_interval_arrays(
    table, start="2022-01-30T02:30:21+00:00", end="2022-01-30T03:30:21+00:00")

print(t_array)
print(enu[0])
