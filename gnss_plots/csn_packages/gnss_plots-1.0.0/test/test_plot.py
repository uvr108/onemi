from gnss_redb.parser import GsofEcefParser
from gnss_plots import GnssPlot

host = '10.54.217.15'
port = 28015
db_name = "collector"
station = "BING"
# ---------------------------------------------------------------------------
g_plot = GnssPlot(host, port, db_name=db_name, parsers=[GsofEcefParser])
g_plot.create_single_plot(station, output_file=f"./fig_{station}.png", dpi=80)
# ---------------------------------------------------------------------------
# NOTA: se puede entregar los parametros "start" y "end" a "create_single_plot"
#       para definir una ventana de tiempo arbitraria

# Consultar estaciones y tablas
print(g_plot.get_stations_set())
print(g_plot.get_tables_set())
