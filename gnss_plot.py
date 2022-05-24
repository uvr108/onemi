from gnss_redb.parser import GsofEcefParser
from gnss_plots import GnssPlot
import configparser
from pathlib import Path
from obspy.core import UTCDateTime
import os

host = '10.54.218.39'
port = 28015
db_name = "collector"

def generaplot(stat):

    ahora = UTCDateTime()
    timestamp = UTCDateTime().timestamp

    config = configparser.ConfigParser()
    config.read("config.ini")

    userinfo = config["USERINFO"]
    serverinfo = config["SERVERCONFIG"]
    directories = config["DIRECTORIES"]

    print("Password is {}".format(userinfo["password"]))
    print("Host is {}".format(serverinfo["host_gnss"]))

    yr = ahora.year
    jl = "%03d" % ahora.julday

    image = f'{timestamp}.png'

    ruta = f'{directories["assets"]}/img/gps/{stat}/{jl}'
    
    print(ruta)
 
    if Path(ruta).is_dir():
        pass
    else:
        os.makedirs(ruta, exist_ok=True)

    g_plot = GnssPlot(host, port, db_name=db_name, parsers=[GsofEcefParser])
    g_plot.create_single_plot(stat, output_file=f"{ruta}/{image}", dpi=80)
    
    return {'image': image , 'jl': jl}



def main():

    print(generaplot('BING'))
    

if __name__ == "__main__":

    main()

