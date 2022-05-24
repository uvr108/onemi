#!/usr/bin/env python
import glob
from obspy.core import read
from obspy.core import UTCDateTime
import configparser
from pathlib import Path
import os


def generagraph(tipo,net,stat):

    ahora = UTCDateTime()
    timestamp = UTCDateTime().timestamp

    config = configparser.ConfigParser()
    config.read("config.ini")

    # userinfo = config["USERINFO"]
    serverinfo = config["DIRECTORIES"]

    yr = ahora.year
    jl = "%03d" % ahora.julday

    image = f'{timestamp}.png'

    pat = f"*{yr}.{jl}"
    
    if tipo == 'CSN':

        st = read(f"/mnt/seiscomp/data/{yr}/{net}/{stat}/HHE.D/{net}*{pat}")
        st += read(f"/mnt/seiscomp/data/{yr}/{net}/{stat}/HHN.D/{net}*{pat}")
        st += read(f"/mnt/seiscomp/data/{yr}/{net}/{stat}/HHZ.D/{net}*{pat}")

    else:

        st = read(f"/mnt/seiscomp/data/{yr}/{net}/{stat}/HNE.D/{net}*{pat}")
        st += read(f"/mnt/seiscomp/data/{yr}/{net}/{stat}/HNN.D/{net}*{pat}")
        st += read(f"/mnt/seiscomp/data/{yr}/{net}/{stat}/HNZ.D/{net}*{pat}")

    tr = st[0]
    dt = st[0].stats.starttime

   
    ruta = f'{serverinfo["assets"]}/img/sismos/{net}/{stat}/{jl}'

    print(ruta)

    if Path(ruta).is_dir():
        pass 
    else:
        os.makedirs(ruta, exist_ok=True) 
    st.plot(outfile=f'{ruta}/{image}', starttime=ahora-310, endtime=ahora-10, size=(1100, 600)) 
    return {'image': image , 'jl': jl}

def main():

     print(generagraph('CSN','CX','PB01'))
    # print(generagraph('','C1','T14A'))

if __name__ == "__main__":

    main()


