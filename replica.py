import click
import os
from enum import Enum
# from graficar import grafica_data
# from conexion_rabbit import AsyncRabbit

"""
Obtener la lista de estaciones
"""
import asyncio 
import csv
from rich import print
from pathlib import Path 
from rethinkdb import RethinkDB
rdb = RethinkDB()
from tasktools.taskloop import TaskLoop
# from slugify import slugify

from data_rdb import Rethink_DBS
from networktools.time import get_datetime_di
from networktools.geo import (radius, deg2rad, ecef2llh, llh2ecef)
from networktools.library import geojson2json
from data_geo import GeoJSONData

import configparser

class DBSend(Enum):
    CREATE = 1
    CONNECT = 2
    SEND = 3  


class DBStep(Enum):
    CREATE = 1
    CONNECT = 2
    COLLECT = 3  

class Mode(Enum):
    GSOF = 1
    SEED = 2


import numpy as np

def get_llh(data):
    data['ECEF_X'] = float(data['ECEF_X'])
    data['ECEF_Y'] = float(data['ECEF_Y'])
    data['ECEF_Z'] = float(data['ECEF_Z'])
    x = data['ECEF_X']
    y = data['ECEF_Y']
    z = data['ECEF_Z']
    (lat, lon, h) = ecef2llh(x, y, z)
    data["POSITION"] = {
        "ECEF"  : {
            "X": data['ECEF_X'],
            "Y": data['ECEF_Y'],
            "Z": data['ECEF_Z']
        },

        'llh': {'lat': lat, 'lon': lon, 'z': h}
    }
    # data["ECEF"] =  {
    #         "X": data['ECEF_X'],
    #         "Y": data['ECEF_Y'],
    #         "Z": data['ECEF_Z']
    #     },



def add_process_instance(station):
    CODE = station['protocol'].upper()
    station_code = station['code']
    kwargs = dict()
    kwargs['code'] = CODE
    kwargs['station'] = station
    kwargs['position'] = station["POSITION"]
    kwargs['log_path'] = str(Path("./tmp/log/db") / 'geo_json_data')
    process_instance = GeoJSONData(**kwargs)
    return process_instance


def read_csv(filename:Path, mode:Mode):
    dataset = {}
    if filename.exists():
        with open(filename,'r') as f:
            reader = csv.DictReader(f, delimiter=';')
            for row in reader:
                code = row.get("code")
                if mode == Mode.GSOF:
                    row["table_name"] = f"{code}_{mode.name}"
                    get_llh(row)
                    process = add_process_instance(row)
                    row["process_instance"] = process
                dataset[code] = row
    else:
        print("No existe", filename)
    return dataset

import aiofiles as aiof
from datetime import datetime
import ujson as json
import aiofiles

async def process_data(
        queue:asyncio.Queue, 
        dataset, control, db_insta, **kwargs):
    now = datetime.utcnow()
    base_path = kwargs.get("base_path")
    print("Getting from queue and send")
    if control == DBSend.CREATE:
        opts = kwargs.get("opts")
        db_insta = Rethink_DBS(**opts)
        control = DBSend.CONNECT
    if control == DBSend.CONNECT:
        fail = False
        opts = kwargs.get("opts")
        dbname = opts.get("dbname")
        try:
            await asyncio.wait_for(
                db_insta.async_connect(),
                timeout=10)
            await db_insta.list_dbs()
            await db_insta.list_tables()
            control = DBSend.SEND
        except asyncio.TimeoutError as e:
            print(f"Falla al intentar conectar, {e}")
            control = DBSend.CONNECT
            fail = True
            await asyncio.sleep(30)
        if not fail:
            control = DBSend.SEND
    if not base_path.exists():
        base_path.mkdir(parents=True, exist_ok=True)
    print("Queue len", queue.qsize())

    if not queue.empty():
        print("Control", control)

        for i in range(queue.qsize()):
            station, info = await queue.get()
            # path_station = base_path / f"{station}_{now.isoformat()}.json"
            # print(now, "Saving data---", path_station, len(dataset))
            # async with aiofiles.open(path_station,"w") as f:
            #     dumps = json.dumps(geodataset, f, default=lambda e:
            #               e.isoformat())
            #     await f.write(dumps)
            try:
                #pngpath, buf = grafica_data(geodataset, station)
                #print(pngpath, buf)
                mensaje = {
                    "station": station,
                    "network": "GPS",
                    "info": info,
                    # "latency": info['delta_time'],
                    #"dataset": geodataset,
                }
                """
                mensaje -> grafico (bytes que son la imagen)
                Cargar blob a imagen en frontend (con javascript")
                var blob = new Blob([bytes], {type: 'image/bmp'});
                // Use createObjectURL to make a URL for the blob
                var image = new Image();
                image.src = URL.createObjectURL(blob);
                """
                # ENVIAR MENSAJE CON GRAFICAS
                if control == DBSend.SEND:
                    #await rabbit.send(mensaje)
                    destiny = kwargs.get("destiny")
                    # print(destiny, "=> ",mensaje)
                    result = await db_insta.save_data(destiny,
                                                      mensaje)
                    #print("Resultado de guardar dato", result)
            except Exception as e:
                print("Error",e)
                await db_insta.close()

    await asyncio.sleep(10)
    return (queue, dataset, control, db_insta), kwargs

"""
Crear, conectar y obtener data entre tramos.
"""

async def db_manage(control, db_insta,  dataset, di, queue, **kwargs):
    delta_time = int(os.environ.get("DELTA_TIME", 10))
    if control == DBStep.CREATE:
        opts = kwargs.get("opts")
        db_insta = Rethink_DBS(**opts)
        control = DBStep.CONNECT
    elif control == DBStep.CONNECT:
        fail = False
        opts = kwargs.get("opts")
        dbname = opts.get("dbname")
        try:
            await asyncio.wait_for(db_insta.async_connect(),timeout=10)
            await db_insta.list_dbs()
            await db_insta.list_tables()
        except asyncio.TimeoutError as e:
            print(f"Falla al intentar conectar, {e}")
            control = DBStep.CONNECT
            fail = True
            await asyncio.sleep(30)
        if not fail:
            control = DBStep.COLLECT
    if control == DBStep.COLLECT:
        key = "DT_GEN"
        filter_opt = {'left_bound': 'open', 'index': key}
        df = rdb.iso8601(get_datetime_di(delta=0))
        try:
            for code, station in dataset.items():
                # try:
                #     print([di,df])
                #     print("Obteniendo...",[di,df],
                #           (df-di))
                #     diff =  await (df-di).run(db_insta.session)
                #     print(diff)
                # except Exception as e:
                #     print("eeor", e)
                #     pass
                table_name = station.get("table_name")
                cursor = await db_insta.get_data_filter(
                        table_name,
                        [di, df],
                        filter_opt,
                        key)
                geodataset = []
                mean = []
                recv = []
                for data in cursor:
                    process_instance = station.get("process_instance")
                    
                    try:
                        if process_instance:
                            if "POSITION_VCV" in data:
                                #geodata = process_instance.manage_data(data)
                                #geojson = geojson2json(geodata,
                                #                       destiny='db')
                                #geojson["time"] = {
                                #    "recv": data.get("DT_RECV"),
                                #    "delta": data.get("DELTA_TIME")
                                #}
                                recv.append(data.get("DT_RECV"))
                                mean.append(data.get("DELTA_TIME"))
                                #geodataset.append(geojson)
                            else:
                                now = datetime.utcnow()
                                recv.append(data.get("DT_RECV"))
                                mean.append(data.get("DELTA_TIME"))
                                print(now, "No VCV",station)
                                print(data)
                        else:
                            print("No hay process instance para...",
                                  station, process_instance)
                            await asyncio.sleep(10)
                    except Exception as e:
                        print("Data error", data)
                        await asyncio.sleep(10)
                        raise e
                if mean and recv:
                    array = np.array(mean)
                    now = datetime.utcnow()
                    start = recv[0]
                    end = recv[-1]
                    delta = np.mean(array)
                    std = np.std(mean)
                    info = {
                        "start": start.isoformat(),
                        "end": end.isoformat(),
                        "delta_time":delta,
                        "std": std,"elements":len(mean)
                    }
                    await queue.put((code, info))
        except asyncio.TimeoutError as te:
            control = DBStep.CONNECT
        di = df
    if control == DBStep.COLLECT:
        await asyncio.sleep(delta_time)
    return (control, db_insta, dataset, di, queue), kwargs

if __name__ == "__main__":

    home_project = '/home/ulises/Documentos/Projectos/Python/onemi'

    config = configparser.ConfigParser()
    config.read(f"{home_project}/config.ini")
    serverinfo = config["SERVERCONFIG"]

    rethink = format(serverinfo["rethinkdb"])
    collector = format(serverinfo["collector"])
    rethinkport = format(serverinfo["rethinkport"])

    print(f'rethink -> {rethink} collector-> {collector}')

    delta_time = int(os.environ.get("DELTA_TIME", 300))

    filename = Path(__file__).parent / "fuentes/station.csv"
    mode = Mode.GSOF
    dataset = read_csv(filename, mode)
    opts = {
        "host": f'{rethink}',
        "port": f'{rethinkport}',
        "dbname": f'{collector}'
    }
    di = rdb.iso8601(get_datetime_di(delta=delta_time))
    queue = asyncio.Queue()

    args = [DBStep.CREATE, None, dataset, di, queue]
    kwargs = {"opts": opts}
    task = TaskLoop(db_manage, coro_args=args, coro_kwargs=kwargs)
    task.create()

    # mq_opts = {
    #     "username": "enugeojson_prod",
    #     "password": "geognss_prod", 
    #     "host": "10.54.217.95",
    #     "vhost": "gpsdata",
    #     "exchange": "graficas",
    #     "queue_name": "stations",
    #     "routing_key": "plots",
    # }
    db_insta = None#AsyncRabbit(**opts)
    control = DBSend.CREATE
    args = [queue, dataset, control, db_insta]
    base_path = Path(__file__).absolute().parent / f"{home_project}/data"
    opts_destiny = {
        "host": "127.0.0.1",
        "port": 28015,
        "dbname": "csn"
    }
    table_destiny = "ULI_SOFT"

    task = TaskLoop(process_data, coro_args=args,
                    coro_kwargs={"base_path":base_path, "opts":
                                 opts_destiny, 
                                 "destiny": table_destiny})
    task.create()

    loop = asyncio.get_event_loop()

    loop.run_forever()
