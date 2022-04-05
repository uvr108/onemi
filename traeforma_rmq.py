from obspy.clients.seedlink import EasySeedLinkClient
# from rethinkdb import RethinkDB
# import requests
from mylibrary import request
from tasktools.taskloop import TaskLoop
# import asyncio
from random import seed
from random import random
import asyncio
from multiprocessing import Pool
# from data_rdb import Rethink_DBS
from rich import print
import pika
import json
import time
import calendar

# r = RethinkDB()
# r.connect( "10.54.217.85", 28015).repl()


class MyClient(EasySeedLinkClient):
    
    def on_data(self, trace):
        # print("Received", trace)
        network = trace.stats.network
        station = trace.stats.station
        channel = trace.stats.channel
        delta = trace.stats.delta
        npts = trace.stats.npts
        starttime = trace.stats.starttime.timestamp
        endtime = trace.stats.endtime.timestamp

        # print(trace.stats.endtime.year)
        # print(trace.stats.endtime.month)
        period = "%d%02d%02d" % (trace.stats.endtime.year, trace.stats.endtime.month, trace.stats.endtime.day)
        latency = calendar.timegm(time.gmtime()) - endtime 
        new_data = {
            # 'period':period, 
            'network':network,
            # 'npts': npts,  
            'station': station,
            # 'delta': delta, 
            # 'channel': channel, 
            # 'starttime': starttime,  
            # 'endtime': endtime
            'latency': latency
        }
        print("New data received", new_data)
        # r.db('csn').table('indicadores').insert(new_data).run()

        connection = pika.BlockingConnection( pika.ConnectionParameters(host='localhost') )
        
        channel = connection.channel()

        channel.basic_publish(
            exchange='logs', 
            routing_key='',
            # properties=pika.BasicProperties( content_type='text/json' ),
            body=json.dumps(new_data))

        connection.close()

async def extract_data(
    create, 
    select, 
    run, 
    stations, 
    component, 
    client, *args, **kwargs):
    if create:
        address = kwargs.get("address")
        print(f"Connecting to {address}")
        client = MyClient(address)
        create = False
        select = True
        print("Cliente creado y conectado", client)
    elif select:
        for network, station in stations:
            # print("Seleccionando estación", station)
            client.select_stream(network ,station, component)
            select = False
            run = True
    elif run:
        # print("Corriendo extractor de datos para", stations)
        try:
            client.run()
        except Exception as e:
            print(f"Falla en cliente {station}, {client} --> {e}")
            create = True
            run = False
            await asyncio.sleep(5)
    else:
        await asyncio.sleep(30)
    return [create, select, run, stations, component, client], kwargs

def main():

    # time.sleep(15)   
 
    stations=[]
    i=0
    for st in  request():
        if st[0] != 'GPS':
            stations.append([st[0],st[1]])
            # print(st[0],st[1].rstrip()) 

    tasks = []
    component = "HHZ"
    #args = [create, select, run, network, component, client]
    args = [True, False, False, stations, component, None]
    kwargs = {
        "address":'10.54.217.12:18000'
    }
    task  = TaskLoop(
        extract_data, 
        coro_args=args, 
        coro_kwargs=kwargs)
    tasks.append(task)
    task.create()
    loop = asyncio.get_event_loop()
    if not loop.is_running():
        loop.run_forever()
    
if __name__ == "__main__":

    main()


