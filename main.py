from typing import Optional
from fastapi import FastAPI
from request_conn import request
import requests
import calendar
import time
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import pika
import threading
import json

from datetime import datetime

from multiprocessing import Process

from random import uniform
# from graficos import generagraph
# from gnss_plot import generaplot

from rethinkdb import RethinkDB

sismologicas = {}

def getstat():
   
    connection = pika.BlockingConnection(
    pika.ConnectionParameters(host='localhost'))
    channel = connection.channel()

    channel.exchange_declare(exchange='logs', exchange_type='fanout')

    result = channel.queue_declare(queue='', exclusive=True)
    queue_name = result.method.queue

    channel.queue_bind(exchange='logs', queue=queue_name)

    print(' [*] Waiting for logs. To exit press CTRL+C')

    def callback(ch, method, properties, body):
        stats = json.loads(body)
        
        try:
            sismologicas[stats['station']] =  [stats['latency'], stats['endtime']]
            # print(sismologicas)
        except Exception as e:
            print("Oops!", e.__class__, "occurred.")
        finally:
            pass  

    channel.basic_consume(
        queue=queue_name, on_message_callback=callback, auto_ack=True)

    channel.start_consuming()
    

def stations():

    return sismologicas 

def graph(net, stat):
    return generagraph(net,stat)

def gps(stat):
    return generaplot(stat)

def recibe():

    r = RethinkDB()
    rethink_conn = r.connect(host='localhost', port=28015)
   
    table_changes = r.db('csn').table('ULI_SOFT').changes()

    cont = 0
    for change in table_changes.run(rethink_conn):
        new = change['new_val']
        try:
            stat = new['station'] 
            retraso = new['info']['delta_time']
            sismologicas[stat] = ['',  retraso ] 
            # sismologicas[stat] = new['info'] 
        except Exception as e:
            print("Oops!", e.__class__, "occurred.")
        finally:
            pass
 

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"Hello": "World"}

@app.get("/getstat")
def mostra_stat():
    thread = threading.Thread(target=getstat(), args=()) 
    thread.start()
    return {"Hello": "World"}

@app.get("/recibe")
def mostra_recibe():
    thread = threading.Thread(target=recibe(), args=()) 
    thread.start()
    return {"Hello": "World"}

@app.get("/stations")
def mostra_stations():
    return stations() 


@app.get("/sismo/{net}/{stat}")
def mostra_graph(net: str, stat: str ):
    return graph(net,stat) 

@app.get("/gps/{stat}")
def gps(stat: str):
    return gps(stat)

