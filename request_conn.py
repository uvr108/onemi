import requests
from random import seed
from random import random
import configparser

def request():
    try:

        home_project = '/home/uvergara/Projectos/Python/onemi'

        config = configparser.ConfigParser()
        config.read(f"{home_project}/config.ini")

        serverinfo = config["SERVERCONFIG"]

        restfull = format(serverinfo["restfull"])


        response = requests.get(f"http://{restfull}:3000/api/estacion/-20/-60/0/true/false/true")
        response.raise_for_status()

        resp = response.json()

        st = []

        for r in resp:
            net = r['Network']['nombre']
            if len(net) == 1:
                net = net + ' ';
            st.append([net.rstrip(), r['codigo']])

        return st

    except requests.exceptions.HTTPError as errh:
        print(errh)
    except requests.exceptions.ConnectionError as errc:
        print(errc)
    except requests.exceptions.Timeout as errt:
        print(errt)
    except requests.exceptions.RequestException as err:
        print(err)

def main():


    stations=[]
    i=0
    for st in  request():
        print(st[0],st[1])

if __name__ == "__main__":

    main()
