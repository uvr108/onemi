import requests
from random import seed
from random import random

def request():
    try:

        response = requests.get("http://10.54.217.85:3000/api/estacion")
        response.raise_for_status()

        resp = response.json()

        st = []

        for r in resp:
            net = r['Network']['nombre']
            if len(net) == 1:
                net = net + ' ';
            st.append([net.rstrip(), r['code']])

        return st

    except requests.exceptions.HTTPError as errh:
        print(errh)
    except requests.exceptions.ConnectionError as errc:
        print(errc)
    except requests.exceptions.Timeout as errt:
        print(errt)
    except requests.exceptions.RequestException as err:
        print(err)

