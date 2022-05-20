from pathlib import Path
from rethinkdb import RethinkDB, errors
import numpy as np
from gnss_redb.db import datetime_window
this_dir = Path(__file__).parent.absolute()

# -----------------
# server_name = "Test"
server_name = "Atlas"
# recipients = "test-geodesia"
end = "2021-12-29T10:00:00Z"
# end = "now-hour"
dt_start, dt_end = datetime_window(end, window_dict=dict(days=1))

server_dicts = {
    "Atlas": {
      "host": "10.54.217.15", "port": 28015,
      "dbs": {
        "collector": ["GSOF-ECEF"]}
    },

    "Test": {
      "host": "10.54.218.39", "port": 28015,
      "dbs": {
        "amqp_collector": ["PPP-ENU"]}
    },

    "Bellaco": {
      "host": "10.54.218.66", "port": 28045,
      "dbs": {
        "collector": ["GSOF-ECEF"]}
    }
}

srv_dict = server_dicts[server_name]


redb = RethinkDB()
connection = redb.connect(srv_dict["host"], srv_dict["port"])

db_names = redb.db_list().run(connection)
print(db_names)
db_name = "collector"
tables = redb.db(db_name).table_list().run(connection)

# *********************
scan_all = True
scan_one = True
n_day = 86400

no_loss = list()
no_data = list()
some_loss = dict()
loss_percentage = list()

if scan_all:
    print(tables)
    print("_______________________________________________________")
    print("                       PERDIDAS")
    #
    for tb_name in tables:
        tb = redb.db(db_name).table(tb_name)
        try:
            n_query = tb.between(dt_start, dt_end, index="DT_GEN"
                                 ).count().run(connection, time_format='raw')
            if n_query == 0:
                no_data.append(tb_name)
                loss = 100.
            elif n_query == n_day:
                no_loss.append(tb_name)
                loss = 0.
            else:
                loss = 100.*(1- n_query/n_day)
                some_loss[tb_name] = loss
                print("{:}  {:.2f} %".format(tb_name, loss))
        except errors.ReqlOpFailedError as err:
            print(err)


    print("_______________________________________________________")
    print("                       SIN DATOS")
    for table_name in no_data:
        print(table_name)

    print("_______________________________________________________")
    print("                       SIN PERDIDAS")
    for table_name in no_loss:
        print(table_name)
if scan_one:
    table_name = "TMCO_GSOF"
    # table_name = "PARC_GSOF"
    # table_name = "CRZL_GSOF"
    tb = redb.db(db_name).table(table_name)
    print(f" Checking {table_name} ...")
    print(f"  {table_name} in tables = {table_name in tables}")

    result = tb.between(dt_start, dt_end, index="DT_GEN").pluck("DT_GEN").run(connection, time_format='raw')
    # result = tb.limit(3).run(connection, time_format='raw')

    ts = np.array(list(doc["DT_GEN"]["epoch_time"] for doc in result))
    print(f"N data ({table_name}): {ts.size}")
    print(f"Loss   ({table_name}): {100*(1. - ts.size/n_day):.2f}")
    if ts.size > 0:
        ts.sort()
        ts -= ts[0]
        aux = np.abs(ts - np.round(ts))
        print(aux)
        print(aux.max())
