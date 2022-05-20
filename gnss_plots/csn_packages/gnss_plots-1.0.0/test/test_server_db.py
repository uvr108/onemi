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
db_names = tuple(x[0] for x in dbs_info)
server_db = ServerDB(host, port)
server_db.do_connect()

for db_name, parsers in dbs_info:
    server_db.new_db_manager(db_name, parsers, update_table_set=False)

server_db.update_table_set()

print(' Active DBs:')
print(server_db.active_dbs)
for db_name in db_names:
    print('\n Active tables: ', db_name)
    print(server_db.active_tables_in_db(db_name))
    print('\n Table-Station pairs: ', db_name)
    for table, station in server_db.table_station_pairs(db_name):
        print('{}  -- {}'.format(table, station))

db_name = db_names[0]
active_tables_0 = tuple(server_db.active_tables_in_db(db_name))
print(f"N tables {len(active_tables_0)}")
table = active_tables_0[9]

gen = server_db.sample_values_generator(db_name, table, 10)
print('\n  DB: {}  --  table: {}'.format(db_name, table))
for x in gen:
    print(x)

print('\n\n --------------------------------------')

end = datetime.now(tz=timezone.utc).timestamp() - 360000
# end = 1602642310
start = end - 10
start = datetime.fromtimestamp(start, tz=timezone.utc)
end = datetime.fromtimestamp(end, tz=timezone.utc)

gen = server_db.time_interval_generator(db_name, table, start, end)
print('\n  DB: {}  --  table: {}'.format(db_name, table))
for x in gen:
    print(x)
