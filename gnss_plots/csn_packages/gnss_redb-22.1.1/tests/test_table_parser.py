# from datetime import datetime, timedelta, timezone
from time import perf_counter
from rethinkdb import RethinkDB
from gnss_redb.parser import GsofEcefParser, GsofEnuParser, PppEnuParser
from gnss_redb.coords import read_ref_coords

addresses = {
    'Atlas': ['10.54.217.15', 28015, [('collector', GsofEcefParser)]],
    # 'Bellaco': ['10.54.218.66', 28045],
    'Test': ['10.54.218.39', 28015, [('enu_data', GsofEnuParser),
                                     ('amqp_collector', PppEnuParser)]]
}
t_start = perf_counter()

ref_coords_dict = read_ref_coords()

for server, aux in addresses.items():
    host, port, dbs = aux
    re_db = RethinkDB()
    connection = re_db.connect(host=host, port=port)
    print('=====================================')
    print(re_db.db_list().run(connection))
    print('=====================================')
    for db, parser in dbs:
        tables = re_db.db(db).table_list().run(connection)
        print(' (+) DB: ', db)
        for table in tables:
            station = parser.table_to_station(table)
            try:
                assert table == parser.station_to_table(station)
                assert parser.is_valid_table_name(table)
            except:
                continue
            print(station, ref_coords_dict.get(station), table)
            cursor = re_db.db(db).table(table).limit(100).pluck(
                parser.relevant_fields).run(connection, time_format='raw')
            for doc in cursor:
                try:
                    dt = parser.parse_dt_gen(doc)
                    coords = parser.parse_coords(doc)
                    cov_coords = parser.parse_errors(doc)
                    try:
                        assert isinstance(coords, (tuple, list))
                        assert all(isinstance(x, float) for x in coords)
                        assert isinstance(cov_coords, (tuple, list))
                        assert all(isinstance(x, float) for x in cov_coords)
                        assert isinstance(dt, (int, float))
                    except AssertionError:
                        print('~~~~~~~~~~~~ AssertionError ~~~~~~~~~~~~')
                        print(coords)
                        print(cov_coords)
                        print(dt)
                except KeyError:
                    print('------------ KeyError ------------')
                    print(db, table)
                    print(doc)

print(' Execution time = {:.3f} s'.format(perf_counter()-t_start))
# dt_end = datetime.now(tz=timezone.utc,) - timedelta(seconds=6000)
# dt_start = dt_end - timedelta(seconds=600)
# start = dt_start.isoformat()
# end = dt_end.isoformat()
# dt_end = datetime.fromtimestamp(1602642312, tz=timezone.utc)
# dt_start = datetime.fromtimestamp(1602642308, tz=timezone.utc)
