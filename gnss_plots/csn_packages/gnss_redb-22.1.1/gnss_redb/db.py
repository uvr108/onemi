# -*- coding: utf-8 -*-
"""
*gnss_redb* is a package that allows to get data from a RethinkDB databases
with GNSS data. Designed for internal use in the CSN.

Copyright (c) 2020 Centro Sismológico Nacional, Universidad de Chile

This file is part of "gnss_redb".

"gnss_redb" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"gnss_redb" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "gnss_redb". If not, see <https://www.gnu.org/licenses/>.
"""
from typing import Union, Optional, Sequence, List, Set
from functools import wraps, partial
from itertools import chain
from datetime import datetime, timezone, timedelta
from dateutil.parser import parse
from asyncio import iscoroutinefunction
from rethinkdb import RethinkDB
from rethinkdb.errors import ReqlOpFailedError
import numpy as np
from geoproj.coords import ecef2geo
from gnss_redb.coords import StationsCoordsConverter, read_ref_coords
__docformat__ = 'reStructuredText en'

StrOrSequence = Optional[Union[str, Sequence]]


_time_strings = ("year", "month", "day", "hour", "minute")
_time_str_to_index = {t_str: (i+1) for i, t_str in enumerate(_time_strings)}


def datetime_window(end="now-hour", window_dict=None):
    """

    :param end: "now" or iso-8601 date-time string
    :param window_dict: time window length {"days": <d>, "hours": <h>, ...}
    :return: datetime objects at start and end of the window
    """
    if end[:3] == "now":
        dt_end = datetime.now(tz=timezone.utc)
        n = _time_str_to_index[end[4:]]
        args = tuple(getattr(dt_end, t_str) for t_str in _time_strings[:n])
        dt_end = datetime(*args, tzinfo=timezone.utc)
    else:
        dt_end = parse(end)
    if window_dict is None:
        window_dict = dict(days=1)
    dt_start = dt_end - timedelta(**window_dict)
    return dt_start, dt_end


def datetime_sequence(dt_start, dt_end, segment_dict=None):
    """

    :param dt_start: datetime of start
    :param segment_dict: time segment length {"days": <d>, "hours": <h>, ...}
    :return: list of datetime objects from start to the end of window
    """
    if segment_dict is None:
        segment_dict = dict(days=1)
    delta_segment = timedelta(**segment_dict)
    dt_list = [dt_start]
    dt_new = dt_start + delta_segment
    while dt_new < dt_end:
        dt_list.append(dt_new)
        dt_new += delta_segment
    dt_list.append(dt_end)
    return dt_list


def segmented_query(func):
    @wraps(func)
    def query_func(*args, segment_dict=None, **kwargs):
        if segment_dict is None:
            return func(*args, **kwargs)
        start, end = args[-2:]
        dt_list = datetime_sequence(start, end, segment_dict=segment_dict)

        def gen():
            for i in range(len(dt_list) - 1):
                for x in func(*chain(args[:-2], dt_list[i:(i+2)]), **kwargs):
                    yield x

        return gen()
    return query_func


def decorate_in_context(context_manager, func):

    if iscoroutinefunction(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            with context_manager(*args, **kwargs):
                return await func(*args, **kwargs)
    else:
        @wraps(func)
        def wrapper(*args, **kwargs):
            with context_manager(*args, **kwargs):
                return func(*args, **kwargs)

    return wrapper


class ServerDB(RethinkDB):

    def __init__(self, host, port, io_loop=None):
        super().__init__()
        self._host = host
        self._port = port
        self._io_loop = io_loop
        self._conn = None
        self._db_managers = dict()
        self._active_dbs = []
        self._station_coords_conv = StationsCoordsConverter()
        if self._io_loop is not None:
            self.set_loop_type(io_loop)

    # ----------------------------
    #      Coordinates ECEF, ENU
    def table_to_station(self, db, table: str) -> str:
        return self._db_managers[db].table_to_station(table)

    def station_to_table(self, db, station: str) -> str:
        return self._db_managers[db].station_to_table(station)

    def get_tables_set(self, db, update: bool = False) -> set:
        return self._db_managers[db].get_tables_set(update=update)

    def get_stations_set(self, db, update: bool = False) -> set:
        return self._db_managers[db].get_stations_set(update=update)

    def get_ref_llh(self, station):
        return self._station_coords_conv.get_ref_llh(station)

    def get_ref_ecef(self, station):
        return self._station_coords_conv.get_ref_ecef(station)

    def pos_cov_ecef_to_pos_std_enu(self, station, xyz, cov_xyz):
        return self._station_coords_conv.ecef_2_enu(
            station, xyz, cov_xyz)

    # ----------------------------
    #      connect
    @property
    def connection(self):
        return self._conn

    @property
    def host(self):
        return self._host

    @property
    def port(self):
        return self._port

    def do_connect(self, **kwargs):
        self._conn = self.connect(host=self._host, port=self._port, **kwargs)

    async def do_connect_async(self, **kwargs):
        self._conn = await self.connect(io_loop=self._io_loop,
                                        host=self._host, port=self._port,
                                        **kwargs)

    def reconnect(self, wait=False):
        self._conn.reconnect(noreply_wait=wait)

    async def reconnect_async(self, wait=False):
        await self._conn.reconnect(noreply_wait=wait)

    # ----------------------------
    #      DB managers
    def new_db_manager(self, db_name, table_parsers, update_table_set=True):
        self._active_dbs.append(db_name)
        self._db_managers[db_name] = DataBaseManager(
            db_name, self, self._conn, table_parsers,
            self._station_coords_conv)
        if update_table_set:
            self._db_managers[db_name].update_tables_set()

    def add_parser_to_db(self, db_name, parser, update_table_set=True):
        self._db_managers[db_name].add_parser(
            parser, update_table_set=update_table_set)

    def update_tables_set(self, db_name=None):
        if db_name is None:
            for db_manager in self._db_managers.values():
                db_manager.update_tables_set()
        else:
            self._db_managers[db_name].update_tables_set()

    update_table_set = update_tables_set

    def update_(self):
        for db_manager in self._db_managers.values():
            stations = db_manager.tabl
        # self._station_coords_conv

    # ----------------------------
    #      Info: DBs, Tables
    @property
    def active_dbs(self) -> List:
        return self._active_dbs

    def active_tables_in_db(self, db_name) -> Set:
        return self._db_managers[db_name].table_set

    def table_station_pairs(self, db_name, as_gen=True):
        return self._db_managers[db_name].table_station_pairs(as_gen=as_gen)

    # ----------------------------
    #      Get all DBs and tables
    def query_db_list(self):
        return self.db_list().run(self._conn)

    def query_table_list(self, db_name):
        return self.db(db_name).table_list().run(self._conn)

    # ----------------------------
    #      Queries: get data
    def count_interval(self, db_name: str, table: str,
                       start: datetime, end: datetime) -> int:
        return self._db_managers[db_name].count_interval(table, start, end)

    def time_interval_generator(self, db_name: str, table: str,
                                start: datetime, end: datetime,
                                only_dt_gen: bool = False,
                                segment_dict=None):
        return self._db_managers[db_name].time_interval_generator(
            table, start, end, only_dt_gen=only_dt_gen,
            segment_dict=segment_dict)

    def sample_values_generator(self, db_name: str, table, n_docs: int,
                                only_dt_gen: bool = False,):
        return self._db_managers[db_name].sample_values_generator(
            table, n_docs, only_dt_gen=only_dt_gen)


class DataBaseManager:

    def __init__(self, db_name: str, re_db: ServerDB,
                 connection, table_parsers: Sequence,
                 stations_coords_converter: StationsCoordsConverter):

        self._name = db_name
        self._se_db = re_db
        self._parsers_dict = {parser.name: parser
                              for parser in table_parsers}
        self._tables_by_parser = {parser.name: set()
                                  for parser in table_parsers}
        self._table_to_parser = dict()
        self._station_to_table = dict()
        self._all_tables = set()
        self._station_coords_conv = stations_coords_converter

    @property
    def db_name(self):
        return self._name

    @property
    def table_set(self):
        return self._all_tables

    def table_station_pairs(self, as_gen=True):

        def gen_func():
            for p_name, parser in self._parsers_dict.items():
                for table in self._tables_by_parser[p_name]:
                    yield (table, parser.table_to_station(table))
        if as_gen:
            return gen_func()
        return tuple(gen_func())

    def get_tables_set(self, update: bool = False) -> set:
        if update:
            self.update_tables_set()
        return self._all_tables

    def get_stations_set(self, update: bool = False) -> set:
        tables = self.get_tables_set(update=update)
        return set(self.table_to_station(table) for table in tables)

    def table_to_station(self, table: str) -> str:
        return self.get_parser(table).table_to_station(table)

    def station_to_table(self, station: str) -> str:
        return self._station_to_table.get(station)

    def update_tables_set(self):
        tables = self._se_db.query_table_list(self._name)
        self._aux_update_table_set(tables)
        ref_coords_dict = read_ref_coords()
        for table in tables:
            try:
                parser = self.get_parser(table)
                code = parser.table_to_station(table)
                rc_not_set = self._station_coords_conv.ref_coords_not_set(code)
                is_llh = True
                ref_coords = None
                if rc_not_set:
                    ref_coords = ref_coords_dict.get(code)
                    if (ref_coords is None) and parser.is_ecef:
                        x_values = []
                        y_values = []
                        z_values = []
                        for value in self.sample_values_generator(
                                table, 10, force_id_mapping=True):
                            x, y, z = value[1]
                            x_values.append(x)
                            y_values.append(y)
                            z_values.append(z)
                        if len(x_values) > 0:
                            ref_coords = (np.median(x_values),
                                          np.median(y_values),
                                          np.median(z_values))
                        else:
                            ref_coords = None
                        # compute
                        is_llh = False
                self._station_coords_conv.add_station(
                    code, ref_coords, is_llh=is_llh, compute_coeff=parser.is_ecef)
            except KeyError:
                pass

    update_table_set = update_tables_set

    async def update_table_set_async(self):
        tables = await self._se_db.query_table_list(self._name)
        self._aux_update_table_set(tables)

    def _aux_update_table_set(self, tables):
        for table in tables:
            if table in self._all_tables:
                continue
            for parser in self._parsers_dict.values():
                if parser.is_valid_table_name(table):
                    parser_name = parser.name
                    self._tables_by_parser[parser_name].add(table)
                    self._table_to_parser[table] = parser_name
                    self._all_tables.add(table)
                    break
        self._station_to_table = {}
        for table, station in self.table_station_pairs():
            self._station_to_table[station] = table

    def add_parser(self, parser, update_table_set=True):
        if parser.name in self._parsers_dict:
            return
        self._parsers_dict[parser.name] = parser
        self._tables_by_parser[parser.name] = set()
        if update_table_set:
            self.update_tables_set()

    def get_parser(self, table):
        return self._parsers_dict[self._table_to_parser.get(table)]

    def count_interval(self, table: str, start: datetime, end: datetime) -> int:
        parser, table_obj = self._aux_query(table)
        try:
            return table_obj.between(
                start, end, index=parser.dt_gen_field).count().run(
                self._se_db.connection, time_format='raw')
        except ReqlOpFailedError:
            return table_obj.between(start, end).count().run(
                self._se_db.connection, time_format='raw')

    @segmented_query
    def time_interval_generator(self, table: str,
                                start: datetime, end: datetime,
                                only_dt_gen: bool = False):
        parser, table_obj = self._aux_query(table)
        return self._aux_query_gen(parser, only_dt_gen, table,
                                   table_obj.between,
                                   start, end, index=parser.dt_gen_field)

    def sample_values_generator(self, table: str, n_docs: int,
                                only_dt_gen: bool = False,
                                force_id_mapping=False):
        parser, table_obj = self._aux_query(table)
        return self._aux_query_gen(parser, only_dt_gen, table,
                                   table_obj.limit, n_docs,
                                   force_id_mapping=force_id_mapping)

    def _aux_query_gen(self, parser, only_dt_gen, table, query,
                       *query_args, force_id_mapping=False, **query_kwargs):
        mapping = None
        if only_dt_gen:
            fds = parser.dt_gen_field
            gen_from_cursor = self._dt_gen_generator_from_cursor
        else:
            fds = parser.relevant_fields
            gen_from_cursor = self._values_generator_from_cursor
            if parser.is_ecef:
                if force_id_mapping:
                    mapping = None
                else:
                    mapping = partial(self._station_coords_conv.ecef_2_enu,
                                      parser.table_to_station(table))
        try:
            cursor = query(*query_args, **query_kwargs).pluck(fds).run(
                self._se_db.connection, time_format='raw')
        except ReqlOpFailedError:
            return []
        return gen_from_cursor(cursor, parser, mapping)

    def _aux_query(self, table):
        """Returns:  parser, table_object"""
        return self.get_parser(table), self._se_db.db(self._name).table(table)

    @classmethod
    def _values_generator_from_cursor(cls, cursor, parser, mapping):
        return cls._generator_from_parser(parser.parse_values, cursor, mapping)

    @classmethod
    def _dt_gen_generator_from_cursor(cls, cursor, parser, mapping):
        return cls._generator_from_parser(parser.parse_dt_gen, cursor, mapping)

    @staticmethod
    def _generator_from_parser(parse_func, cursor, mapping):

        if mapping is None:
            parse = parse_func
        else:
            def parse(doc):
                dt_gen, xyz, cov_xyz = parse_func(doc)
                enu, std_enu = mapping(xyz, cov_xyz)
                return dt_gen, enu, std_enu

        def gen_func():
            for doc in cursor:
                try:
                    yield parse(doc)
                except KeyError:
                    pass
        return gen_func()


def id_func(*x):
    return x
