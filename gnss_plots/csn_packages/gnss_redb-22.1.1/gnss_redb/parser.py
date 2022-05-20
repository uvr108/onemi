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
import abc
from sys import modules
from inspect import getmembers, isclass
from typing import Union, Optional, Sequence
__docformat__ = 'reStructuredText en'
StrOrSequence = Optional[Union[str, Sequence]]


class TableParserInterface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'name') and
                hasattr(subclass, 'is_ecef') and
                hasattr(subclass, 'relevant_fields') and
                hasattr(subclass, 'dt_gen_field') and
                hasattr(subclass, 'table_to_station') and
                callable(subclass.table_to_station) and
                hasattr(subclass, 'station_to_table') and
                callable(subclass.station_to_table) and
                hasattr(subclass, 'is_valid_table_name') and
                callable(subclass.is_valid_table_name) and
                hasattr(subclass, 'parse_dt_gen') and
                callable(subclass.parse_dt_gen) and
                hasattr(subclass, 'parse_coords') and
                callable(subclass.parse_coords) and
                hasattr(subclass, 'parse_errors') and
                callable(subclass.parse_errors))

    @staticmethod
    @abc.abstractmethod
    def table_to_station(table: str) -> str:
        """Get station name"""
        raise NotImplementedError

    @staticmethod
    @abc.abstractmethod
    def station_to_table(station: str) -> str:
        """Get table name"""
        raise NotImplementedError

    @staticmethod
    @abc.abstractmethod
    def is_valid_table_name(station: str) -> bool:
        """Is this string a valid table name?"""
        raise NotImplementedError

    @staticmethod
    @abc.abstractmethod
    def parse_dt_gen(doc):
        """get epoch time of generation time"""
        raise NotImplementedError

    @staticmethod
    @abc.abstractmethod
    def parse_coords(doc):
        """get coordinates"""
        raise NotImplementedError

    @staticmethod
    @abc.abstractmethod
    def parse_errors(doc):
        """get errors of coordinates (covariances or standard deviations)"""
        raise NotImplementedError

    @classmethod
    def parse_values(cls, doc):
        """get (dt_gen, coords, errors)"""
        return (cls.parse_dt_gen(doc),
                cls.parse_coords(doc),
                cls.parse_errors(doc))


class GsofEcefParser(TableParserInterface):
    name = 'GSOF-ECEF'
    is_ecef = True  # ECEF coordinates
    relevant_fields = ['DT_GEN',
                       {'ECEF': ['X_POS', 'Y_POS', 'Z_POS'],
                        'POSITION_VCV': ['VCV_XX', 'VCV_XY', 'VCV_XZ',
                                         'VCV_YY', 'VCV_YZ', 'VCV_ZZ']}]
    dt_gen_field = relevant_fields[0]

    @staticmethod
    def table_to_station(table: str):
        return table[:-5]

    @staticmethod
    def station_to_table(station: str):
        return '{}_GSOF'.format(station)

    @staticmethod
    def is_valid_table_name(table: str):
        return table[-5:] == '_GSOF'

    @staticmethod
    def parse_dt_gen(doc):
        return doc['DT_GEN']['epoch_time']

    @staticmethod
    def parse_coords(doc):
        aux = doc['ECEF']
        return aux['X_POS'], aux['Y_POS'], aux['Z_POS']

    @staticmethod
    def parse_errors(doc):
        """Returns values of covariance matrix of ECEF coordinates

        :param doc:
        :return: C_xx, C_xy, C_xz, C_yy, C_yz, C_zz
        """
        aux = doc['POSITION_VCV']
        return (aux['VCV_XX'], aux['VCV_XY'], aux['VCV_XZ'],
                aux['VCV_YY'], aux['VCV_YZ'], aux['VCV_ZZ'])


class GsofEnuParser(TableParserInterface):
    name = 'GSOF-ENU'
    is_ecef = False  # ENU coordinates
    relevant_fields = ['DT_GEN', {'data': {'E': ['value', 'error'],
                                           'N': ['value', 'error'],
                                           'U': ['value', 'error']}}]
    dt_gen_field = relevant_fields[0]

    @staticmethod
    def table_to_station(table: str):
        return table[:-5]

    @staticmethod
    def station_to_table(station: str):
        return '{}_GSOF'.format(station)

    @staticmethod
    def is_valid_table_name(table: str):
        return table[-5:] == '_GSOF'

    @staticmethod
    def parse_dt_gen(doc):
        return doc['DT_GEN']['epoch_time']

    @staticmethod
    def parse_coords(doc):
        aux = doc['data']
        return aux['E']['value'], aux['N']['value'], aux['U']['value']

    @staticmethod
    def parse_errors(doc):
        aux = doc['data']
        return aux['E']['error'], aux['N']['error'], aux['U']['error']


class PppEnuParser(TableParserInterface):
    name = 'PPP-ENU'
    is_ecef = False  # ENU coordinates
    relevant_fields = [
        'DT_GEN', {'features': {'geometry': 'coordinates',
                                'properties': ['EError', 'NError', 'UError']}}]
    dt_gen_field = relevant_fields[0]

    @staticmethod
    def table_to_station(table: str):
        return table[:-4]

    @staticmethod
    def station_to_table(station: str):
        return '{}_PPP'.format(station)

    @staticmethod
    def is_valid_table_name(table: str):
        return table[-4:] == '_PPP'

    @staticmethod
    def parse_dt_gen(doc):
        return doc['DT_GEN']['epoch_time']

    @staticmethod
    def parse_coords(doc):
        return doc['features'][0]['geometry']['coordinates']

    @staticmethod
    def parse_errors(doc):
        aux = doc['features'][0]['properties']
        return aux['EError'], aux['NError'], aux['UError']


parser_classes = tuple(x[1] for x in getmembers(modules[__name__], isclass)
                       if (x[0][-6:] == 'Parser'))
parser_classes_dict = {parser.name: parser for parser in parser_classes}
