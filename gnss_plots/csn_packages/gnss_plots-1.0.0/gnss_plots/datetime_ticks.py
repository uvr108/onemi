# -*- coding: utf-8 -*-
"""
*gnss-plots* is a package that allows to get data from a RethinkDB databases
with GNSS data. Designed for internal use in the CSN.

Copyright (c) 2022 Centro Sismológico Nacional, Universidad de Chile

This file is part of "gnss-plots".

"gnss-plots" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"gnss-plots" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "gnss-plots". If not, see <https://www.gnu.org/licenses/>.
"""
from datetime import datetime, timezone
from matplotlib.ticker import AutoMinorLocator
import math as mt
import numpy as np
sec_min = 60  # 1 minute in seconds
sec_hour = 3600  # 1 hour in seconds
sec_day = 86400  # 1 day in seconds
sec_month = 2629800  # 1 month in seconds (30.4375 days)
sec_year = 31557600  # 1 year in seconds (365.25 days)
_time_intervals = [(5*sec_min,    'm',  6, '{:%H:%M}'),
                   (10*sec_min,   'm',  0, '{:%H:%M}'),
                   (20*sec_min,  '5m',  5, '{:%H:%M}'),
                   (45*sec_min, '10m',  5, '{:%H:%M}'),
                   (75*sec_min, '15m',  3, '{:%H:%M}'),
                   (3*sec_hour,   'h',  4, '{:%H:00}'),
                   (30*sec_hour,  'h',  0, '{:%H:00}'),
                   (7*sec_day,    'D',  4, '{:%d/%m}'),
                   (45*sec_day,   'D',  0, '{:%d/%m}'),
                   (6*sec_month,  'M',  3, '{:%Y.%m}'),
                   (18*sec_month, 'M',  0, '{:%Y.%m}'),
                   (5*sec_year,   'Y', 12, '{:%Y}'),
                   (20*sec_year,  'Y',  0, '{:%Y}'),
                   (np.inf,     '10Y',  0, '{:%Y}')
                   ]
_ts_aux = datetime(2016, 8, 21, 14, tzinfo=timezone.utc).timestamp()
default_t_interval = (_ts_aux, _ts_aux + 3*sec_hour)

_dict_months = {1: 'Ene', 2: 'Feb', 3: 'Mar', 4: 'Abr', 5: 'May', 6: 'Jun',
                7: 'Jul', 8: 'Ago', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}

_dict_units = {'': '', 's': 'seconds', 'm': 'minutes', 'h': 'hours',
               'D': 'days', 'M': 'months', 'Y': 'years'}
_not_rotated = ('', 'seconds', 'minutes', 'hours')
_rotation_angle = 45


def timestamp_to_isoformat(ts):
    return datetime.fromtimestamp(ts, tz=timezone.utc).isoformat()


def format_time_axis(ax, ts_interval):
    """
    Format x-axis with ticks and tick labels in time units (minutes, hours,
    days, months or years). See :py:func:`datetime_ticks.dt_ticks` for details.

    :param ax: Axes instance
    :param ts_interval: timestamp (unix time) interval
    :return: same as output of :py:func:`dt_ticks`
    :rtype: float, float, str, str
    """
    d_ticks = dt_ticks(*ts_interval)
    set_ticks_on_axis(ax, d_ticks)
    return d_ticks


def set_ticks_on_axis(ax, dict_ticks):
    ax.set_xticks(dict_ticks['major_ticks'])
    ax.set_xticklabels(
        dict_ticks['major_tick_labels'],
        rotation=None if dict_ticks['major_unit'] in _not_rotated
        else _rotation_angle)
    n_minor = dict_ticks['n_minor']
    if n_minor > 0:
        # ax.xaxis.set_minor_locator(FixedLocator(minor_ticks))
        ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor))
    ax.set_xlim(*dict_ticks['t_interval'])


def dt_ticks(ts_initial, ts_final):
    """
    Major and minor ticks and major tick labels around a time interval

          Interval                    Tick units
                                    Major     Minor
     0 min    to  5 min     -->      min    10 sec
     5 min    to 10 min     -->      min      --
    10 min    to 20 min     -->    5 min       min
    20 min    to 45 min     -->   10 min     2 min
    45 min    to 75 min     -->   15 min     5 min
    75 min    to  3 hours   -->      hour   15 min
     3 hours  to 30 hours   -->      hour     --
    30 hours  to  7 days    -->      day     6 hours
     7 days   to 45 days    -->      day      --
    45 days   to  6 months  -->      month     day
     6 months to 18 months  -->      month    --
    18 months to  5 years   -->      years     month
     5 years  to 20 years   -->      years    --
     > 20 years             -->   10 years    --

    :param ts_initial: initial timestamp (unix time)
    :param ts_final: final timestamp (unix time)
    :return: major ticks (unix time), major tick labels, number of minor ticks,
        major tick units, two strings with date and time at the beginning and
        end of the formatted interval (example: '2017-03-14 05:41:42').

    """
    major_ticks = []
    dt1 = datetime.fromtimestamp(ts_initial, timezone.utc)
    dt2 = datetime.fromtimestamp(ts_final-0.1, timezone.utc)
    yy1 = dt1.year
    yy2 = dt2.year
    mm1 = dt1.month
    mm2 = dt2.month
    dd1 = dt1.day
    dd2 = dt2.day
    h1 = dt1.hour
    h2 = dt2.hour
    m1 = dt1.minute
    m2 = dt2.minute

    major, n_minor, label_stencil = _major_minor(ts_final - ts_initial)

    # case: major ticks in minutes
    if major[-1] == 'm':
        int_min = 1 if len(major) == 1 else int(major[:-1])
        ts_1 = (_int_datetime(yy1, mm1, dd1, h1, m1)
                - (m1 % int_min)*sec_min)
        ts_2 = (_int_datetime(yy2, mm2, dd2, h2, m2) +
                (int_min - (m2 % int_min))*sec_min)
        dt_major = int_min*sec_min
        n_min = int(round((ts_2-ts_1)/dt_major))
        major_ticks = _major_ticks_regular(ts_1, n_min, dt_major)
    # case: major ticks in hours
    elif major == 'h':
        ts_1 = datetime(yy1, mm1, dd1, h1, tzinfo=timezone.utc).timestamp()
        ts_2 = datetime(yy2, mm2, dd2, h2, tzinfo=timezone.utc).timestamp()
        n_ticks = int(round((ts_2 - ts_1)/sec_hour)) + 1
        major_ticks = _major_ticks_regular(ts_1, n_ticks, sec_hour)
    # case: major ticks in days
    elif major == 'D':
        ts_1 = datetime(yy1, mm1, dd1, tzinfo=timezone.utc).timestamp()
        ts_2 = datetime(yy2, mm2, dd2, tzinfo=timezone.utc).timestamp()
        n_ticks = int(round((ts_2 - ts_1)/sec_day)) + 1
        major_ticks = _major_ticks_regular(ts_1, n_ticks, sec_day)
    # case: major ticks in months
    elif major == 'M':
        if yy1 == yy2:
            for mm in range(mm1, mm2+2):
                major_ticks.append(
                    datetime(yy1, mm, 1, tzinfo=timezone.utc).timestamp())
        else:
            for mm in range(mm1, 13):
                major_ticks.append(
                    datetime(yy1, mm, 1, tzinfo=timezone.utc).timestamp())
            for mm in range(1, mm2+2):
                major_ticks.append(
                    datetime(yy2, mm, 1, tzinfo=timezone.utc).timestamp())
    # case: major ticks in years
    elif major == 'Y':
        for yy in range(yy1, yy2+2):
            major_ticks.append(
                datetime(yy, 1, 1, tzinfo=timezone.utc).timestamp())
    # case: major == '10Y'
    else:
        for k in range(int(mt.floor(yy1/10)), int(mt.ceil(yy2/10)) + 1):
            yy = 10*k
            major_ticks.append(
                datetime(yy, 1, 1, tzinfo=timezone.utc).timestamp())
    major_tick_labels = [label_stencil.format(
        datetime.fromtimestamp(t, timezone.utc)) for t in major_ticks]
    return {
        'major_ticks': major_ticks, 'major_tick_labels': major_tick_labels,
        'n_minor': n_minor, 'major_unit': major,
        't_interval': (ts_initial, ts_final)}


def _major_minor(lapse):
    for tau, major, n_minor, label_stencil in _time_intervals:
        if lapse < tau:
            return major, n_minor, label_stencil


def _int_datetime(yy, mm, dd, h, m):
    return int(round(
        datetime(yy, mm, dd, h, m, tzinfo=timezone.utc).timestamp()))


def _major_ticks_regular(ts_init, n_major, deta_ts):
    major_ticks = [ts_init]
    ts_aux = ts_init
    for k in range(n_major):
        ts_aux += deta_ts
        major_ticks.append(ts_aux)
    return major_ticks


def _minor_ticks_regular(major_ticks, delta_ts_minor):
    minor_ticks = []
    ts_major_curr = major_ticks[0]
    for ts in major_ticks[1:]:
        ts_aux = ts_major_curr
        n_minor = int(round((ts - ts_major_curr)/delta_ts_minor)) - 1
        for k in range(n_minor):
            ts_aux += delta_ts_minor
            minor_ticks.append(ts_aux)
        ts_major_curr = ts
    return minor_ticks


def _minor_ticks_months(yy1, yy2):
    minor_ticks = []
    for yy in range(yy1, yy2 + 1):
        for mm in range(2, 12):
            minor_ticks.append(
                datetime(yy, mm, 1, tzinfo=timezone.utc).timestamp())
    return minor_ticks
