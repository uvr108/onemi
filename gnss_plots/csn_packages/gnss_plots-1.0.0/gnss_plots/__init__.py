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
from typing import Optional, Union, Iterable, Sequence
from collections import Iterable as c_Iterable
from datetime import datetime, timezone, timedelta
from functools import wraps
import numpy as np
import matplotlib.pyplot as plt
from gnss_redb.db import ServerDB
from gnss_plots.datetime_ticks import set_ticks_on_axis, dt_ticks
__version__ = "1.0.0"

Time = Union[datetime, str, int, float]
_sta_label_color = '#0033b5'
_dt_label_color = '#005f96'



def update_current_db(method):
    @wraps(method)
    def new_method(instance, *args, db_name="", **kwargs):
        if len(db_name) > 0:
            instance.set_current_db(db_name)
        return method(instance, *args, db_name=instance.current_db, **kwargs)
    return new_method


class GnssPlot:

    LABELS = ("Este [m]", "Norte [m]", "Elev. [m]")
    ENU_COLORS = ("red", "blue", "black")
    DEFAULT_INTERVAL_SECS = 300

    def __init__(self, host, port, db_name: str = "",
                 parsers: Optional[Iterable] = None,
                 update_table_set: bool = True,
                 figsize: Sequence = (14, 8),
                 reference_time_at_end: bool = True):
        self._server_db = ServerDB(host, port)
        self._current_db = db_name
        if (len(db_name) > 0) and isinstance(parsers, c_Iterable):
            self._server_db.do_connect()
            self.set_current_db(db_name)
            self.new_db_manager(db_name, parsers,
                                update_table_set=update_table_set)

        self._fig = plt.Figure(figsize=figsize)
        self._ref_time_at_end = reference_time_at_end
        self._ref_datetime_label = ""

        self._ax_legend,self._axes = _add_axes(self._fig, n_plots=3)
        self._ax_legend.axis("off")
        self.legend_text = self._ax_legend.text(
            0.1, 0.5, '', ha='left', va='top', fontsize=15,
            color=_sta_label_color, family='',
            transform=self._ax_legend.transAxes)
        self.legend_time = self._ax_legend.text(
            0.3, 0.5, '', ha='left', va='top', fontsize=15,
            color=_dt_label_color, family='serif',
            transform=self._ax_legend.transAxes)
        self._y_label_kwargs = dict(fontsize=12,
                                    rotation="vertical",
                                    labelpad=10,
                                    va="center")

        self._dt_ticks_dict = None
        self._initial_format()

    def clear(self):
        for ax in self._axes:
            ax.clear()
        self.legend_text.set_text("")
        self.legend_time.set_text("")

    def do_connect(self):
        self._server_db.do_connect()

    def new_db_manager(self, db_name: str, parsers: Iterable,
                       update_table_set: bool = True):
            self._server_db.new_db_manager(
                db_name, parsers, update_table_set=update_table_set)

    def set_current_db(self, db_name: str):
        self._current_db = db_name

    @property
    def current_db(self):
        return self._current_db

    @update_current_db
    def get_tables_set(self, db_name: str = "", update: bool = False):
        return self._server_db.get_tables_set(db_name, update=update)

    @update_current_db
    def get_stations_set(self, db_name: str = "", update: bool = False):
        return self._server_db.get_stations_set(db_name, update=update)

    @update_current_db
    def time_interval_generator(self, table: str,
                                start: Optional[Time] = None,
                                end: Optional[Time] = None,
                                db_name: str = "",
                                only_dt_gen: bool = False,
                                segment_dict=None):

        if end is None:
            t_end = datetime.now(tz=timezone.utc)
        else:
            t_end = self._to_datetime_object(end)
        if start is None:
            t_start = t_end - timedelta(seconds=self.DEFAULT_INTERVAL_SECS)
        else:
            t_start = self._to_datetime_object(start)
        return self._server_db.time_interval_generator(
            db_name, table, t_start, t_end, only_dt_gen=only_dt_gen,
            segment_dict=segment_dict)

    def time_interval_arrays(self, table: str,
                             start: Optional[Time] = None,
                             end: Optional[Time] = None,
                             db_name: Optional[str] = "",
                             only_dt_gen: bool = False,
                             segment_dict=None):
        values = list(self.time_interval_generator(table, start, end,
                                                   db_name=db_name,
                                                   only_dt_gen=only_dt_gen,
                                                   segment_dict=segment_dict))
        if len(values) == 0:
            return None, (None, None, None), (None, None, None)
        n_points = len(values)
        if only_dt_gen:
            values.sort()
            return np.array(values)
        else:
            values.sort(key=lambda x: x[0])
            t_array = np.zeros(n_points)
            displacements = list(np.zeros(n_points) for _ in range(3))
            displacements_std = list(np.zeros(n_points) for _ in range(3))
            for i, val in enumerate(values):
                t_array[i] = val[0]
                for j in range(3):
                    displacements[j][i] = val[1][j]
                    displacements_std[j][i] = val[2][j]
            return t_array, displacements, displacements_std

    @update_current_db
    def create_single_plot(self, station: str, start: Optional[Time] = None,
                           end: Optional[Time] = None, db_name: str = "",
                           output_file="", dpi: int = 150):
        table = self._server_db.station_to_table(db_name, station)
        t_array, enu, std_enu = self.time_interval_arrays(
            table, start=start, end=end, db_name=db_name)
        if t_array is None:
            print(f"No hay datos en el intervalo consultado para la tabla {table}")
            return
        self._initial_format((t_array[0], t_array[-1]))
        self.set_axes_background_color()
        self.legend_text.set_text(station)
        for ax, x_array, color in zip(self._axes, enu, self.ENU_COLORS):
            ax.plot(t_array, x_array, color=color, linewidth=0.5)
        if isinstance(output_file, str) and (len(output_file) > 0):
            self.save_fig(output_file, dpi=dpi)

    def show_fig(self):
        self._fig.show()

    @staticmethod
    def _to_datetime_object(instance):
        if isinstance(instance, datetime):
            return instance
        elif isinstance(instance, (int, float)):
            return datetime.fromtimestamp(instance, tz=timezone.utc)
        elif isinstance(instance, str):
            return datetime.fromisoformat(instance)
        else:
            raise TypeError(
                "  Wrong type:  pass integer, float, string (iso8601) or datetime")

    @staticmethod
    def _to_utc_timestamp(instance):
        if isinstance(instance, (int, float)):
            return instance
        else:
            GnssPlot._to_datetime_object(instance).timestamp()

    def format_time_axes(self):
        """
        Format x-axis with ticks and tick labels in time units (minutes, hours,
        days, months or years) depending on .
        """
        set_ticks_on_axis(self._axes[-1], self._dt_ticks_dict)
        for k in range(2):
            self._axes[k].tick_params(labelbottom=False)
        self._axes[-1].set_xlabel("tiempo",
                                  fontsize=15,
                                  family='serif')

    def reset_t_lims(self):
        self._axes[-1].set_xlim(*self._dt_ticks_dict['t_interval'])

    def format_y_axes(self):
        for index, label in enumerate(self.LABELS):
            self._axes[index].set_ylabel(label, **self._y_label_kwargs)

    def _initial_format(self, ts_range = None):
        self.clear()
        self.set_t_ticks(ts_range)
        self.format_time_axes()
        self.format_y_axes()
        for ax in self._axes:
            # both axes: tick label size
            ax.tick_params(labelsize=10)

    @property
    def default_bg_color(self):
        return "#f4f8ff"

    def set_axes_background_color(self, color: str = "default"):
        bg_color = self.default_bg_color if (color == "default") else color
        for ax in self._axes:
            ax.set_facecolor(bg_color)

    def set_t_ticks(self, t_interval: Optional[Sequence] = None):
        if t_interval is None or not np.isfinite(t_interval[0]):
            t_final_utc = datetime.utcnow().timestamp()
            t_interval = (t_final_utc - 3590, t_final_utc)
        self._dt_ticks_dict = dt_ticks(*t_interval)
        self.set_datetime_label(t_interval)
        return t_interval

    def set_datetime_label(self, timestamp: Optional[Sequence] = None):
        if timestamp is not None:
            dt = datetime.utcfromtimestamp(
                timestamp[int(self._ref_time_at_end)])
            self._ref_datetime_label = dt.strftime(
                r'$t_{:s}$ :    %Y.%m.%d  %H:%M:%S   UTC'.format(
                    '{final}' if self._ref_time_at_end else '{begin}'))
        self.legend_time.set_text(self._ref_datetime_label)

    def save_fig(self, path: str, dpi: int = 300):
        self._fig.savefig(path, dpi=dpi)


def _add_axes(fig, **kwargs):
    # generates rectangles of legend and plot axes
    left = kwargs.get('padding_left', 0.06)
    right = kwargs.get('padding_right', 0.03)
    up = kwargs.get('padding_up', 0.05)
    down = kwargs.get('padding_down', 0.1)
    inter_axes_coeff = 1. + kwargs.get('inter_axes_rel', 0.08)
    height_weights = kwargs.get('height_weights')
    left_legend = kwargs.get('padding_legend_left', 0.03)
    right_legend = kwargs.get('padding_legend_right', 0.03)

    h_weights_is_none = (height_weights is None)
    n_plots = (kwargs.get('n_plots', 1) if h_weights_is_none
               else len(height_weights))
    h_free = (1 - up - down)/inter_axes_coeff

    if h_weights_is_none:
        dh = h_free/n_plots

        def get_height(index):
            return dh
    else:
        dh = h_free/sum(height_weights)

        def get_height(index):
            return dh*height_weights[index]

    # legend axes
    bottom = 1. - up
    ax_legend = fig.add_axes(
        (left_legend, bottom, 1. - (left_legend + right_legend), up))

    # rect's for the plot axes
    width = 1. - (left + right)
    axes = []
    axes_bottom = []
    axes_top = []

    for k in range(n_plots):
        height = get_height(k)
        bottom -= height*inter_axes_coeff
        axes_bottom.append(bottom)
        axes_top.append(bottom + height)
        if k > 0:
            axes.append(fig.add_axes((left, bottom, width, height),
                                     sharex=axes[0]))
        else:
            axes.append(fig.add_axes((left, bottom, width, height)))

    return ax_legend, axes