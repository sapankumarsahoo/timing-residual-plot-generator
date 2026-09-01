#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Publication-quality TEMPO2 timing residual plot.

Compatible with:
    Python 2.7.5
    Older Matplotlib
    TEMPO2

Important MJD check
-------------------
MJD 58849 = 2020-01-01.
MJD 59215 = 2021-01-01.
MJD 60310 = 2024-01-01.

Therefore, for a plot from 2021-01-01 to 2024-01-01,
use MJD 59215 to 60310.

Features
--------
1. Select custom MJD range.
2. TEMPO2 refits using only TOAs inside selected MJD range.
3. Plot boundaries remain EXACTLY at user-entered MJD limits.
4. X-axis displays ONLY integer calendar years.
5. Exact integer-year boundaries are shown inside the plotting
   box (left boundary left-aligned, right boundary right-aligned).
6. For long timing baselines, year labels are automatically
   thinned (every 2, 3, 4, or 5 years) to avoid overlap.
6. Control figure width and height.
6. Control square-marker size.
7. Frequency-dependent colours.
8. Optional frequency legend.
9. Interactive Y-axis limits entered in microseconds.
10. Y-axis is plotted directly in microseconds, with no
    multiplicative scientific-notation factor at the top.
11. The microsecond symbol uses a direct Unicode Greek mu
    for compatibility with Python 2.7 / old Matplotlib.
11. Saves PDF and 600-dpi PNG.

Example
-------

python plot_timing_residuals.py J0248+4230_GWB_sensitiveToAs.tim

or

python plot_timing_residuals.py \
    J0248+4230_GWB.tim \
    J1207-5050_GWB.tim
"""

from __future__ import print_function

import os
import sys
import argparse
import subprocess
import datetime
import tempfile

import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib.ticker import ScalarFormatter

try:
    from matplotlib.ticker import AutoMinorLocator
except:
    AutoMinorLocator = None

try:
    from matplotlib.ticker import MaxNLocator
except:
    MaxNLocator = None


# ============================================================
# PYTHON 2 / PYTHON 3 INPUT
# ============================================================

try:
    input_function = raw_input
except NameError:
    input_function = input


# ============================================================
# DEFAULT SETTINGS
# ============================================================

DEFAULT_MARKER_SIZE = 4.0

DEFAULT_FIGURE_WIDTH = 6.0

DEFAULT_FIGURE_HEIGHT = 4.5

ERROR_LINEWIDTH = 0.55

AXIS_LINEWIDTH = 0.8

ZERO_LINEWIDTH = 0.45


# ============================================================
# AUTOMATIC COLOURS
# ============================================================

AUTO_COLORS = [

    "black",

    "red",

    "blue",

    "green",

    "magenta",

    "orange",

    "cyan",

    "purple",

    "brown",

    "gray"
]


# ============================================================
# MATPLOTLIB STYLE
# ============================================================

plt.rcParams["font.family"] = "serif"


try:

    plt.rcParams["font.serif"] = [

        "Times New Roman",

        "Times",

        "DejaVu Serif"
    ]

except:

    pass


plt.rcParams["font.size"] = 12

plt.rcParams["axes.labelsize"] = 14

plt.rcParams["axes.titlesize"] = 15

plt.rcParams["xtick.labelsize"] = 12

plt.rcParams["ytick.labelsize"] = 12

plt.rcParams["axes.linewidth"] = AXIS_LINEWIDTH


try:

    plt.rcParams["xtick.direction"] = "in"

except:

    pass


try:

    plt.rcParams["ytick.direction"] = "in"

except:

    pass


try:

    plt.rcParams["pdf.fonttype"] = 42

except:

    pass


try:

    plt.rcParams["ps.fonttype"] = 42

except:

    pass


# ============================================================
# FIND CORRESPONDING PAR FILE
# ============================================================

def get_par_file(tim_file):


    basename = os.path.splitext(
        tim_file
    )[0]


    par_file = (
        basename +
        ".par"
    )


    if not os.path.isfile(
            par_file):


        print("")

        print("=" * 72)

        print("ERROR")

        print("=" * 72)

        print("")

        print(
            "Could not find corresponding PAR file:"
        )

        print(
            "    " + par_file
        )

        print("")


        sys.exit(1)


    return par_file


# ============================================================
# GET PULSAR NAME
# ============================================================

def get_psr_name(par_file):


    try:


        f = open(
            par_file,
            "r"
        )


        for line in f:


            line = line.strip()


            if len(line) == 0:

                continue


            if line.startswith("#"):

                continue


            parts = line.split()


            if len(parts) < 2:

                continue


            key = parts[0].upper()


            if (

                key == "PSRJ"

                or key == "PSR"

                or key == "PSRB"

            ):


                name = parts[1]


                f.close()


                return (
                    "PSR " +
                    name
                )


        f.close()


    except:

        pass


    # ========================================================
    # FALLBACK TO FILENAME
    # ========================================================


    name = os.path.basename(
        par_file
    )


    name = os.path.splitext(
        name
    )[0]


    name = name.replace(

        "_GWB_sensitiveToAs",

        ""
    )


    name = name.replace(

        "_sensitiveToAs",

        ""
    )


    name = name.replace(

        "_GWB",

        ""
    )


    return (
        "PSR " +
        name
    )


# ============================================================
# MJD -> DECIMAL YEAR
# ============================================================

def mjd_to_year(mjd_array):


    mjd_zero = datetime.datetime(

        1858,

        11,

        17
    )


    output = []


    for mjd in mjd_array:


        date = (

            mjd_zero

            +

            datetime.timedelta(

                days=float(mjd)

            )

        )


        year_start = datetime.datetime(

            date.year,

            1,

            1

        )


        next_year = datetime.datetime(

            date.year + 1,

            1,

            1

        )


        elapsed = (

            date -

            year_start

        ).total_seconds()


        duration = (

            next_year -

            year_start

        ).total_seconds()


        fraction = (

            elapsed /

            duration

        )


        output.append(

            date.year +

            fraction

        )


    return np.asarray(

        output,

        dtype=float

    )


# ============================================================
# SINGLE MJD -> DECIMAL YEAR
# ============================================================

def single_mjd_to_year(mjd):


    result = mjd_to_year(

        np.asarray(

            [mjd],

            dtype=float

        )

    )


    return float(

        result[0]

    )


# ============================================================
# SINGLE MJD -> CALENDAR DATE STRING
# ============================================================

def mjd_to_date_string(mjd):


    """
    Convert an MJD to YYYY-MM-DD.

    This is printed in the terminal so that the selected MJD
    range can be checked against the intended calendar years.
    """


    mjd_zero = datetime.datetime(

        1858,

        11,

        17

    )


    date = (

        mjd_zero

        +

        datetime.timedelta(

            days=float(mjd)

        )

    )


    return date.strftime(

        "%Y-%m-%d"

    )


# ============================================================
# READ MJDS FROM TIM FILE
# ============================================================

def read_tim_mjds(tim_file):


    mjds = []


    try:


        f = open(

            tim_file,

            "r"

        )


    except IOError:


        print("")

        print(
            "ERROR: Cannot open TIM file:"
        )

        print(
            "    " + tim_file
        )

        print("")


        sys.exit(1)


    for line in f:


        stripped = line.strip()


        if len(stripped) == 0:

            continue


        if stripped.startswith("#"):

            continue


        parts = stripped.split()


        if len(parts) < 4:

            continue


        try:


            candidate = float(

                parts[2]

            )


        except:


            continue


        if (

            candidate > 30000.0

            and

            candidate < 100000.0

        ):


            mjds.append(

                candidate

            )


    f.close()


    if len(mjds) == 0:


        print("")

        print("ERROR:")

        print(
            "No TOA MJDs could be read."
        )

        print("")


        sys.exit(1)


    return np.asarray(

        mjds,

        dtype=float

    )


# ============================================================
# SELECT MJD RANGE
# ============================================================

def choose_mjd_limits(
        psr_name,
        tim_file):


    mjds = read_tim_mjds(

        tim_file

    )


    available_min = np.min(

        mjds

    )


    available_max = np.max(

        mjds

    )


    print("")

    print("=" * 72)

    print(
        "MJD RANGE FOR " +
        psr_name
    )

    print("=" * 72)

    print("")


    print(
        "Available MJD range:"
    )


    print(

        "    %.6f -- %.6f"

        % (

            available_min,

            available_max

        )

    )


    print(
        "Available calendar-date range:"
    )


    print(

        "    %s -- %s"

        % (

            mjd_to_date_string(
                available_min
            ),

            mjd_to_date_string(
                available_max
            )

        )

    )


    print("")


    print(

        "Number of original TOAs: %d"

        % len(mjds)

    )


    print("")


    print(

        "    1 : Complete MJD range and refit"

    )


    print(

        "    2 : Custom MJD range and refit"

    )


    print("")


    while True:


        choice = input_function(

            "Select [1]: "

        ).strip()


        if choice == "":

            choice = "1"


        # ====================================================
        # COMPLETE RANGE
        # ====================================================

        if choice == "1":


            print("")


            print(
                "Requested plot range:"
            )


            print(

                "    %.6f -- %.6f"

                % (

                    available_min,

                    available_max

                )

            )


            print("")


            return (

                available_min,

                available_max

            )


        # ====================================================
        # CUSTOM RANGE
        # ====================================================

        elif choice == "2":


            while True:


                lower_text = input_function(

                    "Minimum MJD: "

                ).strip()


                upper_text = input_function(

                    "Maximum MJD: "

                ).strip()


                try:


                    requested_min = float(

                        lower_text

                    )


                    requested_max = float(

                        upper_text

                    )


                except:


                    print("")

                    print(
                        "Please enter valid MJD values."
                    )

                    print("")


                    continue


                if requested_min >= requested_max:


                    print("")

                    print(
                        "Minimum must be smaller than maximum."
                    )

                    print("")


                    continue


                # =================================================
                # AT LEAST SOME DATA MUST OVERLAP
                # =================================================

                if (

                    requested_max < available_min

                    or

                    requested_min > available_max

                ):


                    print("")

                    print(
                        "Selected interval contains no TOAs."
                    )

                    print("")


                    continue


                # =================================================
                # FIND ACTUAL TOAs WITHIN USER REQUEST
                # =================================================

                mask = (

                    (mjds >= requested_min)

                    &

                    (mjds <= requested_max)

                )


                number = np.sum(

                    mask

                )


                if number == 0:


                    print("")

                    print(
                        "No TOAs lie inside selected interval."
                    )

                    print("")


                    continue


                selected_toas = (

                    mjds[mask]

                )


                actual_min = np.min(

                    selected_toas

                )


                actual_max = np.max(

                    selected_toas

                )


                print("")

                print("=" * 72)

                print(
                    "SELECTED MJD RANGE"
                )

                print("=" * 72)

                print("")


                print(
                    "Requested plot range:"
                )


                print(

                    "    %.6f -- %.6f"

                    % (

                        requested_min,

                        requested_max

                    )

                )

                print(
                    "Requested calendar-date range:"
                )


                print(

                    "    %s -- %s"

                    % (

                        mjd_to_date_string(
                            requested_min
                        ),

                        mjd_to_date_string(
                            requested_max
                        )

                    )

                )


                print(
                    "Requested decimal-year range:"
                )


                print(

                    "    %.6f -- %.6f"

                    % (

                        single_mjd_to_year(
                            requested_min
                        ),

                        single_mjd_to_year(
                            requested_max
                        )

                    )

                )


                print("")


                print(
                    "Actual TOAs inside requested range:"
                )


                print(

                    "    %.6f -- %.6f"

                    % (

                        actual_min,

                        actual_max

                    )

                )


                print("")


                print(

                    "TOAs selected: %d"

                    % number

                )


                print("")


                print(
                    "Plot boundaries remain exactly:"
                )


                print(

                    "    %.6f -- %.6f"

                    % (

                        requested_min,

                        requested_max

                    )

                )


                print("")


                # =================================================
                # IMPORTANT:
                # RETURN EXACT USER VALUES.
                # DO NOT CLAMP TO ACTUAL TOAs.
                # =================================================

                return (

                    requested_min,

                    requested_max

                )


        else:


            print("")

            print(
                "Please select 1 or 2."
            )

            print("")


# ============================================================
# CREATE TEMPORARY TIM FILE
# ============================================================

def create_selected_tim_file(
        original_tim,
        mjd_min,
        mjd_max):


    temp_handle = tempfile.NamedTemporaryFile(

        mode="w",

        suffix=".tim",

        prefix="selected_",

        dir=".",

        delete=False

    )


    temp_tim = (

        temp_handle.name

    )


    f = open(

        original_tim,

        "r"

    )


    selected_count = 0


    for line in f:


        stripped = line.strip()


        # ====================================================
        # BLANK LINES
        # ====================================================

        if len(stripped) == 0:


            temp_handle.write(

                line

            )


            continue


        # ====================================================
        # COMMENTS
        # ====================================================

        if stripped.startswith("#"):


            temp_handle.write(

                line

            )


            continue


        parts = stripped.split()


        is_toa = False

        toa_mjd = None


        if len(parts) >= 4:


            try:


                candidate = float(

                    parts[2]

                )


                if (

                    candidate > 30000.0

                    and

                    candidate < 100000.0

                ):


                    is_toa = True

                    toa_mjd = (

                        candidate

                    )


            except:


                pass


        # ====================================================
        # CONTROL / HEADER LINE
        # ====================================================

        if not is_toa:


            temp_handle.write(

                line

            )


            continue


        # ====================================================
        # SELECTED TOA
        # ====================================================

        if (

            toa_mjd >= mjd_min

            and

            toa_mjd <= mjd_max

        ):


            temp_handle.write(

                line

            )


            selected_count += 1


    f.close()

    temp_handle.close()


    if selected_count == 0:


        try:


            os.remove(

                temp_tim

            )


        except:


            pass


        print("")

        print(
            "ERROR: no selected TOAs."
        )

        print("")


        sys.exit(1)


    print("")


    print(
        "Temporary selected TIM file:"
    )


    print(

        "    " +

        temp_tim

    )


    print(

        "Selected TOAs: %d"

        % selected_count

    )


    print("")


    return temp_tim


# ============================================================
# CREATE TEMPORARY PAR FILE
# ============================================================

def create_selected_par_file(
        original_par):


    temp_handle = tempfile.NamedTemporaryFile(

        mode="w",

        suffix=".par",

        prefix="selected_",

        dir=".",

        delete=False

    )


    temp_par = (

        temp_handle.name

    )


    f = open(

        original_par,

        "r"

    )


    for line in f:


        stripped = line.strip()


        parts = stripped.split()


        if len(parts) > 0:


            key = (

                parts[0].upper()

            )


            # =================================================
            # REMOVE ORIGINAL START / FINISH
            # =================================================

            if (

                key == "START"

                or

                key == "FINISH"

            ):


                continue


        temp_handle.write(

            line

        )


    f.close()

    temp_handle.close()


    return temp_par


# ============================================================
# CLEAN TEMP FILES
# ============================================================

def cleanup_temp_files(
        tim_file,
        par_file):


    try:


        if os.path.isfile(

                tim_file):


            os.remove(

                tim_file

            )


    except:


        pass


    try:


        if os.path.isfile(

                par_file):


            os.remove(

                par_file

            )


    except:


        pass


# ============================================================
# TEMPO2 REFIT
# ============================================================

def get_residuals(
        original_par,
        original_tim,
        mjd_limits,
        nofit=False):


    mjd_min = (

        mjd_limits[0]

    )


    mjd_max = (

        mjd_limits[1]

    )


    selected_tim = (

        create_selected_tim_file(

            original_tim,

            mjd_min,

            mjd_max

        )

    )


    selected_par = (

        create_selected_par_file(

            original_par

        )

    )


    # ========================================================
    # TEMPO2 OUTPUT FORMAT
    # ========================================================

    output_format = (

        "PYRES "

        "{bat} "

        "{post} "

        "{err} "

        "{freq}\\n"

    )


    command = [

        "tempo2"

    ]


    if nofit:


        command.append(

            "-nofit"

        )


    command.extend([

        "-output",

        "general2",

        "-s",

        output_format,

        "-f",

        os.path.abspath(

            selected_par

        ),

        os.path.abspath(

            selected_tim

        )

    ])


    print("")

    print("=" * 72)


    if nofit:


        print(

            "TEMPO2 SELECTED-RANGE CALCULATION"

        )


    else:


        print(

            "REFITTING SELECTED MJD RANGE"

        )


    print("=" * 72)

    print("")


    print(
        "Requested MJD interval:"
    )


    print(

        "    %.6f -- %.6f"

        % (

            mjd_min,

            mjd_max

        )

    )


    print("")


    try:


        process = subprocess.Popen(

            command,

            stdout=subprocess.PIPE,

            stderr=subprocess.PIPE

        )


    except OSError:


        cleanup_temp_files(

            selected_tim,

            selected_par

        )


        print("")

        print(
            "ERROR: TEMPO2 not found."
        )

        print("")


        sys.exit(1)


    stdout_data, stderr_data = (

        process.communicate()

    )


    if not isinstance(

            stdout_data,

            str):


        stdout_data = stdout_data.decode(

            "utf-8",

            "ignore"

        )


    if not isinstance(

            stderr_data,

            str):


        stderr_data = stderr_data.decode(

            "utf-8",

            "ignore"

        )


    if process.returncode != 0:


        cleanup_temp_files(

            selected_tim,

            selected_par

        )


        print("")

        print(

            stdout_data

        )


        print("")

        print(

            stderr_data

        )


        sys.exit(1)


    mjd = []

    residual = []

    error_us = []

    frequency = []


    # ========================================================
    # EXTRACT GENERAL2 DATA
    # ========================================================

    for line in stdout_data.splitlines():


        line = line.strip()


        if not line.startswith(

                "PYRES"):


            continue


        parts = line.split()


        if len(parts) < 5:


            continue


        try:


            this_mjd = float(

                parts[1]

            )


            this_residual = float(

                parts[2]

            )


            this_error = float(

                parts[3]

            )


            this_frequency = float(

                parts[4]

            )


        except:


            continue


        if not np.isfinite(

                this_mjd):


            continue


        if not np.isfinite(

                this_residual):


            continue


        if not np.isfinite(

                this_error):


            continue


        if not np.isfinite(

                this_frequency):


            continue


        mjd.append(

            this_mjd

        )


        residual.append(

            this_residual

        )


        error_us.append(

            this_error

        )


        frequency.append(

            this_frequency

        )


    # ========================================================
    # REMOVE TEMP FILES
    # ========================================================

    cleanup_temp_files(

        selected_tim,

        selected_par

    )


    if len(mjd) == 0:


        print("")

        print(
            "ERROR: no residuals extracted."
        )

        print("")


        sys.exit(1)


    # ========================================================
    # NUMPY ARRAYS
    # ========================================================

    mjd = np.asarray(

        mjd,

        dtype=float

    )


    residual = np.asarray(

        residual,

        dtype=float

    )


    error_sec = (

        np.asarray(

            error_us,

            dtype=float

        )

        *

        1.0e-6

    )


    frequency = np.asarray(

        frequency,

        dtype=float

    )


    # ========================================================
    # SORT BY MJD
    # ========================================================

    sort_index = np.argsort(

        mjd

    )


    mjd = mjd[

        sort_index

    ]


    residual = residual[

        sort_index

    ]


    error_sec = error_sec[

        sort_index

    ]


    frequency = frequency[

        sort_index

    ]


    # ========================================================
    # RMS
    # ========================================================

    mean_residual = np.mean(

        residual

    )


    rms = np.sqrt(

        np.mean(

            (

                residual -

                mean_residual

            ) ** 2

        )

    )


    print("")


    print(

        "Number of fitted TOAs : %d"

        % len(mjd)

    )


    print(

        "Actual TOA MJD range  : %.6f -- %.6f"

        % (

            np.min(mjd),

            np.max(mjd)

        )

    )


    print(

        "Post-fit RMS          : %.3f microseconds"

        % (

            rms *

            1.0e6

        )

    )


    print("")


    return (

        mjd,

        residual,

        error_sec,

        frequency,

        rms

    )


# ============================================================
# FIGURE SIZE
# ============================================================

def choose_figure_size():


    print("")

    print("=" * 72)

    print(
        "PLOT / FIGURE SIZE"
    )

    print("=" * 72)

    print("")


    width_text = input_function(

        "Figure width [6.0]: "

    ).strip()


    height_text = input_function(

        "Figure height [4.5]: "

    ).strip()


    try:


        width = (

            DEFAULT_FIGURE_WIDTH

            if width_text == ""

            else float(width_text)

        )


    except:


        width = (

            DEFAULT_FIGURE_WIDTH

        )


    try:


        height = (

            DEFAULT_FIGURE_HEIGHT

            if height_text == ""

            else float(height_text)

        )


    except:


        height = (

            DEFAULT_FIGURE_HEIGHT

        )


    if width <= 0:


        width = (

            DEFAULT_FIGURE_WIDTH

        )


    if height <= 0:


        height = (

            DEFAULT_FIGURE_HEIGHT

        )


    print("")


    print(
        "Using figure size:"
    )


    print(

        "    %.2f x %.2f inch"

        % (

            width,

            height

        )

    )


    print("")


    return (

        width,

        height

    )


# ============================================================
# MARKER SIZE
# ============================================================

def choose_marker_size():


    print("")

    print("=" * 72)

    print(
        "DATA POINT / SQUARE BOX SIZE"
    )

    print("=" * 72)

    print("")


    while True:


        text = input_function(

            "Square-marker size [4.0]: "

        ).strip()


        if text == "":


            return (

                DEFAULT_MARKER_SIZE

            )


        try:


            value = float(

                text

            )


        except:


            print("")

            print(
                "Enter a valid number."
            )

            print("")


            continue


        if value <= 0:


            print("")

            print(
                "Marker size must be > 0."
            )

            print("")


            continue


        return value


# ============================================================
# Y AXIS LIMITS
# ============================================================

def choose_y_limits(
        psr_name,
        residual,
        error):


    lower_us = (

        np.min(

            residual -

            error

        )

        *

        1.0e6

    )


    upper_us = (

        np.max(

            residual +

            error

        )

        *

        1.0e6

    )


    print("")

    print("=" * 72)

    print(

        "Y-AXIS RANGE FOR " +

        psr_name

    )

    print("=" * 72)

    print("")


    print(
        "Data + error-bar range:"
    )


    print(

        "    %.3f to %.3f microseconds"

        % (

            lower_us,

            upper_us

        )

    )


    print("")


    print(

        "    1 : Automatic"

    )


    print(

        "    2 : Symmetric +/- limit"

    )


    print(

        "    3 : Custom lower/upper limits"

    )


    print("")


    while True:


        choice = input_function(

            "Select [1]: "

        ).strip()


        if choice == "":


            choice = "1"


        if choice == "1":


            return None


        elif choice == "2":


            try:


                value = abs(

                    float(

                        input_function(

                            "Enter +/- limit in microseconds: "

                        )

                    )

                )


            except:


                continue


            if value <= 0:


                continue


            return (

                -value *

                1.0e-6,

                value *

                1.0e-6

            )


        elif choice == "3":


            try:


                low = float(

                    input_function(

                        "Lower limit in microseconds: "

                    )

                )


                high = float(

                    input_function(

                        "Upper limit in microseconds: "

                    )

                )


            except:


                continue


            if low >= high:


                continue


            return (

                low *

                1.0e-6,

                high *

                1.0e-6

            )


# ============================================================
# PRINT FREQUENCIES
# ============================================================

def print_available_frequencies(
        frequencies):


    freq = np.asarray(

        frequencies,

        dtype=float

    )


    unique = np.unique(

        np.round(

            freq,

            3

        )

    )


    print("")

    print(
        "Observed frequencies:"
    )

    print("")


    if len(unique) <= 40:


        for value in unique:


            print(

                "    %.3f MHz"

                % value

            )


    else:


        print(

            "    %.3f -- %.3f MHz"

            % (

                np.min(freq),

                np.max(freq)

            )

        )


    print("")


# ============================================================
# AUTOMATIC FREQUENCY GROUPING
# ============================================================

def make_automatic_frequency_groups(
        frequencies,
        tolerance):


    freq = np.sort(

        np.asarray(

            frequencies,

            dtype=float

        )

    )


    if len(freq) == 0:


        return []


    clusters = []


    current = [

        freq[0]

    ]


    for value in freq[1:]:


        current_mean = np.mean(

            current

        )


        if abs(

                value -

                current_mean

        ) <= tolerance:


            current.append(

                value

            )


        else:


            clusters.append(

                current

            )


            current = [

                value

            ]


    clusters.append(

        current

    )


    bands = []


    for i in range(

            len(clusters)):


        cluster = np.asarray(

            clusters[i],

            dtype=float

        )


        center = np.mean(

            cluster

        )


        bands.append(

            {

                "min":

                    np.min(cluster)

                    -

                    tolerance,


                "max":

                    np.max(cluster)

                    +

                    tolerance,


                "center":

                    center,


                "color":

                    AUTO_COLORS[

                        i %

                        len(AUTO_COLORS)

                    ],


                "label":

                    "%.0f MHz"

                    % center

            }

        )


    return bands


# ============================================================
# FREQUENCY COLOUR CHOICE
# ============================================================

def choose_frequency_bands(
        frequencies):


    print("")

    print("=" * 72)

    print(
        "FREQUENCY-DEPENDENT COLOURS"
    )

    print("=" * 72)


    print_available_frequencies(

        frequencies

    )


    print(

        "    1 : Automatic groups"

    )


    print(

        "    2 : Manual groups"

    )


    print(

        "    3 : Single colour"

    )


    print("")


    while True:


        choice = input_function(

            "Select [1]: "

        ).strip()


        if choice == "":


            choice = "1"


        # ====================================================
        # AUTOMATIC
        # ====================================================

        if choice == "1":


            text = input_function(

                "Grouping tolerance in MHz [20]: "

            ).strip()


            try:


                tolerance = (

                    20.0

                    if text == ""

                    else abs(

                        float(text)

                    )

                )


            except:


                tolerance = 20.0


            if tolerance <= 0:


                tolerance = 20.0


            return (

                make_automatic_frequency_groups(

                    frequencies,

                    tolerance

                )

            )


        # ====================================================
        # MANUAL
        # ====================================================

        elif choice == "2":


            try:


                number = int(

                    input_function(

                        "Number of frequency groups: "

                    )

                )


            except:


                continue


            if number <= 0:


                continue


            bands = []


            for i in range(

                    number):


                print("")


                print(

                    "Frequency group %d"

                    % (

                        i + 1

                    )

                )


                while True:


                    try:


                        low = float(

                            input_function(

                                "Minimum frequency (MHz): "

                            )

                        )


                        high = float(

                            input_function(

                                "Maximum frequency (MHz): "

                            )

                        )


                    except:


                        continue


                    if low < high:


                        break


                default_label = (

                    "%.0f-%.0f MHz"

                    % (

                        low,

                        high

                    )

                )


                label = input_function(

                    "Legend label [%s]: "

                    % default_label

                ).strip()


                if label == "":


                    label = (

                        default_label

                    )


                default_color = (

                    AUTO_COLORS[

                        i %

                        len(AUTO_COLORS)

                    ]

                )


                color = input_function(

                    "Colour [%s]: "

                    % default_color

                ).strip()


                if color == "":


                    color = (

                        default_color

                    )


                bands.append(

                    {

                        "min":

                            low,


                        "max":

                            high,


                        "center":

                            0.5 *

                            (

                                low +

                                high

                            ),


                        "color":

                            color,


                        "label":

                            label

                    }

                )


            return bands


        # ====================================================
        # SINGLE COLOUR
        # ====================================================

        elif choice == "3":


            color = input_function(

                "Colour [red]: "

            ).strip()


            if color == "":


                color = "red"


            return [

                {

                    "min":

                        -1.0e30,


                    "max":

                        1.0e30,


                    "center":

                        0.0,


                    "color":

                        color,


                    "label":

                        "Timing residuals"

                }

            ]


# ============================================================
# LEGEND ON/OFF
# ============================================================

def choose_show_labels():


    print("")

    print("=" * 72)

    print(
        "FREQUENCY LABELS / LEGEND"
    )

    print("=" * 72)

    print("")


    print(

        "    1 : Keep frequency legend"

    )


    print(

        "    2 : Remove frequency legend"

    )


    print("")


    while True:


        choice = input_function(

            "Select [1]: "

        ).strip()


        if choice == "":


            choice = "1"


        if choice == "1":


            return True


        elif choice == "2":


            return False


        else:


            print("")

            print(
                "Please select 1 or 2."
            )

            print("")


# ============================================================
# AUTOMATIC YEAR-LABEL SPACING
# ============================================================

def get_year_tick_step(
        plot_start_year,
        plot_end_year):


    """
    Choose a clean integer-year spacing so X-axis labels
    remain readable and do not overlap.

    Only the visible YEAR LABELS are thinned.
    The exact user-selected MJD plotting range is unchanged.

    Typical behaviour:
        <= 6 years   -> every 1 year
        <= 12 years  -> every 2 years
        <= 18 years  -> every 3 years
        <= 24 years  -> every 4 years
        > 24 years   -> every 5 years
    """


    span = (

        plot_end_year -

        plot_start_year

    )


    if span <= 6.0:


        return 1


    elif span <= 12.0:


        return 2


    elif span <= 18.0:


        return 3


    elif span <= 24.0:


        return 4


    else:


        return 5


# ============================================================
# FORMAT AXIS
# ============================================================

def format_axis(
        ax,
        data):


    residual = data[

        "residual"

    ]


    error = data[

        "error"

    ]


    y_limits = data[

        "y_limits"

    ]


    # ========================================================
    # Y AXIS IN MICROSECONDS
    #
    # Residuals are stored internally in seconds, but the
    # plotted Y-axis is expressed directly in microseconds.
    #
    # This removes scientific-notation factors such as
    # x10^-5 or x10^-4 from the top of the panel.
    # ========================================================

    residual_us = (

        residual *

        1.0e6

    )


    error_us = (

        error *

        1.0e6

    )


    yformatter = ScalarFormatter(

        useMathText=False

    )


    try:


        yformatter.set_scientific(

            False

        )


    except:


        pass


    try:


        yformatter.set_useOffset(

            False

        )


    except:


        pass


    ax.yaxis.set_major_formatter(

        yformatter

    )


    if MaxNLocator is not None:


        try:


            ax.yaxis.set_major_locator(

                MaxNLocator(5)

            )


        except:


            pass


    if AutoMinorLocator is not None:


        try:


            ax.yaxis.set_minor_locator(

                AutoMinorLocator(2)

            )


        except:


            pass


    # ========================================================
    # EXACT USER-SELECTED MJD PLOT BOUNDARIES
    #
    # IMPORTANT:
    #
    # The physical plot begins and ends exactly where
    # the user requested, even if no TOA exists there.
    # ========================================================

    plot_start_year = single_mjd_to_year(

        data[

            "selected_mjd_min"

        ]

    )


    plot_end_year = single_mjd_to_year(

        data[

            "selected_mjd_max"

        ]

    )


    print(
        "Plot X-axis calendar range: %s -- %s"
        % (
            mjd_to_date_string(
                data[
                    "selected_mjd_min"
                ]
            ),
            mjd_to_date_string(
                data[
                    "selected_mjd_max"
                ]
            )
        )
    )


    print(
        "Plot X-axis decimal-year range: %.6f -- %.6f"
        % (
            plot_start_year,
            plot_end_year
        )
    )


    ax.set_xlim(

        plot_start_year,

        plot_end_year

    )


    # ========================================================
    # X AXIS MAJOR TICKS
    #
    # CRITICAL CHANGE:
    #
    # Display ONLY COMPLETE INTEGER YEARS.
    #
    # Example:
    #
    # Plot range:
    #      2020.62 -------- 2024.15
    #
    # Labels:
    #
    #       2021   2022   2023   2024
    #
    # NO:
    #
    #       2020.62
    #
    # The box still begins at 2020.62.
    # ========================================================

    first_integer_year = int(

        np.ceil(

            plot_start_year

        )

    )


    last_integer_year = int(

        np.floor(

            plot_end_year

        )

    )


    major_ticks = []


    major_labels = []


    # ========================================================
    # AUTOMATIC YEAR-LABEL THINNING
    #
    # Keep only integer years, but for long time spans show
    # every 2nd / 3rd / 4th / 5th year as needed.
    #
    # The exact plot boundaries are NOT changed.
    # ========================================================

    year_step = get_year_tick_step(

        plot_start_year,

        plot_end_year

    )


    if first_integer_year <= last_integer_year:


        # ----------------------------------------------------
        # Align labels to a clean multiple of year_step.
        #
        # Example:
        #   step = 2 -> 2014, 2016, 2018, ...
        #   step = 3 -> 2013, 2016, 2019, ...
        # ----------------------------------------------------

        first_label_year = (

            int(

                np.ceil(

                    float(first_integer_year)

                    /

                    float(year_step)

                )

            )

            *

            year_step

        )


        for yy in range(

                first_label_year,

                last_integer_year + 1,

                year_step):


            yy_float = float(

                yy

            )


            if (

                yy_float >=

                plot_start_year

                and

                yy_float <=

                plot_end_year

            ):


                major_ticks.append(

                    yy_float

                )


                major_labels.append(

                    "%d"

                    % yy

                )


    # ========================================================
    # FALLBACK FOR VERY SHORT TIME RANGE
    # ========================================================

    if len(major_ticks) == 0:


        middle_year = (

            0.5 *

            (

                plot_start_year +

                plot_end_year

            )

        )


        rounded_middle = int(

            round(

                middle_year

            )

        )


        if (

            rounded_middle >=

            plot_start_year

            and

            rounded_middle <=

            plot_end_year

        ):


            major_ticks = [

                float(

                    rounded_middle

                )

            ]


            major_labels = [

                "%d"

                % rounded_middle

            ]


    # ========================================================
    # APPLY INTEGER YEAR TICKS
    # ========================================================

    if len(major_ticks) > 0:


        ax.set_xticks(

            major_ticks

        )


        ax.set_xticklabels(

            major_labels

        )


    # ========================================================
    # X AXIS MINOR TICKS
    #
    # Half-year minor ticks.
    #
    # These are NOT labelled.
    # ========================================================

    try:


        # For short spans keep half-year minor ticks.
        # For longer spans use one-year minor ticks.

        if (

            plot_end_year -

            plot_start_year

        ) <= 6.0:


            minor_step = 0.5


        else:


            minor_step = 1.0


        minor_start = (

            np.floor(

                plot_start_year /

                minor_step

            )

            *

            minor_step

        )


        minor_end = (

            np.ceil(

                plot_end_year /

                minor_step

            )

            *

            minor_step

        )


        minor_ticks = np.arange(

            minor_start,

            minor_end +

            0.5 * minor_step,

            minor_step

        )


        minor_ticks = minor_ticks[

            (

                minor_ticks >

                plot_start_year

            )

            &

            (

                minor_ticks <

                plot_end_year

            )

        ]


        # Remove minor ticks coinciding
        # with major integer-year ticks.

        clean_minor_ticks = []


        for tick in minor_ticks:


            is_major = False


            for major in major_ticks:


                if abs(

                    tick -

                    major

                ) < 0.01:


                    is_major = True

                    break


            if not is_major:


                clean_minor_ticks.append(

                    tick

                )


        ax.set_xticks(

            clean_minor_ticks,

            minor=True

        )


    except:


        pass


    # ========================================================
    # Y LIMITS IN MICROSECONDS
    # ========================================================

    if y_limits is None:


        lower = np.min(

            residual_us -

            error_us

        )


        upper = np.max(

            residual_us +

            error_us

        )


        maximum = max(

            abs(lower),

            abs(upper)

        )


        if maximum <= 0:


            maximum = (

                1.0

            )


        maximum = (

            maximum *

            1.10

        )


        ax.set_ylim(

            -maximum,

            maximum

        )


    else:


        # choose_y_limits() stores the selected limits in
        # seconds, so convert them here to microseconds.

        ax.set_ylim(

            y_limits[0] *

            1.0e6,

            y_limits[1] *

            1.0e6

        )


    # ========================================================
    # TICKS ON ALL FOUR SIDES
    # ========================================================

    try:


        ax.xaxis.set_ticks_position(

            "both"

        )


        ax.yaxis.set_ticks_position(

            "both"

        )


    except:


        pass


    # ========================================================
    # MAJOR TICKS
    # ========================================================

    try:


        ax.tick_params(

            axis="both",

            which="major",

            direction="in",

            length=6,

            width=0.8

        )


    except:


        pass


    # ========================================================
    # MINOR TICKS
    # ========================================================

    try:


        ax.tick_params(

            axis="both",

            which="minor",

            direction="in",

            length=3,

            width=0.6

        )


    except:


        pass


    # ========================================================
    # AXIS BOX
    # ========================================================

    try:


        ax.spines[

            "left"

        ].set_linewidth(

            AXIS_LINEWIDTH

        )


        ax.spines[

            "right"

        ].set_linewidth(

            AXIS_LINEWIDTH

        )


        ax.spines[

            "top"

        ].set_linewidth(

            AXIS_LINEWIDTH

        )


        ax.spines[

            "bottom"

        ].set_linewidth(

            AXIS_LINEWIDTH

        )


    except:


        pass


# ============================================================
# HIDE X-AXIS LABELS THAT CROSS THE PLOTTING BOX
# ============================================================

def hide_crossing_xlabels(
        fig,
        ax):


    """
    Keep exact integer-year boundary labels INSIDE the plotting
    box, and hide only any other X-axis label that still extends
    outside the box.

    Example for an exact 2020--2024 range:

        2020                                      2024
        |------------------------------------------|

    The 2020 label is left-aligned at the left boundary and the
    2024 label is right-aligned at the right boundary. Therefore
    neither label extends outside the plotting box.

    Nothing about the X-axis limits is changed.
    """


    try:


        # ----------------------------------------------------
        # CURRENT PHYSICAL X LIMITS
        # ----------------------------------------------------

        xlim = ax.get_xlim()

        xleft = float(
            xlim[0]
        )

        xright = float(
            xlim[1]
        )


        ticks = ax.get_xticks()

        labels = ax.get_xticklabels()


        # Tolerance in year units for identifying a tick that
        # is exactly at the selected boundary.

        boundary_tolerance = 1.0e-6


        # ----------------------------------------------------
        # FIRST: KEEP TRUE BOUNDARY YEAR LABELS INSIDE
        # ----------------------------------------------------

        for i in range(
                min(
                    len(ticks),
                    len(labels)
                )):


            tick = float(
                ticks[i]
            )

            label = labels[i]


            if label.get_text() == "":

                continue


            if abs(
                    tick -
                    xleft
            ) <= boundary_tolerance:


                # Text begins at the left edge and extends
                # inward to the right.

                label.set_horizontalalignment(
                    "left"
                )

                label.set_visible(
                    True
                )


            elif abs(
                    tick -
                    xright
            ) <= boundary_tolerance:


                # Text ends at the right edge and extends
                # inward to the left.

                label.set_horizontalalignment(
                    "right"
                )

                label.set_visible(
                    True
                )


            else:


                label.set_horizontalalignment(
                    "center"
                )


        # ----------------------------------------------------
        # DRAW ONCE TO GET REAL TEXT EXTENTS
        # ----------------------------------------------------

        fig.canvas.draw()


        renderer = (

            fig.canvas.get_renderer()

        )


        axes_bbox = (

            ax.get_window_extent(

                renderer

            )

        )


        tolerance_pixels = 1.0


        changed = False


        # ----------------------------------------------------
        # SECOND: HIDE ONLY NON-BOUNDARY LABELS THAT CROSS
        # ----------------------------------------------------

        for i in range(
                min(
                    len(ticks),
                    len(labels)
                )):


            tick = float(
                ticks[i]
            )

            label = labels[i]


            if not label.get_visible():

                continue


            if label.get_text() == "":

                continue


            # Exact left/right boundary labels were aligned
            # inward above, so preserve them.

            if (
                abs(
                    tick -
                    xleft
                ) <= boundary_tolerance
                or
                abs(
                    tick -
                    xright
                ) <= boundary_tolerance
            ):


                continue


            label_bbox = (

                label.get_window_extent(

                    renderer

                )

            )


            crosses_left = (

                label_bbox.x0 <

                axes_bbox.x0 -

                tolerance_pixels

            )


            crosses_right = (

                label_bbox.x1 >

                axes_bbox.x1 +

                tolerance_pixels

            )


            if (

                crosses_left

                or

                crosses_right

            ):


                label.set_visible(

                    False

                )


                changed = True


        if changed:


            fig.canvas.draw()


    except:


        # Keep compatibility with old Matplotlib.
        pass


# ============================================================
# FREQUENCY MASK
# ============================================================

def get_frequency_mask(
        frequency,
        band,
        last_band=False):


    if last_band:


        return (

            (frequency >= band["min"])

            &

            (frequency <= band["max"])

        )


    return (

        (frequency >= band["min"])

        &

        (frequency < band["max"])

    )


# ============================================================
# DRAW RESIDUAL PANEL
# ============================================================

def plot_residual_panel(
        ax,
        data,
        frequency_bands,
        marker_size,
        show_labels):


    year = data[

        "year"

    ]


    residual = data[

        "residual"

    ]


    error = data[

        "error"

    ]


    frequency = data[

        "frequency"

    ]


    # ========================================================
    # CONVERT RESIDUALS AND UNCERTAINTIES TO MICROSECONDS
    #
    # TEMPO2 residuals are retained internally in seconds for
    # calculations, but plotting is performed in microseconds.
    # ========================================================

    residual_plot = (

        residual *

        1.0e6

    )


    error_plot = (

        error *

        1.0e6

    )


    plotted = np.zeros(

        len(year),

        dtype=bool

    )


    # ========================================================
    # PLOT EACH FREQUENCY GROUP
    # ========================================================

    for i in range(

            len(frequency_bands)):


        band = frequency_bands[

            i

        ]


        mask = get_frequency_mask(

            frequency,

            band,

            i == (

                len(

                    frequency_bands

                ) - 1

            )

        )


        if np.sum(

                mask) == 0:


            continue


        plotted = (

            plotted |

            mask

        )


        ax.errorbar(

            year[

                mask

            ],

            residual_plot[

                mask

            ],

            yerr=error_plot[

                mask

            ],

            fmt="s",

            markersize=marker_size,

            markeredgewidth=0.0,

            color=band[

                "color"

            ],

            markerfacecolor=band[

                "color"

            ],

            markeredgecolor=band[

                "color"

            ],

            ecolor=band[

                "color"

            ],

            elinewidth=ERROR_LINEWIDTH,

            capsize=0,

            alpha=0.90,

            label=band[

                "label"

            ],

            zorder=3

        )


    # ========================================================
    # OTHER FREQUENCIES
    # ========================================================

    other_mask = (

        plotted == False

    )


    if np.sum(

            other_mask) > 0:


        ax.errorbar(

            year[

                other_mask

            ],

            residual_plot[

                other_mask

            ],

            yerr=error_plot[

                other_mask

            ],

            fmt="s",

            markersize=marker_size,

            markeredgewidth=0.0,

            color="gray",

            markerfacecolor="gray",

            markeredgecolor="gray",

            ecolor="gray",

            elinewidth=ERROR_LINEWIDTH,

            capsize=0,

            alpha=0.90,

            label="Other frequency",

            zorder=3

        )


    # ========================================================
    # ZERO RESIDUAL LINE
    # ========================================================

    ax.axhline(

        0.0,

        color="0.60",

        linestyle=":",

        linewidth=ZERO_LINEWIDTH,

        zorder=1

    )


    # ========================================================
    # AXIS LABEL
    # ========================================================

    ax.set_ylabel(

        u"Post-fit Residual (\u03bcs)"

    )


    # ========================================================
    # TITLE
    # ========================================================

    ax.set_title(

        data[

            "name"

        ]

    )


    try:


        ax.title.set_position(

            (

                0.5,

                1.02

            )

        )


    except:


        pass


    # ========================================================
    # AXIS FORMAT
    # ========================================================

    format_axis(

        ax,

        data

    )


    # ========================================================
    # OPTIONAL FREQUENCY LEGEND
    # ========================================================

    if show_labels:


        try:


            legend = ax.legend(

                loc="best",

                prop={

                    "size": 9

                },

                numpoints=1

            )


            try:


                legend.get_frame().set_linewidth(

                    0.6

                )


            except:


                pass


        except:


            pass


# ============================================================
# SAVE SINGLE PULSAR
# ============================================================

def save_single(
        data,
        frequency_bands,
        marker_size,
        figure_width,
        figure_height,
        show_labels):


    fig = plt.figure(

        figsize=(

            figure_width,

            figure_height

        )

    )


    ax = fig.add_subplot(

        111

    )


    plot_residual_panel(

        ax,

        data,

        frequency_bands,

        marker_size,

        show_labels

    )


    ax.set_xlabel(

        "Year"

    )


    fig.subplots_adjust(

        left=0.16,

        right=0.97,

        bottom=0.18,

        top=0.87

    )


    # ========================================================
    # REMOVE ONLY YEAR LABELS THAT CROSS THE AXES BOX
    # ========================================================

    hide_crossing_xlabels(

        fig,

        ax

    )


    clean_name = (

        data[

            "name"

        ]

        .replace(

            " ",

            "_"

        )

        .replace(

            "/",

            "_"

        )

    )


    pdf_name = (

        clean_name +

        "_timing_residuals.pdf"

    )


    png_name = (

        clean_name +

        "_timing_residuals.png"

    )


    # ========================================================
    # VECTOR PDF
    # ========================================================

    fig.savefig(

        pdf_name,

        bbox_inches="tight"

    )


    # ========================================================
    # 600 DPI PNG
    # ========================================================

    fig.savefig(

        png_name,

        dpi=600,

        bbox_inches="tight"

    )


    plt.close(

        fig

    )


    print("")

    print(
        "Saved:"
    )


    print(

        "    " +

        pdf_name

    )


    print(

        "    " +

        png_name

    )


    print("")


# ============================================================
# SAVE COMBINED FIGURE
# ============================================================

def save_combined(
        datasets,
        frequency_bands,
        marker_size,
        figure_width,
        figure_height,
        show_labels):


    number = len(

        datasets

    )


    fig = plt.figure(

        figsize=(

            figure_width,

            figure_height *

            number

        )

    )


    panel_labels = [

        "(a)",

        "(b)",

        "(c)",

        "(d)"

    ]


    for i in range(

            number):


        ax = fig.add_subplot(

            number,

            1,

            i + 1

        )


        data = datasets[

            i

        ]


        plot_residual_panel(

            ax,

            data,

            frequency_bands,

            marker_size,

            show_labels

        )


        ax.set_xlabel(

            "Year"

        )


        if i < len(

                panel_labels):


            try:


                ax.text(

                    0.018,

                    0.92,

                    panel_labels[

                        i

                    ],

                    transform=ax.transAxes,

                    horizontalalignment="left",

                    verticalalignment="top",

                    fontsize=12

                )


            except:


                pass


    fig.subplots_adjust(

        left=0.16,

        right=0.97,

        bottom=0.08,

        top=0.95,

        hspace=0.42

    )


    # ========================================================
    # REMOVE ONLY YEAR LABELS THAT CROSS EACH AXES BOX
    # ========================================================

    for ax in fig.axes:


        hide_crossing_xlabels(

            fig,

            ax

        )


    fig.savefig(

        "timing_residuals_both.pdf",

        bbox_inches="tight"

    )


    fig.savefig(

        "timing_residuals_both.png",

        dpi=600,

        bbox_inches="tight"

    )


    plt.close(

        fig

    )


    print("")

    print(
        "Saved combined figure:"
    )


    print(

        "    timing_residuals_both.pdf"

    )


    print(

        "    timing_residuals_both.png"

    )


    print("")


# ============================================================
# MAIN
# ============================================================

def main():


    parser = argparse.ArgumentParser(

        description=(

            "Selected-MJD TEMPO2 refit "

            "and publication timing-residual plot."

        )

    )


    parser.add_argument(

        "timfiles",

        nargs="+",

        help=(

            "Input TEMPO2 .tim file(s)"

        )

    )


    parser.add_argument(

        "--nofit",

        action="store_true",

        help=(

            "Do not refit timing model. "

            "Normally leave this option OFF."

        )

    )


    args = parser.parse_args()


    datasets = []


    # ========================================================
    # PROCESS EVERY INPUT PULSAR
    # ========================================================

    for tim_file in args.timfiles:


        if not os.path.isfile(

                tim_file):


            print("")

            print(
                "ERROR: TIM file not found:"
            )


            print(

                "    " +

                tim_file

            )


            print("")


            sys.exit(1)


        # ====================================================
        # PAR FILE
        # ====================================================

        par_file = get_par_file(

            tim_file

        )


        # ====================================================
        # PSR NAME
        # ====================================================

        psr_name = get_psr_name(

            par_file

        )


        # ====================================================
        # USER MJD RANGE
        # ====================================================

        mjd_limits = choose_mjd_limits(

            psr_name,

            tim_file

        )


        # ====================================================
        # PRESERVE EXACT REQUESTED LIMITS
        # ====================================================

        selected_mjd_min = (

            mjd_limits[0]

        )


        selected_mjd_max = (

            mjd_limits[1]

        )


        # ====================================================
        # REFIT SELECTED DATA
        # ====================================================

        (

            mjd,

            residual,

            error,

            frequency,

            rms

        ) = get_residuals(

            par_file,

            tim_file,

            mjd_limits,

            args.nofit

        )


        # ====================================================
        # STORE
        # ====================================================

        data = {


            "name":

                psr_name,


            "mjd":

                mjd,


            "year":

                mjd_to_year(

                    mjd

                ),


            "residual":

                residual,


            "error":

                error,


            "frequency":

                frequency,


            "selected_rms":

                rms,


            # ================================================
            # EXACT USER PLOT BOUNDARIES
            # ================================================

            "selected_mjd_min":

                selected_mjd_min,


            "selected_mjd_max":

                selected_mjd_max,


            "y_limits":

                None

        }


        datasets.append(

            data

        )


    # ========================================================
    # FIGURE SIZE
    # ========================================================

    (

        figure_width,

        figure_height

    ) = choose_figure_size()


    # ========================================================
    # MARKER SIZE
    # ========================================================

    marker_size = (

        choose_marker_size()

    )


    # ========================================================
    # ALL FREQUENCIES
    # ========================================================

    all_frequencies = np.concatenate(

        [

            data[

                "frequency"

            ]

            for data in datasets

        ]

    )


    # ========================================================
    # FREQUENCY COLOUR GROUPS
    # ========================================================

    frequency_bands = (

        choose_frequency_bands(

            all_frequencies

        )

    )


    # ========================================================
    # LEGEND
    # ========================================================

    show_labels = (

        choose_show_labels()

    )


    # ========================================================
    # Y LIMITS FOR EACH PULSAR
    # ========================================================

    for data in datasets:


        data[

            "y_limits"

        ] = choose_y_limits(


            data[

                "name"

            ],


            data[

                "residual"

            ],


            data[

                "error"

            ]

        )


    # ========================================================
    # SAVE INDIVIDUAL FIGURES
    # ========================================================

    for data in datasets:


        save_single(

            data,

            frequency_bands,

            marker_size,

            figure_width,

            figure_height,

            show_labels

        )


    # ========================================================
    # SAVE COMBINED FIGURE
    # ========================================================

    if len(

            datasets) > 1:


        save_combined(

            datasets,

            frequency_bands,

            marker_size,

            figure_width,

            figure_height,

            show_labels

        )


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":


    main()
