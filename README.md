# Pulsar timing-residual plot generator

`plot_timing_residuals_v4.py` is an interactive, publication-oriented plotting tool for TEMPO2 timing residuals. It was written for the software environment on the **tgss69 server** and remains compatible with Python 2.7.5 and older Matplotlib releases.

The script can process one pulsar or several pulsars in one run. For each input `.tim` file it finds the matching `.par` file, optionally filters to a selected MJD interval, asks TEMPO2 to refit that selected data, and writes PDF and 600-dpi PNG plots.

## What the program does

For every input pulsar, the script:

1. Requires a `.tim` file and a same-basename `.par` file in the same location.
2. Reads TOA MJDs from column 3 of ordinary TEMPO2-format data lines.
3. Lets the user select the complete baseline or an exact custom MJD interval.
4. Creates temporary `.tim` and `.par` files in the current directory. The temporary `.tim` contains only TOAs inside the chosen interval; the temporary `.par` omits existing `START` and `FINISH` lines.
5. Runs TEMPO2 with the `general2` output plugin and refits by default.
6. Extracts barycentric arrival MJD, post-fit residual, uncertainty, and observing frequency.
7. Converts MJDs to decimal calendar years and residuals to microseconds for plotting.
8. Lets the user choose figure dimensions, square-marker size, frequency colours/groups, legend visibility, and Y-axis limits.
9. Saves an individual PDF and PNG for each pulsar. With multiple inputs, it also saves a vertically stacked combined figure.

The physical X-axis limits remain exactly at the chosen MJD boundaries. Only complete integer calendar years are labelled; labels are thinned automatically for long baselines.

## tgss69 requirements

The intended tgss69 environment is:

- Linux shell
- Python 2.7.5 (Python 3 also works with the current source)
- NumPy
- Matplotlib with a non-interactive `Agg` backend
- TEMPO2, including the `general2` output plugin
- A writable working directory

Before running, load the normal tgss69 TEMPO2/Python environment used for pulsar timing work. Confirm that the executables and environment are visible in the same shell:

```bash
python --version
which tempo2
tempo2 -h
echo "$TEMPO2"
```

If tgss69 provides TEMPO2 through a login profile or module, activate that site-specific setup first. This repository deliberately does not guess a module name or hard-code a tgss69 installation path.

Confirm the Python packages:

```bash
python -c "import numpy, matplotlib; print(numpy.__version__); print(matplotlib.__version__)"
```

## Input-file rule

Every `.tim` file must have a corresponding `.par` file with exactly the same path and basename:

```text
J0248+4230_GWB.tim
J0248+4230_GWB.par
```

For example, passing `data/J0248+4230_GWB.tim` makes the script look for `data/J0248+4230_GWB.par`.

The parser expects the MJD in the third whitespace-separated field of a TOA line, consistent with the files for which this tool was designed. Non-TOA control/header lines and comments are retained in the temporary selected `.tim` file.

## Run on tgss69

Clone the repository and enter it:

```bash
git clone https://github.com/sapankumarsahoo/timing-residual-plot-generator.git
cd timing-residual-plot-generator
```

Copy or link the required `.tim`/`.par` pairs into a writable working directory. Run the script from that directory so temporary and output files are created there.

One pulsar:

```bash
python /path/to/timing-residual-plot-generator/plot_timing_residuals_v4.py \
    J0248+4230_GWB.tim
```

Two or more pulsars:

```bash
python /path/to/timing-residual-plot-generator/plot_timing_residuals_v4.py \
    J0248+4230_GWB.tim \
    J1207-5050_GWB.tim
```

To calculate residuals without refitting the timing model:

```bash
python /path/to/timing-residual-plot-generator/plot_timing_residuals_v4.py \
    --nofit J0248+4230_GWB.tim
```

Normally leave `--nofit` off.

You may also make the script executable and invoke it directly:

```bash
chmod +x plot_timing_residuals_v4.py
./plot_timing_residuals_v4.py J0248+4230_GWB.tim
```

## Interactive choices

The prompts appear in this order:

1. **MJD range for each pulsar** — complete range or a custom lower/upper MJD. Custom boundaries are kept exactly even when the first/last TOA lies inside them.
2. **Figure width and height** — defaults are 6.0 × 4.5 inches per panel.
3. **Square-marker size** — default is 4.0.
4. **Frequency-dependent colours** — automatic grouping (20 MHz default tolerance), manually specified bands, or one colour.
5. **Frequency legend** — keep or remove it.
6. **Y-axis range for each pulsar** — automatic, symmetric ± limit, or custom lower/upper limits in microseconds.

Press Enter to accept a displayed default.

Useful MJD checks already documented in the code:

```text
MJD 58849 = 2020-01-01
MJD 59215 = 2021-01-01
MJD 60310 = 2024-01-01
```

Thus a 2021-01-01 through 2024-01-01 plot uses 59215 through 60310.

## Outputs

For each pulsar, the program writes:

```text
PSR_<name>_timing_residuals.pdf
PSR_<name>_timing_residuals.png
```

When more than one pulsar is supplied, it also writes:

```text
timing_residuals_both.pdf
timing_residuals_both.png
```

Despite the historical `both` filename, the combined plot can contain more than two input pulsars. The PNG files are saved at 600 dpi; PDFs are vector output. Existing files with the same names are overwritten.

## Important implementation details

- TEMPO2 residuals are retained internally in seconds and plotted in microseconds.
- TEMPO2-reported uncertainties are interpreted as microseconds and converted to seconds internally.
- Post-fit RMS is computed after subtracting the mean residual and printed in microseconds.
- Automatic frequency grouping sorts frequencies and builds adjacent clusters within the selected tolerance of the current cluster mean.
- Data outside manually defined frequency bands are plotted in grey as `Other frequency`.
- Temporary files are named `selected_*.tim` and `selected_*.par` in the current directory and normally removed after TEMPO2 finishes.
- If the process is interrupted, leftover `selected_*` files can be removed after confirming that no run is active.

## Troubleshooting

### `ERROR: TEMPO2 not found`

TEMPO2 is not on `PATH` in the current shell. Activate the tgss69 timing environment and verify `which tempo2` before retrying.

### `ERROR: no residuals extracted`

Check that TEMPO2 succeeds for the same `.par`/`.tim` pair and that its `general2` plugin supports `{bat}`, `{post}`, `{err}`, and `{freq}`. Also inspect any TEMPO2 messages printed by the script.

### Corresponding PAR file not found

Rename or copy the parameter file so its path and basename exactly match the `.tim` input.

### No TOA MJDs could be read

The script expects numeric MJDs in field 3 and accepts values between 30000 and 100000. Check the `.tim` layout.

### Cannot create temporary files or plots

Run from a writable directory. The input files may live elsewhere, but the current directory must permit creation of temporary and output files.

### Font looks different on tgss69

The preferred serif fonts are Times New Roman, Times, and DejaVu Serif, in that order. Matplotlib falls back to an available serif font when needed.

## Notes and limitations

- The tool is interactive and has no non-interactive configuration-file mode.
- The X-axis uses a Gregorian calendar conversion from the MJD epoch and displays decimal years.
- The combined output filename is fixed.
- Panel labels are defined for the first four panels only.
- Broad exception handling is retained for compatibility with old Python/Matplotlib versions on tgss69.

## License

No license has been assigned. Add one before redistributing or accepting outside contributions.
