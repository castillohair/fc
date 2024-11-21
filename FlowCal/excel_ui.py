"""
``FlowCal``'s Microsoft Excel User Interface.

This module contains functions to read, gate, and transform data from a set
of FCS files, as specified by an input Microsoft Excel file. This file
should contain the following tables:
    - **Instruments**: Describes the instruments used to acquire the
      samples listed in the other tables. Each instrument is specified by a
      row containing at least the following fields:

      - **ID**: Short string identifying the instrument. Will be referenced
        by samples in the other tables.
      - **Forward Scatter Channel**: Name of the forward scatter channel,
        as specified by the ``$PnN`` keyword in the associated FCS files.
      - **Side Scatter Channel**: Name of the side scatter channel, as
        specified by the ``$PnN`` keyword in the associated FCS files.
      - **Fluorescence Channels**: Name of the fluorescence channels in a
        comma-separated list, as specified by the ``$PnN`` keyword in the
        associated FCS files.
      - **Time Channel**: Name of the time channel, as specified by the
        ``$PnN`` keyword in the associated FCS files.

    - **Beads** (optional): Describes the calibration beads samples that will
      be used to calibrate cell samples in the **Samples** table. The
      following information should be available for each beads sample:

      - **ID**: Short string identifying the beads sample. Will be
        referenced by cell samples in the **Samples** table.
      - **Instrument ID**: ID of the instrument used to acquire the sample.
        Must match one of the rows in the **Instruments** table.
      - **File Path**: Path of the FCS file containing the sample's data.
      - **<Fluorescence Channel Name> MEF Values**: The fluorescence in MEF
        of each bead subpopulation, as given by the manufacturer, as a
        comma-separated list of numbers. Any element of this list can be
        replaced with the word ``None``, in which case the corresponding
        subpopulation will not be used when fitting the beads fluorescence
        model. Note that the number of elements in this list (including
        the elements equal to ``None``) are the number of subpopulations
        that ``FlowCal`` will try to find.
      - **Gate fraction**: The fraction of events to keep from the sample
        after density-gating in the forward/side scatter channels.
      - **Clustering Channels**: The fluorescence channels used to identify
        the different bead subpopulations.

    - **Samples** (optional): Describes the biological samples to be
      processed. The following information should be available for each
      sample:

      - **ID**: Short string identifying the sample. Will be used as part
        of the plot's filenames and in the **Histograms** table in the
        output Excel file.
      - **Instrument ID**: ID of the instrument used to acquire the sample.
        Must match one of the rows in the **Instruments** table.
      - **Beads ID**: ID of the beads sample used to convert data to
        calibrated MEF.
      - **File Path**: Path of the FCS file containing the sample's data.
      - **<Fluorescence Channel Name> Units**: Units to which the event
        list in the specified fluorescence channel should be converted, and
        all the subsequent plots and statistics should be reported. Should
        be one of the following: "Channel" (raw units), "a.u." or "RFI"
        (arbitrary units) or "MEF" (calibrated Molecules of Equivalent
        Fluorophore). If "MEF" is specified, the **Beads ID** should be
        populated, and should correspond to a beads sample with the
        **MEF Values** specified for the same channel.
      - **Gate fraction**: The fraction of events to keep from the sample
        after density-gating in the forward/side scatter channels.

Any columns other than the ones specified above can be present, but will be
ignored by ``FlowCal``.

"""

import collections
import sys
import os
import os.path
import platform
import re
import six
import subprocess
import time
import warnings

# Tkinter is imported differently depending on which version of python we're
# using.
# six.PY2 is True when the code is running in python 2, False otherwise.
# six.PY3 is the equivalent for python 3.
if six.PY2:
    from Tkinter import Tk
    from tkFileDialog import askopenfilename
elif six.PY3:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
try:
    import xlrd
except ImportError:
    pass
import zipfile

import FlowCal.io
import FlowCal.plot
import FlowCal.gate
import FlowCal.transform
import FlowCal.stats
import FlowCal.mef

# Regular expressions for headers that specify some fluorescence channel
re_mef_values = re.compile(r'^\s*(\S(?:.*\S)?)\s+MEF\s+Values\s*$')
re_units      = re.compile(r'^\s*(\S(?:.*\S)?)\s+Units\s*$')

class ExcelUIException(Exception):
    """
    FlowCal Excel UI Error.

    """
    pass

def read_table(filename, sheetname, index_col=None, engine=None):
    """
    Return the contents of an Excel table as a pandas DataFrame.

    Parameters
    ----------
    filename : str
        Name of the Excel file to read.
    sheetname : str or int
        Name or index of the sheet inside the Excel file to read.
    index_col : str, optional
        Column name or index to be used as row labels of the DataFrame. If
        None, default index will be used.
    engine : str, optional
        Engine used by `pd.read_excel()` to read Excel file. If None, try
        'openpyxl' then 'xlrd'.

    Returns
    -------
    table : DataFrame
        A DataFrame containing the data in the specified Excel table. If
        `index_col` is not None, rows in which their `index_col` field
        is empty will not be present in `table`.

    Raises
    ------
    ValueError
        If `index_col` is specified and two rows contain the same
        `index_col` field.

    """
    # Catch sheetname as list or None
    if sheetname is None or \
            (hasattr(sheetname, '__iter__') \
            and not isinstance(sheetname, six.string_types)):
        raise TypeError("sheetname should specify a single sheet")

    # Load excel table using pandas. (pandas >= v1.2.0 delays closing files in
    # some situations, so open and close the file here.)
    with open(filename, 'rb') as f:
        file_in_mem = six.BytesIO(f.read())

    read_excel_kwargs = {'io'         : file_in_mem,
                         'sheet_name' : sheetname,
                         'index_col'  : index_col}

    if engine is None:
        # try reading Excel file using openpyxl engine first, then xlrd
        try:
            read_excel_kwargs['engine'] = 'openpyxl'
            table = pd.read_excel(**read_excel_kwargs)
        except ValueError as e:
            if not('openpyxl' in str(e).lower()
                   and 'unknown' in str(e).lower()):
                raise
            else:
                # pandas does not recognize openpyxl (e.g., pandas
                # version < 0.25.0), try xlrd
                read_excel_kwargs['engine'] = 'xlrd'
                table = pd.read_excel(**read_excel_kwargs)
        except ImportError:
            # pandas recognizes openpyxl but encountered an ImportError, try
            # xlrd. Possible scenarios: openpyxl version is less than what
            # pandas requires, openpyxl is missing (shouldn't happen)
            read_excel_kwargs['engine'] = 'xlrd'
            table = pd.read_excel(**read_excel_kwargs)
        except openpyxl.utils.exceptions.InvalidFileException:
            # unsupported file type (e.g., .xls), try xlrd
            #
            # (note: openpyxl's InvalidFileException has been stable at that
            # location since v2.2.0)
            read_excel_kwargs['engine'] = 'xlrd'
            table = pd.read_excel(**read_excel_kwargs)
        except zipfile.BadZipFile:
            # pandas >= 1.2.0 opens the file and passes the file buffer to
            # openpyxl, which may determine the file is not a zip file (e.g.,
            # if it's a .xls file), try xlrd
            read_excel_kwargs['engine'] = 'xlrd'
            table = pd.read_excel(**read_excel_kwargs)
    else:
        read_excel_kwargs['engine'] = engine
        table = pd.read_excel(**read_excel_kwargs)

    # Eliminate rows whose index are null
    if index_col is not None:
        table = table[pd.notnull(table.index)]
    # Check for duplicated rows
    if table.index.has_duplicates:
        raise ValueError("sheet {} on file {} contains duplicated values "
                         "for column {}".format(sheetname, filename, index_col))

    return table

def write_workbook(filename, table_list, column_width=None):
    """
    Write an Excel workbook from a list of tables.

    Parameters
    ----------
    filename : str
        Name of the Excel file to write.
    table_list : list of ``(str, DataFrame)`` tuples
        Tables to be saved as individual sheets in the Excel table. Each
        tuple contains two values: the name of the sheet to be saved as a
        string, and the contents of the table as a DataFrame.
    column_width: int or float, optional
        The column width to use when saving the spreadsheet. If None,
        calculate width automatically from the maximum number of characters
        in each column.

    """
    # Modify default header format
    # Pandas' default header format is bold text with thin borders. Here we
    # use bold text only, without borders.
    # Also, wrap in a try-except block in case style structure is not found.
    header_format_modified = False
    try:
        # Get format module
        import pandas.io.formats.excel as format_module

        # Save previous style, replace, and indicate that previous style should
        # be restored at the end
        old_header_style = format_module.header_style
        format_module.header_style = {"font": {"bold": True}}
        header_format_modified = True
    except AttributeError as e:
        pass

    # Generate output writer object
    writer = pd.ExcelWriter(filename, engine='openpyxl')

    # Write tables
    for sheet_name, df in table_list:
        # Convert index names to regular columns
        df = df.reset_index()
        # Write to an Excel sheet
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        # Set column width
        for i, (col_name, column) in enumerate(six.iteritems(df)):
            if column_width is None:
                # Get the maximum number of characters in a column
                max_chars_col = column.astype(str).str.len().max()
                max_chars_col = max(len(col_name), max_chars_col)
                width = float(max_chars_col)
            else:
                width = float(column_width)

            # Write width
            col_letter = openpyxl.utils.get_column_letter(i+1)
            writer.sheets[sheet_name].column_dimensions[col_letter].width = width

    # Save and close file
    writer.close()

    # Restore previous header format
    if header_format_modified:
        format_module.header_style = old_header_style

def process_beads_table(beads_table,
                        instruments_table,
                        base_dir=".",
                        verbose=False,
                        plot=False,
                        plot_dir=None,
                        full_output=False,
                        get_transform_fxn_kwargs={}):
    """
    Process calibration bead samples, as specified by an input table.

    This function processes the entries in `beads_table`. For each row, the
    function does the following:
        - Load the FCS file specified in the field "File Path".
        - Transform the forward scatter/side scatter and fluorescence
          channels to RFI
        - Remove the 250 first and 100 last events.
        - Remove saturated events in the forward scatter and side scatter
          channels.
        - Apply density gating on the forward scatter/side scatter
          channels.
        - Generate a standard curve transformation function, for each
          fluorescence channel in which the associated MEF values are
          specified.
        - Generate forward/side scatter density plots and fluorescence
          histograms, and plots of the clustering and fitting steps of
          standard curve generation, if `plot` = True.

    Names of forward/side scatter and fluorescence channels are taken from
    `instruments_table`.

    Parameters
    ----------
    beads_table : DataFrame
        Table specifying beads samples to be processed. For more
        information about the fields required in this table, please consult
        the module's documentation.
    instruments_table : DataFrame
        Table specifying instruments. For more information about the fields
        required in this table, please consult the module's documentation.
    base_dir : str, optional
        Directory from where all the other paths are specified.
    verbose : bool, optional
        Whether to print information messages during the execution of this
        function.
    plot : bool, optional
        Whether to generate and save density/histogram plots of each
        sample, and each beads sample.
    plot_dir : str, optional
        Directory relative to `base_dir` into which plots are saved. If
        `plot` is False, this parameter is ignored. If ``plot==True`` and
        ``plot_dir is None``, plot without saving.
    full_output : bool, optional
        Flag indicating whether to include an additional output, containing
        intermediate results from the generation of the MEF transformation
        functions.
    get_transform_fxn_kwargs : dict, optional
        Additional parameters passed directly to internal
        ``mef.get_transform_fxn()`` function call.

    Returns
    -------
    beads_samples : OrderedDict
        Processed, gated, and transformed samples, indexed by
        ``beads_table.index``.
    mef_transform_fxns : OrderedDict
        MEF transformation functions, indexed by ``beads_table.index``.
    mef_outputs : OrderedDict, only if ``full_output==True``
        Intermediate results from the generation of the MEF transformation
        functions. For every entry in `beads_table`,
        :func:`FlowCal.mef.get_transform_fxn()` is called on the
        corresponding processed and gated beads sample with
        ``full_output=True``, and the full output (a `MEFOutput`
        ``namedtuple``) is added to `mef_outputs`. `mef_outputs` is indexed
        by ``beads_table.index``. Refer to the documentation for
        :func:`FlowCal.mef.get_transform_fxn()` for more information.

    """
    # Initialize output variables
    beads_samples      = collections.OrderedDict()
    mef_transform_fxns = collections.OrderedDict()
    mef_outputs        = collections.OrderedDict()

    # Return empty structures if beads table is empty
    if beads_table.empty:
        if full_output:
            return beads_samples, mef_transform_fxns, mef_outputs
        else:
            return beads_samples, mef_transform_fxns

    if verbose:
        msg = "Processing Beads table ({} entries)".format(len(beads_table))
        print("")
        print(msg)
        print("="*len(msg))

    # Check that plotting directory exist, create otherwise
    if plot and plot_dir is not None \
            and not os.path.exists(os.path.join(base_dir, plot_dir)):
        os.makedirs(os.path.join(base_dir, plot_dir))

    # Extract header and channel names for which MEF values are specified.
    headers = list(beads_table.columns)
    mef_headers_all = [h for h in headers if re_mef_values.match(h)]
    mef_channels_all = [re_mef_values.match(h).group(1)
                        for h in mef_headers_all]

    # Iterate through table
    # We will look for a ExcelUIException on each iteration. If an exception
    # is caught, it will be stored in beads_samples.
    for beads_id, beads_row in beads_table.iterrows():
        try:
            ###
            # Instrument Data
            ###
            # Get the appropriate row in the instrument table
            instruments_row = instruments_table.loc[beads_row['Instrument ID']]
            # Scatter channels: Foward Scatter, Side Scatter
            sc_channels = [instruments_row['Forward Scatter Channel'],
                           instruments_row['Side Scatter Channel'],
                           ]
            # Fluorescence channels is a comma-separated list
            fl_channels = instruments_row['Fluorescence Channels'].split(',')
            fl_channels = [s.strip() for s in fl_channels]

            ###
            # Beads Data
            ###
            if verbose:
                print("\nBeads ID {}...".format(beads_id))
                print("Loading file \"{}\"...".format(beads_row['File Path']))

            # Attempt to open file
            filename = os.path.join(base_dir, beads_row['File Path'])
            try:
                beads_sample = FlowCal.io.FCSData(filename)
            except IOError:
                raise ExcelUIException("file \"{}\" not found".format(
                    beads_row['File Path']))
            # Check that the number of events is greater than 400
            if beads_sample.shape[0] < 400:
                raise ExcelUIException("number of events is lower than 400")

            ###
            # Transform
            ###
            if verbose:
                print("Performing data transformation...")
            # Transform FSC/SSC and fluorescence channels to linear scale
            beads_sample = FlowCal.transform.to_rfi(beads_sample,
                                                    sc_channels + fl_channels)

            # Parse clustering channels data
            cluster_channels = beads_row['Clustering Channels'].split(',')
            cluster_channels = [cc.strip() for cc in cluster_channels]

            ###
            # Gate
            ###
            if verbose:
                print("Performing gating...")
            # Remove first and last events. Transients in fluidics can make the
            # first few and last events slightly different from the rest.
            beads_sample_gated = FlowCal.gate.start_end(beads_sample,
                                                        num_start=250,
                                                        num_end=100)
            # Remove saturating events in forward/side scatter, if the FCS data
            # type is integer. The value of a saturating event is taken
            # automatically from `beads_sample_gated.range`.
            if beads_sample_gated.data_type == 'I':
                beads_sample_gated = FlowCal.gate.high_low(
                    beads_sample_gated,
                    channels=sc_channels)
            # Density gating
            try:
                density_gate_output = FlowCal.gate.density2d(
                    data=beads_sample_gated,
                    channels=sc_channels,
                    gate_fraction=beads_row['Gate Fraction'],
                    xscale='logicle',
                    yscale='logicle',
                    sigma=5.,
                    full_output=True)
                beads_sample_gated = density_gate_output.gated_data
                gate_contour       = density_gate_output.contour
            except ValueError as ve:
                raise ExcelUIException(ve.message)

            # Plot forward/side scatter density plot and fluorescence histograms
            if plot:
                if verbose:
                    print("Plotting density plot and histogram...")
                # Density plot parameters
                density_params = {}
                density_params['mode'] = 'scatter'
                density_params["title"] = "{} ({:.1f}% retained)".format(
                    beads_id,
                    beads_sample_gated.shape[0] * 100. / beads_sample.shape[0])
                density_params['xscale'] = 'logicle'
                density_params['yscale'] = 'logicle'
                # Beads have a tight distribution, so axis limits will be set
                # from 0.75 decades below the 5th percentile to 0.75 decades
                # above the 95th percentile.
                density_params['xlim'] = \
                    (np.percentile(beads_sample_gated[:, sc_channels[0]],
                                   5) / (10**0.75),
                     np.percentile(beads_sample_gated[:, sc_channels[0]],
                                   95) * (10**0.75),
                     )
                density_params['ylim'] = \
                    (np.percentile(beads_sample_gated[:, sc_channels[1]],
                                   5) / (10**0.75),
                     np.percentile(beads_sample_gated[:, sc_channels[1]],
                                   95) * (10**0.75),
                     )
                # Beads have a tight distribution, so less smoothing should be
                # applied for visualization
                density_params['sigma'] = 5.
                # Histogram plot parameters
                hist_params = {'xscale': 'logicle'}
                # Plot
                if plot_dir is not None:
                    figname = os.path.join(
                        base_dir,
                        plot_dir,
                        "density_hist_{}.png".format(beads_id))
                else:
                    figname = None
                plt.figure(figsize=(6,4))
                FlowCal.plot.density_and_hist(
                    beads_sample,
                    beads_sample_gated,
                    density_channels=sc_channels,
                    hist_channels=cluster_channels,
                    gate_contour=gate_contour,
                    density_params=density_params,
                    hist_params=hist_params,
                    savefig=figname)

            ###
            # Process MEF values
            ###
            # For each fluorescence channel, check whether a list of known MEF
            # values of the bead subpopulations is provided in `beads_row`. This
            # involves checking that a column named "[channel] MEF Values"
            # exists and is not empty. If so, store the name of the channel in
            # `mef_channels`, and the specified MEF values in `mef_values`.
            ###
            mef_values = []
            mef_channels = []
            for fl_channel in fl_channels:
                if fl_channel in mef_channels_all:
                    # Get header from channel name
                    mef_header = \
                        mef_headers_all[mef_channels_all.index(fl_channel)]
                    # Extract text. If empty, ignore.
                    mef_str = beads_row[mef_header]
                    if pd.isnull(mef_str):
                        continue
                    # Save channel name
                    mef_channels.append(fl_channel)
                    # Parse list of values
                    mef = mef_str.split(',')
                    mef = [int(e) if e.strip().isdigit() else np.nan
                           for e in mef]
                    mef_values.append(mef)

            # Ensure matching number of `mef_values` for all channels (this
            # implies that the calibration beads have the same number of
            # subpopulations for all channels).
            if mef_values:
                if not np.all([len(mef_values_channel)==len(mef_values[0])
                               for mef_values_channel in mef_values]):
                    raise ExcelUIException("Must specify the same number of"
                                           + " MEF Values for each channel."
                                           + " Use 'None' to instruct FlowCal"
                                           + " to ignore a detected"
                                           + " subpopulation.")
            mef_values = np.array(mef_values)

            # Obtain standard curve transformation
            if mef_channels:
                if verbose:
                    if len(mef_channels) == 1:
                        print("Calculating standard curve for channel {}..." \
                            .format(mef_channels[0]))
                    else:
                        print("Calculating standard curve for channels {}..." \
                            .format(", ".join(mef_channels)))

                mef_output = FlowCal.mef.get_transform_fxn(
                    beads_sample_gated,
                    mef_values,
                    mef_channels=mef_channels,
                    clustering_channels=cluster_channels,
                    verbose=False,
                    plot=plot,
                    plot_filename=beads_id,
                    plot_dir=os.path.join(base_dir, plot_dir) \
                             if plot_dir is not None else None,
                    full_output=full_output,
                    **get_transform_fxn_kwargs)

                if full_output:
                    mef_transform_fxn = mef_output.transform_fxn
                else:
                    mef_transform_fxn = mef_output

            else:
                mef_transform_fxn = None
                mef_output = None

        except ExcelUIException as e:
            # Print Exception message
            if verbose:
                print("ERROR: {}".format(str(e)))
            # Add exception to beads_samples dictionary, and None to
            # everything else
            beads_samples[beads_id] = e
            mef_transform_fxns[beads_id] = None
            if full_output:
                mef_outputs[beads_id] = None

        else:
            # If no errors were found, store results
            beads_samples[beads_id] = beads_sample_gated
            mef_transform_fxns[beads_id] = mef_transform_fxn
            if full_output:
                mef_outputs[beads_id] = mef_output


    if full_output:
        return beads_samples, mef_transform_fxns, mef_outputs
    else:
        return beads_samples, mef_transform_fxns

def process_samples_table(samples_table,
                          instruments_table,
                          mef_transform_fxns=None,
                          beads_table=None,
                          base_dir=".",
                          verbose=False,
                          plot=False,
                          plot_dir=None):
    """
    Process flow cytometry samples, as specified by an input table.

    The function processes each entry in `samples_table`, and does the
    following:
        - Load the FCS file specified in the field "File Path".
        - Transform the forward scatter/side scatter to RFI.
        - Transform the fluorescence channels to the units specified in the
          column "<Channel name> Units".
        - Remove the 250 first and 100 last events.
        - Remove saturated events in the forward scatter and side scatter
          channels.
        - Apply density gating on the forward scatter/side scatter
          channels.
        - Plot combined forward/side scatter density plots and fluorescence
          historgrams, if `plot` = True.

    Names of forward/side scatter and fluorescence channels are taken from
    `instruments_table`.

    Parameters
    ----------
    samples_table : DataFrame
        Table specifying samples to be processed. For more information
        about the fields required in this table, please consult the
        module's documentation.
    instruments_table : DataFrame
        Table specifying instruments. For more information about the fields
        required in this table, please consult the module's documentation.
    mef_transform_fxns : dict or OrderedDict, optional
        Dictionary containing MEF transformation functions. If any entry
        in `samples_table` requires transformation to MEF, a key: value
        pair must exist in mef_transform_fxns, with the key being equal to
        the contents of field "Beads ID".
    beads_table : DataFrame, optional
        Table specifying beads samples used to generate
        `mef_transform_fxns`. This is used to check if a beads sample was
        taken at the same acquisition settings as a sample to be
        transformed to MEF. For any beads sample and channel for which a
        MEF transformation function has been generated, the following
        fields should be populated: ``<channel> Amp. Type`` and
        ``<channel> Detector Volt``. If `beads_table` is not specified, no
        checking will be performed.
    base_dir : str, optional
        Directory from where all the other paths are specified.
    verbose : bool, optional
        Whether to print information messages during the execution of this
        function.
    plot : bool, optional
        Whether to generate and save density/histogram plots of each
        sample, and each beads sample.
    plot_dir : str, optional
        Directory relative to `base_dir` into which plots are saved. If
        `plot` is False, this parameter is ignored. If ``plot==True`` and
        ``plot_dir is None``, plot without saving.

    Returns
    -------
    samples : OrderedDict
        Processed, gated, and transformed samples, indexed by
        ``samples_table.index``.

    """
    # Initialize output variable
    samples = collections.OrderedDict()

    # Return empty dictionary if samples table is empty
    if samples_table.empty:
        return samples

    if verbose:
        msg = "Processing Samples table ({} entries)".format(len(samples_table))
        print("")
        print(msg)
        print("="*len(msg))

    # Check that plotting directory exist, create otherwise
    if plot and plot_dir is not None \
            and not os.path.exists(os.path.join(base_dir, plot_dir)):
        os.makedirs(os.path.join(base_dir, plot_dir))

    # Extract header and channel names for which units are specified.
    headers = list(samples_table.columns)
    report_headers_all = [h for h in headers if re_units.match(h)]
    report_channels_all = [re_units.match(h).group(1)
                           for h in report_headers_all]

    # Iterate through table
    # We will look for a ExcelUIException on each iteration. If an exception
    # is caught, it will be stored in beads_samples.
    for sample_id, sample_row in samples_table.iterrows():
        try:
            ###
            # Instrument Data
            ###
            # Get the appropriate row in the instrument table
            instruments_row = instruments_table.loc[sample_row['Instrument ID']]
            # Scatter channels: Foward Scatter, Side Scatter
            sc_channels = [instruments_row['Forward Scatter Channel'],
                           instruments_row['Side Scatter Channel'],
                           ]
            # Fluorescence channels is a comma-separated list
            fl_channels = instruments_row['Fluorescence Channels'].split(',')
            fl_channels = [s.strip() for s in fl_channels]

            ###
            # Sample Data
            ###
            if verbose:
                print("\nSample ID {}...".format(sample_id))
                print("Loading file \"{}\"...".format(sample_row['File Path']))

            # Attempt to open file
            filename = os.path.join(base_dir, sample_row['File Path'])
            try:
                sample = FlowCal.io.FCSData(filename)
            except IOError:
                raise ExcelUIException("file \"{}\" not found".format(
                    sample_row['File Path']))
            # Check that the number of events is greater than 400
            if sample.shape[0] < 400:
                raise ExcelUIException("number of events is lower than 400")

            ###
            # Transform
            ###
            if verbose:
                print("Performing data transformation...")
            # Transform FSC/SSC to linear scale
            sample = FlowCal.transform.to_rfi(sample, sc_channels)

            # Parse fluorescence channels in which to transform
            report_channels = []
            report_units = []
            for fl_channel in fl_channels:
                if fl_channel in report_channels_all:
                    # Get header from channel name
                    report_header = report_headers_all[
                        report_channels_all.index(fl_channel)]
                    # Extract text. If empty, ignore.
                    units_str = sample_row[report_header]
                    if pd.isnull(units_str):
                        continue
                    # Decide what transformation to perform
                    units = units_str.strip()
                    if units.lower() == 'channel':
                        units_label = "Channel Number"

                    elif units.lower() == 'rfi':
                        units_label = "Relative Fluorescence Intensity, RFI"
                        sample = FlowCal.transform.to_rfi(sample, fl_channel)
                    elif units.lower() == 'a.u.' or units.lower() == 'au':
                        units_label = "Arbitrary Units, a.u."
                        sample = FlowCal.transform.to_rfi(sample, fl_channel)

                    elif units.lower() == 'mef':
                        units_label = "Molecules of Equivalent Fluorophore, MEF"
                        # Check if transformation function is available
                        if mef_transform_fxns is None or \
                                mef_transform_fxns[sample_row['Beads ID']] \
                                is None:
                            raise ExcelUIException("MEF transformation "
                                "function not available")

                        # If beads_table is available, check if the same
                        # settings have been used to acquire the corresponding
                        # beads sample
                        if beads_table is not None:
                            beads_row = beads_table.loc[sample_row['Beads ID']]
                            # Instrument
                            beads_iid = beads_row['Instrument ID']
                            if beads_iid != sample_row['Instrument ID']:
                                raise ExcelUIException("Instruments for "
                                    "acquisition of beads and samples are not "
                                    "the same (beads {}'s instrument: {}, "
                                    "sample's instrument: {})".format(
                                        sample_row['Beads ID'],
                                        beads_iid,
                                        sample_row['Instrument ID']))
                            # Amplification type
                            beads_at = beads_row['{} Amp. Type'. \
                                format(fl_channel)]
                            if sample.amplification_type(fl_channel)[0]:
                                sample_at = "Log"
                            else:
                                sample_at = "Linear"
                            if beads_at != sample_at:
                                raise ExcelUIException("Amplification type for "
                                    "acquisition of beads and samples in "
                                    "channel {} are not the same (beads {}'s "
                                    "amplification: {}, sample's "
                                    "amplification: {})".format(
                                        fl_channel,
                                        sample_row['Beads ID'],
                                        beads_at,
                                        sample_at))
                            # Detector voltage
                            beads_dv = beads_row['{} Detector Volt.'. \
                                format(fl_channel)]
                            if sample.detector_voltage(fl_channel) is not None \
                                    and beads_dv != sample.detector_voltage(
                                        fl_channel):
                                raise ExcelUIException("Detector voltage for "
                                    "acquisition of beads and samples in "
                                    "channel {} are not the same (beads {}'s "
                                    "detector voltage: {}, sample's "
                                    "detector voltage: {})".format(
                                        fl_channel,
                                        sample_row['Beads ID'],
                                        beads_dv,
                                        sample.detector_voltage(fl_channel)))

                        # First, transform to RFI
                        sample = FlowCal.transform.to_rfi(sample, fl_channel)
                        # Attempt to transform to MEF
                        # Transformation function raises a ValueError if a
                        # standard curve does not exist for a channel
                        try:
                            sample = mef_transform_fxns[sample_row['Beads ID']](
                                sample,
                                fl_channel)
                        except ValueError:
                            raise ExcelUIException("no standard curve for "
                                "channel {}".format(fl_channel))
                    else:
                        raise ExcelUIException("units \"{}\" not recognized". \
                            format(units, sample_id))

                    # Register that reporting in this channel must be done
                    report_channels.append(fl_channel)
                    report_units.append(units_label)

            ###
            # Gate
            ###
            if verbose:
                print("Performing gating...")
            # Remove first and last events. Transients in fluidics can make the
            # first few and last events slightly different from the rest.
            sample_gated = FlowCal.gate.start_end(sample,
                                                  num_start=250,
                                                  num_end=100)
            # Remove saturating events in forward/side scatter, and fluorescent
            # channels to report, if the FCS data type is integer. The value of
            # a saturating event is taken automatically from
            # `sample_gated.range`.
            if sample_gated.data_type == 'I':
                sample_gated = FlowCal.gate.high_low(
                    sample_gated,
                    sc_channels + report_channels)
            # Density gating
            try:
                density_gate_output = FlowCal.gate.density2d(
                    data=sample_gated,
                    channels=sc_channels,
                    gate_fraction=sample_row['Gate Fraction'],
                    xscale='logicle',
                    yscale='logicle',
                    full_output=True)
                sample_gated = density_gate_output.gated_data
                gate_contour = density_gate_output.contour
            except ValueError as ve:
                raise ExcelUIException(ve.message)

            # Plot forward/side scatter density plot and fluorescence histograms
            if plot:
                if verbose:
                    print("Plotting density plot and histogram...")
                # Density plot parameters
                density_params = {}
                density_params['mode'] = 'scatter'
                density_params["title"] = "{} ({:.1f}% retained)".format(
                    sample_id,
                    sample_gated.shape[0] * 100. / sample.shape[0])
                density_params['xscale'] = 'logicle'
                density_params['yscale'] = 'logicle'
                # Histogram plot parameters
                hist_params = []
                for rc, ru in zip(report_channels, report_units):
                    param = {}
                    param['xlabel'] = '{} ({})'.format(rc, ru)
                    # Only channel numbers are plotted in linear scale
                    if (ru == 'Channel Number'):
                        param['xscale'] = 'linear'
                    else:
                        param['xscale'] = 'logicle'
                    hist_params.append(param)

                # Plot
                if plot_dir is not None:
                    figname = os.path.join(
                        base_dir,
                        plot_dir,
                        "{}.png".format(sample_id))
                else:
                    figname = None
                FlowCal.plot.density_and_hist(
                    sample,
                    sample_gated,
                    gate_contour=gate_contour,
                    density_channels=sc_channels,
                    density_params=density_params,
                    hist_channels=report_channels,
                    hist_params=hist_params,
                    savefig=figname)

        except ExcelUIException as e:
            # Print Exception message
            if verbose:
                print("ERROR: {}".format(str(e)))
            # Add exception to samples dictionary
            samples[sample_id] = e

        else:
            # If no errors were found, store results
            samples[sample_id] = sample_gated

    return samples

def add_beads_stats(beads_table, beads_samples, mef_outputs=None):
    """
    Add stats fields to beads table.

    The following information is added to each row:

        - Notes (warnings, errors) resulting from the analysis
        - Number of Events
        - Acquisition Time (s)

    The following information is added for each row, for each channel in
    which MEF values have been specified:

        - Detector voltage (gain)
        - Amplification type
        - Bead model fitted parameters

    Parameters
    ----------
    beads_table : DataFrame
        Table specifying bead samples to analyze. For more information
        about the fields required in this table, please consult the
        module's documentation.
    beads_samples : dict or OrderedDict
        FCSData objects from which to calculate statistics.
        ``beads_samples[id]`` should correspond to ``beads_table.loc[id,:]``.
    mef_outputs : dict or OrderedDict, optional
        Intermediate results from the generation of the MEF transformation
        functions, as given by ``mef.get_transform_fxn()``. This is used to
        populate the fields ``<channel> Beads Model``,
        ``<channel> Beads Params. Names``, and
        ``<channel> Beads Params. Values``. If specified,
        ``mef_outputs[id]`` should correspond to ``beads_table.loc[id,:]``.

    """
    # The index name is not preserved if beads_table is empty.
    # Save the index name for later
    beads_table_index_name = beads_table.index.name

    # Add per-row info
    notes = []
    n_events = []
    acq_time = []
    for row_id in beads_table.index:
        # Check if sample is an exception, otherwise assume it's an FCSData
        if isinstance(beads_samples[row_id], ExcelUIException):
            # Print error message
            notes.append("ERROR: {}".format(str(beads_samples[row_id])))
            n_events.append(np.nan)
            acq_time.append(np.nan)
        else:
            notes.append('')
            n_events.append(beads_samples[row_id].shape[0])
            acq_time.append(beads_samples[row_id].acquisition_time)

    beads_table['Analysis Notes'] = notes
    beads_table['Number of Events'] = n_events
    beads_table['Acquisition Time (s)'] = acq_time

    # List of channels that require stats columns
    headers = list(beads_table.columns)
    stats_headers = [h for h in headers if re_mef_values.match(h)]
    stats_channels = [re_mef_values.match(h).group(1) for h in stats_headers]

    # Iterate through channels
    for header, channel in zip(stats_headers, stats_channels):
        # Add empty columns to table
        beads_table[channel + ' Detector Volt.'] = np.nan
        beads_table[channel + ' Amp. Type'] = ""
        if mef_outputs:
            beads_table[channel + ' Beads Model'] = ""
            beads_table[channel + ' Beads Params. Names'] = ""
            beads_table[channel + ' Beads Params. Values'] = ""

        # Iterate
        for row_id in beads_table.index:
            # If error, skip
            if isinstance(beads_samples[row_id], ExcelUIException):
                continue
            # If MEF values are specified, calculate stats. If not, leave empty.
            if pd.notnull(beads_table[header][row_id]):

                # Detector voltage
                beads_table.at[row_id, channel + ' Detector Volt.'] = \
                    beads_samples[row_id].detector_voltage(channel)


                # Amplification type
                if beads_samples[row_id].amplification_type(channel)[0]:
                    amplification_type = "Log"
                else:
                    amplification_type = "Linear"
                beads_table.at[row_id, channel + ' Amp. Type'] = \
                    amplification_type


                # Bead model and parameters
                # Only populate if mef_outputs has been provided
                if mef_outputs:
                    # Try to find the current channel among the mef'd channels.
                    # If successful, extract bead fitted parameters.
                    try:
                        mef_channel_index = mef_outputs[row_id]. \
                            mef_channels.index(channel)
                    except ValueError:
                        pass
                    else:
                        # Bead model
                        beads_model_str = mef_outputs[row_id]. \
                            fitting['beads_model_str'][mef_channel_index]
                        beads_table.at[row_id, channel + ' Beads Model'] = \
                            beads_model_str

                        # Bead parameter names
                        params_names = mef_outputs[row_id]. \
                            fitting['beads_params_names'][mef_channel_index]
                        params_names_str = ", ".join([str(p)
                                                      for p in params_names])
                        beads_table.at[
                            row_id,
                            channel + ' Beads Params. Names'] = \
                                params_names_str

                        # Bead parameter values
                        params = mef_outputs[row_id]. \
                            fitting['beads_params'][mef_channel_index]
                        params_str = ", ".join([str(p) for p in params])
                        beads_table.at[
                            row_id,
                            channel + ' Beads Params. Values'] = \
                                params_str

    # Restore index name if table is empty
    if len(beads_table) == 0:
        beads_table.index.name = beads_table_index_name

def add_samples_stats(samples_table, samples):
    """
    Add stats fields to samples table.

    The following information is added to each row:

        - Notes (warnings, errors) resulting from the analysis
        - Number of Events
        - Acquisition Time (s)

    The following information is added for each row, for each channel in
    which fluorescence units have been specified:

        - Detector voltage (gain)
        - Amplification type
        - Mean
        - Geometric Mean
        - Median
        - Mode
        - Standard Deviation
        - Coefficient of Variation (CV)
        - Geometric Standard Deviation
        - Geometric Coefficient of Variation
        - Inter-Quartile Range
        - Robust Coefficient of Variation (RCV)

    Parameters
    ----------
    samples_table : DataFrame
        Table specifying samples to analyze. For more information about the
        fields required in this table, please consult the module's
        documentation.
    samples : dict or OrderedDict
        FCSData objects from which to calculate statistics. ``samples[id]``
        should correspond to ``samples_table.loc[id,:]``.

    Notes
    -----
    Geometric statistics (geometric mean, standard deviation, and geometric
    coefficient of variation) are defined only for positive data. If there
    are negative events in any relevant channel of any member of `samples`,
    geometric statistics will only be calculated on the positive events,
    and a warning message will be written to the "Analysis Notes" field.

    """
    # The index name is not preserved if samples_table is empty.
    # Save the index name for later
    samples_table_index_name = samples_table.index.name

    # Add per-row info
    notes = []
    n_events = []
    acq_time = []
    for row_id in samples_table.index:
        # Check if sample is an exception, otherwise assume it's an FCSData
        if isinstance(samples[row_id], ExcelUIException):
            # Print error message
            notes.append("ERROR: {}".format(str(samples[row_id])))
            n_events.append(np.nan)
            acq_time.append(np.nan)
        else:
            notes.append('')
            n_events.append(samples[row_id].shape[0])
            acq_time.append(samples[row_id].acquisition_time)

    samples_table['Analysis Notes'] = notes
    samples_table['Number of Events'] = n_events
    samples_table['Acquisition Time (s)'] = acq_time

    # List of channels that require stats columns
    headers = list(samples_table.columns)
    stats_headers = [h for h in headers if re_units.match(h)]
    stats_channels = [re_units.match(h).group(1) for h in stats_headers]

    # Iterate through channels
    for header, channel in zip(stats_headers, stats_channels):
        # Add empty columns to table
        samples_table[channel + ' Detector Volt.'] = np.nan
        samples_table[channel + ' Amp. Type'] = ""
        samples_table[channel + ' Mean'] = np.nan
        samples_table[channel + ' Geom. Mean'] = np.nan
        samples_table[channel + ' Median'] = np.nan
        samples_table[channel + ' Mode'] = np.nan
        samples_table[channel + ' Std'] = np.nan
        samples_table[channel + ' CV'] = np.nan
        samples_table[channel + ' Geom. Std'] = np.nan
        samples_table[channel + ' Geom. CV'] = np.nan
        samples_table[channel + ' IQR'] = np.nan
        samples_table[channel + ' RCV'] = np.nan
        for row_id in samples_table.index:
            # If error, skip
            if isinstance(samples[row_id], ExcelUIException):
                continue
            # If units are specified, calculate stats. If not, leave empty.
            if pd.notnull(samples_table[header][row_id]):
                # Acquisition settings

                # Detector voltage
                samples_table.at[row_id, channel + ' Detector Volt.'] = \
                    samples[row_id].detector_voltage(channel)

                # Amplification type
                if samples[row_id].amplification_type(channel)[0]:
                    amplification_type = "Log"
                else:
                    amplification_type = "Linear"
                samples_table.at[row_id, channel + ' Amp. Type'] = \
                    amplification_type

                # Statistics from event list
                samples_table.at[row_id, channel + ' Mean'] = \
                    FlowCal.stats.mean(samples[row_id], channel)
                samples_table.at[row_id, channel + ' Median'] = \
                    FlowCal.stats.median(samples[row_id], channel)
                samples_table.at[row_id, channel + ' Mode'] = \
                    FlowCal.stats.mode(samples[row_id], channel)
                samples_table.at[row_id, channel + ' Std'] = \
                    FlowCal.stats.std(samples[row_id], channel)
                samples_table.at[row_id, channel + ' CV'] = \
                    FlowCal.stats.cv(samples[row_id], channel)
                samples_table.at[row_id, channel + ' IQR'] = \
                    FlowCal.stats.iqr(samples[row_id], channel)
                samples_table.at[row_id, channel + ' RCV'] = \
                    FlowCal.stats.rcv(samples[row_id], channel)

                # For geometric statistics, first check for non-positive events.
                # If found, throw a warning and calculate statistics on positive
                # events only.
                if np.any(samples[row_id][:, channel] <= 0):
                    # Separate positive events
                    sample_positive = \
                        samples[row_id][samples[row_id][:, channel] > 0]
                    # Throw warning
                    msg = "Geometric statistics for channel" + \
                        " {} calculated on positive events".format(channel) + \
                        " only ({:.1f}%). ".format(
                            100.*sample_positive.shape[0]/samples[row_id].shape[0])
                    warnings.warn("On sample {}: {}".format(row_id, msg))
                    # Write warning message to table
                    if samples_table.loc[row_id, 'Analysis Notes']:
                        msg = samples_table.loc[row_id, 'Analysis Notes'] + msg
                    samples_table.at[row_id, 'Analysis Notes'] = msg
                else:
                    sample_positive = samples[row_id]

                # Calculate and write geometric statistics
                samples_table.at[row_id, channel + ' Geom. Mean'] = \
                    FlowCal.stats.gmean(sample_positive, channel)
                samples_table.at[row_id, channel + ' Geom. Std'] = \
                    FlowCal.stats.gstd(sample_positive, channel)
                samples_table.at[row_id, channel + ' Geom. CV'] = \
                    FlowCal.stats.gcv(sample_positive, channel)

    # Restore index name if table is empty
    if len(samples_table) == 0:
        samples_table.index.name = samples_table_index_name

def generate_histograms_table(samples_table, samples, max_bins=1024):
    """
    Generate a table of histograms as a DataFrame.

    Parameters
    ----------
    samples_table : DataFrame
        Table specifying samples to analyze. For more information about the
        fields required in this table, please consult the module's
        documentation.
    samples : dict or OrderedDict
        FCSData objects from which to calculate statistics. ``samples[id]``
        should correspond to ``samples_table.loc[id,:]``.
    max_bins : int, optional
        Maximum number of bins to use.

    Returns
    -------
    hist_table : DataFrame
        A multi-indexed DataFrame. Rows contain the histogram bins and
        counts for every sample and channel specified in samples_table.
        `hist_table` is indexed by the sample's ID, the channel name,
        and whether the row corresponds to bins or counts.

    """
    # Extract channels that require stats histograms
    headers = list(samples_table.columns)
    hist_headers = [h for h in headers if re_units.match(h)]
    hist_channels = [re_units.match(h).group(1) for h in hist_headers]

    # The number of columns in the DataFrame has to be set to the maximum
    # number of bins of any of the histograms about to be generated.
    # The following iterates through these histograms and finds the
    # largest.
    n_columns = 0
    for sample_id in samples_table.index:
        if isinstance(samples[sample_id], ExcelUIException):
            continue
        for header, channel in zip(hist_headers, hist_channels):
            if pd.notnull(samples_table[header][sample_id]):
                if n_columns < samples[sample_id].resolution(channel):
                    n_columns = samples[sample_id].resolution(channel)
    # Saturate at max_bins
    if n_columns > max_bins:
        n_columns = max_bins

    # Declare multi-indexed DataFrame
    index = pd.MultiIndex.from_arrays([[],[],[]],
                                      names = ['Sample ID', 'Channel', ''])
    columns = ['Bin {}'.format(i + 1) for i in range(n_columns)]
    hist_table = pd.DataFrame([], index=index, columns=columns)

    # Generate histograms
    for sample_id in samples_table.index:
        if isinstance(samples[sample_id], ExcelUIException):
            continue
        for header, channel in zip(hist_headers, hist_channels):
            if pd.notnull(samples_table[header][sample_id]):
                # Get units in which bins are being reported
                unit = samples_table[header][sample_id]
                # Decide which scale to use
                # Channel units result in linear scale. Otherwise, use logicle.
                if unit == 'Channel':
                    scale = 'linear'
                else:
                    scale = 'logicle'
                # Define number of bins
                nbins = min(samples[sample_id].resolution(channel), max_bins)
                # Calculate bin edges and centers
                # We generate twice the necessary number of bins. We then take
                # every other value as the proper bin edges, and the remaining
                # values as the bin centers.
                bins_extended = samples[sample_id].hist_bins(channel,
                                                             2*nbins,
                                                             scale)
                bin_edges = bins_extended[::2]
                bin_centers = bins_extended[1::2]
                # Store bin centers
                hist_table.loc[(sample_id,
                                channel,
                                'Bin Centers ({})'.format(unit)),
                                columns[0:len(bin_centers)]] = bin_centers
                # Calculate and store histogram counts
                hist, __ = np.histogram(samples[sample_id][:,channel],
                                        bins=bin_edges)
                hist_table.loc[(sample_id, channel, 'Counts'),
                               columns[0:len(bin_centers)]] = hist

    return hist_table

def generate_about_table(extra_info={}):
    """
    Make a table with information about FlowCal and the current analysis.

    Parameters
    ----------
    extra_info : dict, optional
        Additional keyword:value pairs to include in the table.

    Returns
    -------
    about_table : DataFrame
        Table with information about FlowCal and the current analysis, as
        keyword:value pairs. The following keywords are included: FlowCal
        version, and date and time of analysis. Keywords and values from
        `extra_info` are also included.

    """
    # Make keyword and value arrays
    keywords = []
    values = []
    # FlowCal version
    keywords.append('FlowCal version')
    values.append(FlowCal.__version__)
    # Analysis date and time
    keywords.append('Date of analysis')
    values.append(time.strftime("%Y/%m/%d"))
    keywords.append('Time of analysis')
    values.append(time.strftime("%I:%M:%S%p"))
    # Add additional keyword:value pairs
    for k, v in six.iteritems(extra_info):
        keywords.append(k)
        values.append(v)

    # Make table as data frame
    about_table = pd.DataFrame(values, index=keywords)

    # Set column names
    about_table.columns = ['Value']
    about_table.index.name = 'Keyword'

    return about_table

def show_open_file_dialog(filetypes):
    """
    Show an open file dialog and return the path of the file selected.

    Parameters
    ----------
    filetypes : list of tuples
        Types of file to show on the dialog. Each tuple on the list must
        have two elements associated with a filetype: the first element is
        a description, and the second is the associated extension.

    Returns
    -------
    filename : str
        The path of the filename selected, or an empty string if no file
        was chosen.

    """

    # initialize tkinter root window
    root = Tk()
    # remove main root window (will cause kernel panic on OSX if not present)
    root.withdraw()
    # link askopenfilename window to root window
    filename = askopenfilename(parent = root, filetypes=filetypes)
    # refresh root window to remove askopenfilename window
    root.update()
    
    return filename

def run(input_path=None,
        output_path=None,
        verbose=True,
        plot=True,
        hist_sheet=False):
    """
    Run the MS Excel User Interface.

    This function performs the following:

     1. If `input_path` is not specified, show a dialog to choose an input
        Excel file.
     2. Read the Instruments table from the Instruments sheet.
     3. If a Beads sheet is specified, read the Beads table, process beads
        samples, and calculate statistics.
     4. If a Samples sheet is specified, read the Samples table, process
        samples, calculate statistics, and (if requested) generate a
        histogram table describing each fluorescence channel of each sample.
     5. Generate a table describing the run (e.g. with run time, date,
        FlowCal version, etc.).
     6. Save tables and calculated statistics to an output Excel file.

    Parameters
    ----------
    input_path : str
        Path to the Excel file to use as input. If None, show a dialog to
        select an input file.
    output_path : str
        Path to which to save the output Excel file. If None, use
        "<input_path>_output".
    verbose : bool, optional
        Whether to print information messages during the execution of this
        function.
    plot : bool, optional
        Whether to generate and save density/histogram plots of each
        sample, and each beads sample.
    hist_sheet : bool, optional
        Whether to generate a sheet in the output Excel file specifying
        histogram bin information. Ignored if Samples sheet is not specified.

    """
    # If input file has not been specified, show open file dialog
    if input_path is None:
        input_path = show_open_file_dialog(filetypes=[('Excel files',
                                                       '*.xlsx')])
        if not input_path:
            if verbose:
                print("No input file selected.")
            return

    # Extract directory, filename, and filename with no extension from path
    input_dir, input_filename = os.path.split(input_path)
    input_filename_no_ext, __ = os.path.splitext(input_filename)

    # Process tables
    if verbose:
        print("Reading {}...".format(input_filename))
    table_list = []  # List of (str, DataFrame) tuples, one for each
                     # sheet to be written to the output Excel file.

    instruments_table = read_table(input_path,
                                   sheetname='Instruments',
                                   index_col='ID')
    table_list.append(('Instruments', instruments_table))

    # Process beads samples
    try:
        beads_table = read_table(input_path,
                                 sheetname='Beads',
                                 index_col='ID')
    except KeyError as e:
        if 'Beads' in str(e):
            # no Beads tab (openpyxl)
            beads_table = None
        else:
            raise
    except Exception as e:
        if 'xlrd' in sys.modules and isinstance(e, xlrd.biffh.XLRDError) \
                and 'Beads' in str(e):
            # no Beads tab (xlrd)
            beads_table = None
        else:
            raise

    if beads_table is not None:
        beads_samples, mef_transform_fxns, mef_outputs = process_beads_table(
            beads_table,
            instruments_table,
            base_dir=input_dir,
            verbose=verbose,
            plot=plot,
            plot_dir='plot_beads',
            full_output=True)

        # Add stats to beads table
        if verbose:
            print("")
            print("Calculating statistics for beads...")
        add_beads_stats(beads_table, beads_samples, mef_outputs)
        table_list.append(('Beads', beads_table))
    else:
        beads_samples = mef_transform_fxns = mef_outputs = None

    # Process samples
    try:
        samples_table = read_table(input_path,
                                   sheetname='Samples',
                                   index_col='ID')
    except KeyError as e:
        if 'Samples' in str(e):
            # no Samples tab (openpyxl)
            samples_table = None
        else:
            raise
    except Exception as e:
        if 'xlrd' in sys.modules and isinstance(e, xlrd.biffh.XLRDError) \
                and 'Samples' in str(e):
            # no Samples tab (xlrd)
            samples_table = None
        else:
            raise

    if samples_table is not None:
        samples = process_samples_table(
            samples_table,
            instruments_table,
            mef_transform_fxns=mef_transform_fxns,
            beads_table=beads_table,
            base_dir=input_dir,
            verbose=verbose,
            plot=plot,
            plot_dir='plot_samples')

        # Add stats to samples table
        if verbose:
            print("")
            print("Calculating statistics for all samples...")
        add_samples_stats(samples_table, samples)
        table_list.append(('Samples', samples_table))

        # Generate histograms
        if hist_sheet:
            if verbose:
                print("Generating histograms table...")
            histograms_table = generate_histograms_table(samples_table, samples)
            table_list.append(('Histograms', histograms_table))
    else:
        samples = None

    # Generate about table
    about_table = generate_about_table({'Input file path': input_path})
    table_list.append(('About Analysis', about_table))

    # Write output excel file
    if verbose:
        print("Saving output Excel file...")
    if output_path is None:
        output_filename = "{}_output.xlsx".format(input_filename_no_ext)
        output_path = os.path.join(input_dir, output_filename)
    write_workbook(output_path, table_list)

    if verbose:
        print("\nDone.")


def run_command_line(args=None):
    """
    Entry point for the FlowCal and flowcal console scripts.

    Parameters
    ----------
    args: list of strings, optional
        Command line arguments. If None or not specified, get arguments
        from ``sys.argv``.

    See Also
    --------
    FlowCal.excel_ui.run

    References
    ----------
    http://amir.rachum.com/blog/2017/07/28/python-entry-points/

    """
    # Get arguments from ``sys.argv`` if necessary.
    # ``sys.argv`` has the name of the script as its first element. We remove
    # this element because it will break ``parser.parse_args()`` later. In fact,
    # ``parser.parse_args()``, if provided with no arguments, will also use
    # ``sys.argv`` after removing the first element.
    if args is None:
        args = sys.argv[1:]

    import argparse
    # Read command line arguments
    parser = argparse.ArgumentParser(
        description="process flow cytometry files with FlowCal's Excel UI.")
    parser.add_argument(
        "-i",
        "--inputpath",
        type=str,
        nargs='?',
        help="input Excel file name. If not specified, show open file window")
    parser.add_argument(
        "-o",
        "--outputpath",
        type=str,
        nargs='?',
        help="output Excel file name. If not specified, use [INPUTPATH]_output")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="print information about individual processing steps")
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="generate and save density plots/histograms of beads and samples")
    parser.add_argument(
        "-H",
        "--histogram-sheet",
        action="store_true",
        help="generate sheet in output Excel file specifying histogram bins")
    args = parser.parse_args(args=args)

    # Run Excel UI
    run(input_path=args.inputpath,
        output_path=args.outputpath,
        verbose=args.verbose,
        plot=args.plot,
        hist_sheet=args.histogram_sheet)

if __name__ == '__main__':
    run_command_line()
