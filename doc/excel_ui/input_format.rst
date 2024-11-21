Format of the Input Excel File
==============================

``FlowCal``'s Excel interface requires a properly formatted Excel file that depicts the samples to be analyzed and the data processing parameters. The Excel input file must have an **Instruments** sheet and typically also has **Beads** and **Samples** sheets. Other sheets may also be present, but ``FlowCal`` will ignore them.

.. warning:: Sheet and column names are case-sensitive.

An example of a properly formatted Excel input file is provided in the ``examples`` folder of ``FlowCal``. The following sections describe the format of the input Excel file, while using the example file as a guide.

Instruments sheet
-----------------

This sheet must be filled with basic information about the flow cytometer used to acquire the samples. Each row represents an instrument. Typically, the user would only need to specify one instrument. However, ``FlowCal`` allows the simultaneous processing of samples taken with different instruments. The figure below shows an example of an **Instruments** sheet.

.. image:: /_static/img/excel_ui/input_instruments.png

For each row, the following columns must be filled.

1. **ID** (column A in the figure above) used to reference the instrument from the other sheets. Each row must have a unique ID.
2. **Forward Scatter Channel** (C) and **Side Scatter Channel** (D): the names of these channels exactly as they appear in the acquisition software.
3. **Fluorescence channels** (E): The names of the relevant fluorescence channels as a comma-separated list, exactly as they appear in the acquisition software.
4. **Time Channel** (F): The name of the channel registering the time of each event. The FCS standard dictates that this should be called “Time”, but some non-standard files may use a different name. This can be found in the acquisition software.

Additional columns, like **Description** (B in the figure above), can be added in any place for the user’s records, and will be copied unmodified to the output Excel file by ``FlowCal``.

Beads sheet
-----------

This sheet contains details about calibration microbeads and how to process them. Each row represents a different sample of beads. The figure below shows an example of an **Beads** sheet.

.. image:: /_static/img/excel_ui/input_beads.png

For each row, the following columns must be filled:

1. **ID** (column A in the figure above): used to reference the beads sample from the Samples sheet, and to name the figures produced by ``FlowCal``. Each row must have a unique ID.
2. **Instrument ID** (B): The ID of the instrument used to take the sample.
3. **File Path** (C): the name of the corresponding FCS file.
4. **<Channel name> MEF Values** (E): MEF values provided by the manufacturer, for each channel in which a standard curve must be calculated. If MEF values are provided for a channel, the corresponding instrument should include this channel name in the **Fluorescence Channels** field. More **<Channel name> MEF Values** columns can be added if needed, or removed if not used.
5. **Gate Fraction** (F): a gate fraction parameter used for :doc:`density gating</fundamentals/density_gate>`.
6. **Clustering Channels** (G): the fluorescence channels used for clustering, as a comma separated list.

Additional columns, like **Beads Lot** (column D), can be added in any place for the user’s records, and will be copied unmodified to the output Excel file by ``FlowCal``.

Samples sheet
-------------

In this sheet, the user specifies cell samples and tells ``FlowCal`` how to process them. Each row contains the information used in the analysis of one FCS file. One file can be analyzed several times with different options (e.g. gating fractions or fluorescence units) by adding more rows that reference the same file. The figure below shows an example of a **Samples** sheet.

.. image:: /_static/img/excel_ui/input_samples.png

For each row, the following columns must be filled:

1. **ID** (column A in the figure above): used to reference the sample while generating figures, and in the output Excel file. Each row must have a unique ID.
2. **Instrument ID** (B): The ID of the instrument used to take the sample.
3. **Beads ID** (C): The ID of the beads sample that will be used to perform the MEF transformation. Can be left blank if MEF units are not desired.
4. **File Path** (D): the name of the corresponding FCS file.
5. **<Channel name> Units** (E): The units in which to report statistics and make plots, for each fluorescence channel. If left blank, no statistics or plots will be made for that channel. More of these columns can be added or removed if necessary. If this field is specified for a channel, the corresponding instrument should include this channel in its **Fluorescence Channels** field. The available options are:

    a. **Channel**: Raw “Channel Number” units, exactly as they are stored in the FCS file.
    b. **RFI** or **a.u.**: Relative Fluorescence Intensity units, also known as Arbitrary Units. 
    c. **MEF**: MEF units.
6. **Gate Fraction** (F): Fraction of samples to keep when performing :doc:`density gating</fundamentals/density_gate>`.

Additional columns, such as **Strain**, **Plasmid**, and **DAPG (uM)** (columns G, H, and I), can be added in any place for the user’s records, and will be copied unmodified to the output Excel file by ``FlowCal``.

.. warning:: If MEF units are requested for a fluorescence channel of a sample, an FCS file with calibration beads data should be specified in the **Beads ID** column. Both beads and samples should have been acquired at the same settings for the specified fluorescence channel, otherwise ``FlowCal`` will throw an error.
