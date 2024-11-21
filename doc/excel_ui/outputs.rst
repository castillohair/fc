Outputs of the Excel UI
=======================

The ``FlowCal`` Excel UI creates analysis plots and an output Excel file. If a Beads sheet is specified in the input Excel file, a `plot_beads` folder is created containing relevant plots. Similarly, if a Samples sheet is specified, a `plot_samples` folder is created with relevant plots. Beads and Samples sheets are also populated in the output Excel file if specified in the input Excel file.

.. _excel-ui-outputs-plots:

Plots
-----

Note: `<ID>` refers to the unique ID of a sample as labeled in the ID column of the input Excel file.

1. The folder ``plot_beads`` contains plots of the individual steps of processing of the calibration particle samples:

    a. ``density_hist_<ID>.png``: A forward/side scatter 2D density diagram of the calibration particle sample, and a histogram for each relevant fluorescence channel.

    .. image:: /_static/img/excel_ui/output_beads_density.png

    b. ``clustering_<ID>.png``: A plot of the sub-populations identified during the clustering step, where the different sub-populations are shown in different colors. Depending on the number of channels used for clustering, this plot is a histogram (when using only one channel), a 2D scatter plot (when using two channels), or a 3D scatter plot with three 2D projections (when using three channels or more).

    .. image:: /_static/img/excel_ui/output_beads_clustering.png

    .. note:: It is normally easy to distinguish the different bead populations in this plot, and the different colors should correspond to this expectation. If the populations have been identified incorrectly, changing the number of channels used for clustering or the density gate fraction can improve the results. These two parameters can be changed in the **Beads** sheet of the input Excel file.

    c. ``populations_<channel>_<ID>.png``: A histogram showing the identified microbead sub-populations in different colors, for each fluorescence channel in which a MEF standard curve is to be calculated. In addition, a vertical line is shown representing the median of each population, which is later used to calculate the standard curve. Sub-populations that were not used to generate the standard curve are shown in gray.

    .. image:: /_static/img/excel_ui/output_beads_populations.png

    .. note:: All populations should be unimodal. Bimodal populations indicate incorrect clustering. This can be fixed by changing the number of channels used for clustering or the density gate fraction in the **Beads** sheet of the input Excel file.

    d. ``std_crv_<channel>_<ID>.png``: A plot of the fitted standard curve, for each channel in which MEF values were specified.

    .. image:: /_static/img/excel_ui/output_beads_sc.png

    .. note:: All the blue dots should line almost perfectly on the green line, otherwise the estimation of the standard curve might not be good. If this is not the case, you should make sure that clustering is being performed correctly by looking at the previous plots. If one dot differs significantly from the curve despite perfect clustering, you might want to manually remove it. This can be done by replacing its MEF value with the word "None" in the **Beads** sheet of the input Excel file.

2. The folder ``plot_samples`` contains plots of the experimental cell samples. Each experimental sample of name “ID” as specified in the Excel input sheet results in a file named ``<ID>.png``. This image contains a forward/side scatter 2D density diagram with the gated region indicated, and a histogram for each user-specified fluorescence channel.

.. image:: /_static/img/excel_ui/output_sample.png

.. _excel-ui-outputs-excel:

Output Excel File
-----------------

The file ``<Name of the input Excel file>_output.xlsx`` contains calculated statistics for beads and samples. To produce this file, ``FlowCal`` copies the **Instruments**, **Beads**, and **Samples** sheets from the :doc:`input Excel<input_format>` file, unmodified, to the output file, and adds columns to the **Beads** and **Samples** sheet with statistics.

In both sheets, the number of events after gating and the acquisition time are reported for each sample. In addition, a column named **Analysis Notes** indicates the user about any errors that occurred during processing.

Statistics per beads file, per fluorescence channel include: the channel gain, the amplifier type, the equation of the beads fluorescence model used, and the values of the fitted parameters.

.. image:: /_static/img/excel_ui/output_spreadsheet_beads.png

Statistics per cell sample, per fluorescence channel include: channel gain, mean, geometric mean, median, mode, arithmetic and geometric standard deviation, arithmetic and geometric coefficient of variation (CV), interquartile range (IQR), and robust coefficient of variation (RCV). Note that if an error has been found, the **Analysis Notes** field will be populated, and statistics and plots will not be reported.

.. image:: /_static/img/excel_ui/output_spreadsheet_samples.png

In addition, a **Histograms** tab is generated, with bin/counts pairs for each sample and relevant fluorescence channel in the specified units.

.. image:: /_static/img/excel_ui/output_spreadsheet_histograms.png

One last tab named **About Analysis** is added with information about the corresponding input Excel file, the date and time of the run, and the FlowCal version used.

.. image:: /_static/img/excel_ui/output_spreadsheet_about.png
