Analysis Performed by the Excel UI
==================================

The analysis that FlowCal's Excel UI performs is divided roughly in two phases: processing of calibration beads and processing of samples. We will now describe the steps involved in each.

Processing of Calibration Beads
-------------------------------
If a **Beads** sheet is specified, the following steps are performed for each calibration beads sample:

1. :doc:`Density gating</fundamentals/density_gate>` is applied in the forward/side scatter channels. This is an automated procedure that eliminates microbead aggregates and debris.
2. The individual microbead subpopulations are identified using automated clustering.
3. For each subpopulation, the median fluorescence is calculated. 
4. Microbead subpopulations are discarded if they are found to be close to the saturation limits of the detector. Only populations that are not saturating are retained.
5. Using the fluorescence values of the retained populations in channel units and the corresponding MEF values provided by the user, a standard curve is generated. This standard curve is used to transform cell fluorescence from raw units to MEF.

:ref:`Plots<excel-ui-outputs-plots>` are generated for each one of these steps, and some intermediate results are saved to the :ref:`output Excel file<excel-ui-outputs-excel>`.

For an introductory discussion of flow cytometry calibration, go to the :doc:`fundamentals of calibration</fundamentals/calibration>` section.

Processing of Cell Samples
--------------------------
If a **Samples** sheet is specified, the following steps are performed for each sample:

1. :doc:`Density gating</fundamentals/density_gate>` is applied in the forward/side scatter channels.
2. Fluorescence data for each specified fluorescence channel is transformed to the units specified in the **Units** column of the :doc:`input Excel file<input_format>`.
3. :ref:`Statistics<excel-ui-outputs-excel>` of the specified fluorescence channels are calculated, including mean, standard deviation, and others. A histogram of each fluorescence channel is also generated.

Statistics and histograms are saved to the :ref:`output Excel file<excel-ui-outputs-excel>`.
