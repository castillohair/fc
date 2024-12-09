"""
Functions for transforming flow cytometry data

All transformations are of the following form::

    data_t = transform(data, channels, *args, **kwargs):

where `data` and `data_t` are NxD FCSData objects or numpy arrays,
representing N events with D channels, `channels` indicate the channels in
which to apply the transformation, and `args` and `kwargs` are
transformation-specific parameters. Each transformation function can apply
its own restrictions or defaults on `channels`.

If `data` is an FCSData object, `transform` should rescale ``data.range``
if necessary.

"""

import six
import numpy as np

def transform(data, channels, transform_fxn, def_channels = None):
    """
    Apply some transformation function to flow cytometry data.

    This function is a template transformation function, intended to be
    used by other specific transformation functions. It performs basic
    checks on `channels` and `data`. It then applies `transform_fxn` to the
    specified channels. Finally, it rescales  ``data.range`` and if
    necessary.

    Parameters
    ----------
    data : FCSData or numpy array
        NxD flow cytometry data where N is the number of events and D is
        the number of parameters (aka channels).
    channels : int, str, list of int, list of str, optional
        Channels on which to perform the transformation. If `channels` is
        None, use def_channels.
    transform_fxn : function
        Function that performs the actual transformation.
    def_channels : int, str, list of int, list of str, optional
        Default set of channels in which to perform the transformation.
        If `def_channels` is None, use all channels.

    Returns
    -------
    data_t : FCSData or numpy array
        NxD transformed flow cytometry data.

    """
    # Copy data array
    data_t = data.copy().astype(np.float64)

    # Default
    if channels is None:
        if def_channels is None:
            channels = range(data_t.shape[1])
        else:
            channels = def_channels

    # Convert channels to iterable
    if not (hasattr(channels, '__iter__') \
            and not isinstance(channels, six.string_types)):
        channels = [channels]

    # Apply transformation
    data_t[:,channels] = transform_fxn(data_t[:,channels])

    # Apply transformation to ``data.range``
    if hasattr(data_t, '_range'):
        for channel in channels:
            # Transform channel name to index if necessary
            channel_idx = data_t._name_to_index(channel)
            if data_t._range[channel_idx] is not None:
                data_t._range[channel_idx] = \
                    transform_fxn(data_t._range[channel_idx])

    return data_t

def to_rfi(data,
           channels=None,
           amplification_type=None,
           amplifier_gain=None,
           resolution=None):
    """
    Transform flow cytometry data to Relative Fluorescence Units (RFI).

    If ``amplification_type[0]`` is different from zero, data has been
    taken using a log amplifier. Therefore, to transform to RFI, the
    following operation is applied::

        y = a[1]*10^(a[0] * (x/r))

    Where ``x`` and ``y`` are the original and transformed data,
    respectively; ``a`` is `amplification_type` argument, and ``r`` is
    `resolution`. This will transform flow cytometry data taken with a log
    amplifier and an ADC of range ``r`` to linear RFIs, such
    that it covers ``a[0]`` decades of signal with a minimum value of
    ``a[1]``.

    If ``amplification_type[0]==0``, however, a linear amplifier has been
    used and the following operation is applied instead::

        y = x/g

    Where ``g`` is `amplifier_gain`. This will transform flow cytometry
    data taken with a linear amplifier of gain ``g`` back to RFIs.

    Parameters
    ----------
    data : FCSData or numpy array
        NxD flow cytometry data where N is the number of events and D is
        the number of parameters (aka channels).
    channels : int, str, list of int, list of str, optional
        Channels on which to perform the transformation. If `channels` is
        None, perform transformation in all channels.
    amplification_type : tuple or list of tuple
        The amplification type of the specified channel(s). This should be
        reported as a tuple, in which the first element indicates how many
        decades the logarithmic amplifier covers, and the second indicates
        the linear value that corresponds to a channel value of zero. If
        the first element is zero, the amplification type is linear. This
        is similar to the $PnE keyword from the FCS standard. If None, take
        `amplification_type` from ``data.amplification_type(channel)``.
    amplifier_gain : float or list of floats, optional
        The linear amplifier gain of the specified channel(s). Only used if
        ``amplification_type[0]==0`` (linear amplifier). If None,
        take `amplifier_gain` from ``data.amplifier_gain(channel)``. If
        `data` does not contain ``amplifier_gain()``, use 1.0.
    resolution : int, float, or list of int or float, optional
        Maximum range, for each specified channel. Only needed if
        ``amplification_type[0]!=0`` (log amplifier). If None, take
        `resolution` from ``len(data.domain(channel))``.

    Returns
    -------
    FCSData or numpy array
        NxD transformed flow cytometry data.

    """
    # Default: all channels
    if channels is None:
        channels = range(data.shape[1])

    if not (hasattr(channels, '__iter__') \
            and not isinstance(channels, six.string_types)):
        # If channels is not an iterable, convert it, along with resolution,
        # amplification_type, and amplifier_gain.
        channels = [channels]
        amplification_type = [amplification_type]
        amplifier_gain = [amplifier_gain]
        resolution = [resolution]
    else:
        # If channels is an iterable, check that the other attributes are either
        # None, or iterables of the same length.

        if amplification_type is None:
            # If None, propagate None for all channels
            amplification_type = [None]*len(channels)
        elif hasattr(amplification_type, '__iter__'):
            # If it's a list, it should be the same length as channels
            if len(amplification_type) != len(channels):
                raise ValueError("channels and amplification_type should have "
                    "the same length")
        else:
            # If it's not a list or None, raise error
            raise ValueError("channels and amplification_type should have the "
                    "same length")

        if amplifier_gain is None:
            # If None, propagate None for all channels
            amplifier_gain = [None]*len(channels)
        elif hasattr(amplifier_gain, '__iter__'):
            # If it's a list, it should be the same length as channels
            if len(amplifier_gain) != len(channels):
                raise ValueError("channels and amplifier_gain should have "
                    "the same length")
        else:
            # If it's not a list or None, raise error
            raise ValueError("channels and amplifier_gain should have the "
                    "same length")

        if resolution is None:
            # If None, propagate None for all channels
            resolution = [None]*len(channels)
        elif hasattr(resolution, '__iter__'):
            # If it's a list, it should be the same length as channels
            if len(resolution) != len(channels):
                raise ValueError("channels and resolution should have "
                    "the same length")
        else:
            # If it's not a list or None, raise error
            raise ValueError("channels and resolution should have the "
                    "same length")

    # Convert channels to integers
    if hasattr(data, '_name_to_index'):
        channels = data._name_to_index(channels)
    else:
        channels = channels

    # Copy data array
    data_t = data.copy().astype(np.float64)

    # Iterate over channels
    for channel, r, at, ag in \
            zip(channels, resolution, amplification_type, amplifier_gain):
        # If amplification type is None, try to obtain from data
        if at is None:
            if hasattr(data, 'amplification_type'):
                at = data.amplification_type(channel)
            else:
                raise ValueError('amplification_type should be specified')
        # Define transformation, depending on at[0]
        if at[0]==0:
            # Linear amplifier
            # If no amplifier gain has been specified, try to obtain from data,
            # otherwise assume one
            if ag is None:
                if hasattr(data, 'amplifier_gain') and \
                        hasattr(data.amplifier_gain, '__call__'):
                    ag = data.amplifier_gain(channel)
                    # If the linear gain has not been specified, it should be
                    # assumed to be one.
                    if ag is None or ag == 0:
                        ag = 1.
                else:
                    ag = 1.
            tf = lambda x: x/ag
        else:
            # Log amplifier
            # If no range has been specified, try to obtain from data.
            if r is None:
                if hasattr(data, 'resolution'):
                    r = data.resolution(channel)
                else:
                    raise ValueError('range should be specified')
            tf = lambda x: at[1] * 10**(at[0]/float(r) * x)
        # Apply transformation to event list
        data_t[:,channel] = tf(data_t[:,channel])
        # Apply transformation to range
        if hasattr(data_t, '_range') and data_t._range[channel] is not None:
            data_t._range[channel] = [tf(data_t._range[channel][0]),
                                      tf(data_t._range[channel][1])]

    return data_t


def to_mef(data, channels, sc_list, sc_channels = None):
    """
    Transform flow cytometry data using a standard curve function.

    This function accepts a list of standard curves (`sc_list`) and a list
    of channels to which those standard curves should be applied
    (`sc_channels`). `to_mef` automatically checks whether a standard curve
    is available for each channel specified in `channels`, and throws an
    error otherwise.

    This function is intended to be reduced to the following signature::

        to_mef_reduced(data, channels)

    by using ``functools.partial`` once a list of standard curves and their
    respective channels is available.

    Parameters
    ----------
    data : FCSData or numpy array
        NxD flow cytometry data where N is the number of events and D is
        the number of parameters (aka channels).
    channels : int, str, list of int, list of str
        Channels on which to perform the transformation. If `channels` is
        None, perform transformation in all channels specified on
        `sc_channels`.
    sc_list : list of functions
        Functions implementing the standard curves for each channel in
        `sc_channels`.
    sc_channels : list of int or list of str, optional
        List of channels corresponding to each function in `sc_list`. If
        None, use all channels in `data`.

    Returns
    -------
    FCSData or numpy array
        NxD transformed flow cytometry data.

    Raises
    ------
    ValueError
        If any channel specified in `channels` is not in `sc_channels`.

    """
    # Default sc_channels
    if sc_channels is None:
        if data.ndim == 1:
            sc_channels = range(data.shape[0])
        else:
            sc_channels = range(data.shape[1])
    # Check that sc_channels and sc_list have the same length
    if len(sc_channels) != len(sc_list):
        raise ValueError("sc_channels and sc_list should have the same length")
    # Convert sc_channels to indices
    if hasattr(data, '_name_to_index'):
        sc_channels = data._name_to_index(sc_channels)

    # Default channels
    if channels is None:
        channels = sc_channels
    # Convert channels to iterable
    if not (hasattr(channels, '__iter__') \
            and not isinstance(channels, six.string_types)):
        channels = [channels]
    # Convert channels to index
    if hasattr(data, '_name_to_index'):
        channels_ind = data._name_to_index(channels)
    else:
        channels_ind = channels
    # Check if every channel is in sc_channels
    for chi, chs in zip(channels_ind, channels):
        if chi not in sc_channels:
            raise ValueError("no standard curve for channel {}".format(chs))

    # Copy data array
    data_t = data.copy().astype(np.float64)

    # Iterate over channels
    for chi, sc in zip(sc_channels, sc_list):
        if chi not in channels_ind:
            continue
        # Apply transformation
        data_t[:,chi] = sc(data_t[:,chi])
        # Apply transformation to range
        if hasattr(data_t, '_range') and data_t._range[chi] is not None:
            data_t._range[chi] = [sc(data_t._range[chi][0]),
                                  sc(data_t._range[chi][1])]

    return data_t
