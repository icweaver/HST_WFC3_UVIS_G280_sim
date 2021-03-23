# UVIS_precision_calculator.py
"""
This program was created by H.R.Wakeford 2020 Feb 6

It has been checked against the HAT-P-41b data obtained through GO-15288 (PIs: D.K.Sing &
N.K.Lewis) when compared to the ETC. These measurements reached photon noise limits which
is what the ETC estimates.

Credit: H.R. Wakeford, when using ETC inputs please also credit the ETC in your proposals.

Citation: Wakeford et al. (2020), AJ

This simulator will produce precisions for a transmission spectrum in custom bins from
200-800nm and the broadband transit depth precision and the specified broadband depth
precision.
"""
import numpy as np
import pandas as pd


def find_nearest(array, value):
    """
    Args:
        array (array-like): Collection of values.
        value (int, float): New, sample value.

    Returns:
        idx (int): The index of `array` nearest to `value`

    Examples:
        >>> find_nearest([6, 7, 8], 6.6)
        1
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def transit_depth_precision(wavelength, spectrum, wav_range, exptransit):
    """Computes the error ``\sigma_w`` in the integrated flux (in ppm):
    .. math::
        \sigma_w = \\frac{1}{S_r t} \quad,

    where ``S_r`` is the count rate and ``t`` is the exposure time.

    Args:
        wavelength (array-like): Collection of wavelengths in the spectrum.
        spectrum (array-like): Corresponding collection of fluxes.
        wav_range (array-like): 2-element collection of lower and upper wavelength range.
        exptransit (int, float): Exposure time, nominally in seconds.

    Returns:
        depth_err_ppm (array-like):
    """
    wav_range = np.where((wavelength > wav_range[0]) & (wavelength < wav_range[1]))[0]

    total_sig = np.sum(spectrum[wav_range])
    precision = 1 / np.sqrt(total_sig)
    depth_err_ppm = (precision / np.sqrt(exptransit)) * 1e6

    return depth_err_ppm


def header(i, wl, binlen, file_name, exposure_time, orbits_in_transit, no_of_transits):
    print(i, wl, binlen, file_name, exposure_time, orbits_in_transit, no_of_transits)
    return (
        f"UVIS G280 Simulation, broadband = {wl[0]/10}-{wl[1]/10}nm, spectroscopic in {binlen/10}nm bins"
        "\nCredit: H.R. Wakeford, Citation: Wakeford et al. (2020, AJ)"
        f"\n{file_name[i]} with an exposure time of {exposure_time[i]} and "
        f"{orbits_in_transit[i]} orbits in transit/eclipse for {no_of_transits[i]} transits/eclipses"
        "\nWavelength (nm), bin size (nm), transit/eclipse precision (ppm), program precision (ppm)"
    )


def UVIS_simulation(
    data_folder,
    file_name,
    exposure_time,
    orbits_in_transit,
    no_of_transits,
    wl=[2000, 8000],
    startw=2000,
    endw=8000,
    binlen=100,
    out_name="transmission",
):
    """
    INPUT:
    data_folder = '..Data/'
        The place where the ETC output file located
        The ETC file needs to contain the following information in the correct structure.
        This file can be self generated if also in this structure.
        (wavelength (A), readnoise, dark_counts, sky_counts, target_counts, total_counts)
        Columns 0 and 4 are used in this calculation, please maintain these column
        positions for self input files
        Wavelength is input from this file in Angstroms

    file_name = ['H41_ETC']
        This can be an array of filenames
        make sure it is always in square brackets

    exposure_time = [190]
        This is an array of exposuretimes corresponding to the
        filenames listed in file_name
        make sure it is always in square brackets

    orbits_in_transit = [2]
        This is an array of the number of HST orbits inside a single
        transit (not the total observation) corresponding to the
        filenames listed in file_name
        make sure it is always in square brackets

    no_of_transits = [2]
        This is an array of the number of transits of each targets
        in your program corresponding to the filenames listed in
        file_name make sure it is always in square brackets

    wl = [2000,8000]
        default value spans the maximum range of the grism. This will
        calculate the broadband signal that will be the first line of
        the produced file.

    startw = 2000
        This is the starting range for your spectral bins.
        Default is at minimum = 2000 angstroms.
        Where this might be changed is for small stars with little
        flux in the short wavelengths. Or if simulating the -1 order

    endw = 8000
        This is the end of the spectral binning range
        Default set to the maximum, 8000 angstroms

    binlen = 100
        This is in angstroms and specifies the size of the
        spectral bins in wavlength space.
        Minimum reccomended = 60 angstorms (Wakeford et al 2020 = 100)

    out_name = 'string'
        This is what you want the output specified as in addition
        to the input name.

    OUTPUT:
    data file for each input planet file:
    output name = {data_folder}{file_name}.UVIS_{out_name}_sim.txt'
        contains the broadband {wl} and spectroscopic precisions ({binlen} bins)
        wavelength (nm), bin size (nm), precision per transit (ppm), precision per target (ppm)

    """
    # UVIS spectral bins for simulation
    # Note Wakeford et al. (2020, AJ) used bins of 100 angstrom by
    # combining two transits of the hot Jupiter HAT-P-41b
    nbins = int((endw - startw) / binlen)
    wav_topbot = np.linspace(startw, endw, nbins)
    wav_bot = wav_topbot[:-1]
    wav_top = wav_topbot[1:]

    wav_array = wav_bot + ((wav_top - wav_bot) / 2)
    wav_err = (wav_top - wav_bot) / 2

    # Set up empty arrays to save the calculations into
    UVIS_spec_precision_per_transit = np.empty(len(wav_array), dtype=float)
    UVIS_spec_precision_total = np.empty(len(wav_array), dtype=float)

    # Convert the input arrays in np.arrays
    exposure_time = np.array(exposure_time)
    orbits_in_transit = np.array(orbits_in_transit)
    no_of_transits = np.array(no_of_transits)

    # Conservative estimate of time between exposures
    overhead = 60  # (fixed)

    # Determine approximately how many exposures to expect per transit
    HST_orbit_len = 48  # minutes (range often 45-52 mins)
    exp_in_transit = np.round(
        (orbits_in_transit * HST_orbit_len * 60) / (exposure_time + overhead)
    )

    # Calculate exposures per HST orbit
    exp_per_orbit = np.round((HST_orbit_len * 60) / (exposure_time + overhead))

    # Compute the UVIS G280 precisions for each file
    for i in range(0, len(file_name)):
        print(file_name[i])

        data_name = "{}{}.csv".format(data_folder, file_name[i])
        wave, signal = pd.read_csv(
            data_name, comment="#", usecols=["wavelength", "target_counts"]
        ).values.T

        wav_range = wl  # [2000,8000]
        exptransit = exp_in_transit[i]

        result = transit_depth_precision(wave, signal, wav_range, exptransit)
        UVIS_broadband_per_transit = result
        UVIS_broadband_total = result / np.sqrt(no_of_transits[i])

        print("exp time = ", exposure_time[i])
        print("exp in transit = ", exp_in_transit[i])
        print("exp per orbit = ", exp_per_orbit[i])
        print("precision/transit = ", UVIS_broadband_per_transit)

        # Compute the precisions for the UVIS G280 spectrum
        for j in range(0, len(wav_array)):
            wav_range = [wav_bot[j], wav_top[j]]
            exptransit = exp_in_transit[i]

            result = transit_depth_precision(wave, signal, wav_range, exptransit)
            UVIS_spec_precision_per_transit[j] = result
            UVIS_spec_precision_total[j] = result / np.sqrt(no_of_transits[i])

        # Create arrays of the results to save to txt file
        wl_cent = wl[0] + ((wl[1] - wl[0]) / 2)
        wl_err = (wl[1] - wl[0]) / 2
        UVIS_wav = np.insert(wav_array, 0, wl_cent) / 10
        UVIS_waverr = np.insert(wav_err, 0, wl_err) / 10
        UVIS_precision_per_transit = np.insert(
            UVIS_spec_precision_per_transit, 0, UVIS_broadband_per_transit
        )
        UVIS_precision_total = np.insert(
            UVIS_spec_precision_total, 0, UVIS_broadband_total
        )

        all_data_to_save = np.array(
            [UVIS_wav, UVIS_waverr, UVIS_precision_per_transit, UVIS_precision_total]
        ).T

        file_save_name = "{}{}.UVIS_{}_sim.csv".format(
            data_folder, file_name[i], out_name
        )

        print(all_data_to_save.shape)
        np.savetxt(
            file_save_name,
            all_data_to_save,
            fmt=["%10.0f", "%10.5f", "%10.5f", "%10.5f"],
            header=header(
                i, wl, binlen, file_name, exposure_time, orbits_in_transit, no_of_transits
            ),
            delimiter=",",
            comments="# ",
        )

        print("Data has been saved as: ", file_save_name)
        print()


def read_UVIS(fpath):
    """Load data saved from `UVIS_simulation`.

    In additon to loading this data into a pandas DataFrame, this method will also
    automatically read the column names stored in the header.

    Args:
        fpath (string): Path to the csv file output by `UVIS_simulation`.

    Returns:
        df (DataFrame): Pandas DataFrame of data output by `UVIS_simulation`.

    """
    df = pd.read_csv(
        fpath,
        header=3,
        escapechar="#",
    )
    df.rename(columns=str.strip, inplace=True)
    return df


def plot_prec(ax, data, **kwargs):
    """Plot the calculated observation precision (in ppm) as a function of wavelength.

    Args:
        ax (matplotlib Axis): Axis to plot onto.
        data (DataFrame): Output data from `UVIS_simulation`.
        **kwargs: Additional options passed to plot command.

    Returns:
        ax (matplotlib Axis): Matplotlib axis after observation precision data plotted.

    """
    wav = data["Wavelength (nm)"]
    prec = data["program precision (ppm)"]

    ax.plot(wav[1:], prec[1:], drawstyle="steps-mid", **kwargs)


def plot_prec_wakeford(ax, data, kwargs_marg=None, kwargs_jit=None):
    """Plot the precision (in ppm) as a function of wavelength from Wakeford+2020.

    This will plot the curves resulting from both the Marginilization and Jitter method.

    Args:
        ax (matplotlib Axis): Axis to plot onto.
        data (DataFrame): Output data from `UVIS_simulation`.
        **kwargs: Additional options passed to plot command.

    Returns:
        ax (matplotlib Axis): Matplotlib axis after observation precision data plotted.

    """
    if kwargs_marg is None:
        kwargs_marg = {}
    if kwargs_jit is None:
        kwargs_jit = {}

    wav = data[:, 0]

    h41_depth_marg = data[:, 1]
    h41_deptherr_marg = data[:, 2] * 1e4  # convert from % to ppm

    h41_depth_jitter = data[:, 3]
    h41_deptherr_jitter = data[:, 4] * 1e4  # convert from % to ppm

    ax.plot(
        wav,
        h41_deptherr_marg,
        label="Marginalization method G280 (Wakeford+2020)",
        drawstyle="steps-mid",
        **kwargs_marg,
    )
    ax.plot(
        wav,
        h41_deptherr_jitter,
        label="Jitter method (Wakeford+2020)",
        drawstyle="steps-mid",
        **kwargs_jit,
    )

    return ax


if __name__ == "__main__":

    data_folder = ""
    file_name = ["H41_ETC"]

    exposure_time = [190]
    orbits_in_transit = [2]
    no_of_transits = [2]

    wl = [2000, 8000]
    startw = 2000
    endw = 8000
    binlen = 100

    out_name = "10nm_transmission"

    UVIS_simulation(
        data_folder,
        file_name,
        exposure_time,
        orbits_in_transit,
        no_of_transits,
        wl,
        startw,
        endw,
        binlen,
        out_name,
    )
