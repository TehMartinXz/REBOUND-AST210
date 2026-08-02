Quick disclaimer:
This README is basically just a copy-paste from the main comment/docstring in the monte_carlo.py file, but translated from Spanish (quick translation done with DeepL). Just in case anyone finds this repo useful or interesting at some point.

Most of the project was coded manually by me, except for the multi-threading usage, which was made with the assistance of AI.

You can also see our methods, results, and analysis in the [Proyecto_TRAPPIST_1_AST210.pdf file](https://github.com/TehMartinXz/REBOUND-AST210/blob/main/Proyecto_TRAPPIST_1_AST210.pdf)

# Description:
    This code uses the Monte Carlo method to determine the number
    of unstable configurations for TRAPPIST-1 by simulating the system's evolution
    over X years in REBOUND, verifying whether all the planets remain stable during
    this integration based on the variation of different orbital parameters within their
    uncertainty ranges.

    Instability criteria (only one needs to be met):
    - Distance between planets less than their mutual Hill radius.
    - Eccentricity greater than 1 (scattering).

    How it works:
    - We obtain the orbital parameters from the NASA Exoplanet Archive.
    - We list the orbital parameters of each planet along with their errors in a dictionary.
    - Then we have two options:
    * We vary only one orbital parameter (either for one of the bodies or for all the
    planets at once) within its error range.
    * We vary all orbital parameters of each planet randomly and uniformly
    within their error ranges.
    - The system is then simulated for X years using REBOUND. During the simulation,
    we verify whether the system remains stable or not based on instability criteria.
    - The process is repeated N times to generate statistics on stability based on
    different orbital parameters.

    Notes:
    - When a parameter is held constant, its reported nominal value will be used.
    - X will be 100 years by default; this is justified in the report.
    - Planetary parameters obtained from Agol et al. (2021) in the NASA Exoplanet Archive.
    - For eccentricity, which was not reported in Agol et al. (2021), the
    nominal value reported in Grimm et al. (2018) will be used when it remains constant. When it varies,
    it will be randomly generated between 0 and 0.9999. We hope this analysis will help narrow down the
    possible values for eccentricity.
    - According to the REBOUND documentation, Jacobi’s orbital elements are always used by default
    .
    - There are some parameters that we do not use directly, but we include them for completeness.
    For example, the semi-major axis (calculated from Kepler’s third law), or the eccentricity
    (calculated from ecosomega and esinomega).
    - We do not consider correlations between parameters, since the covariance matrix is not
    explicitly reported, and we were unable to download the Markov chain data
    to calculate it on our own. If we were able to obtain the MCMC results for the masses,
    periods, and t0, we could perform a more comprehensive analysis by considering these correlations.
    (The MCMC results by Agol et al. 2021 are the physical transit parameters.)
    - There are many conversions to float(). This is to prevent the data from being saved in
    NumPy format within the results files.

Translated with DeepL.com (free version)
