"""
    Este código se encarga de utilizar el método de Monte Carlo para obtener el número
    de configuraciones inestables de TRAPPIST-1 basándose en simular la evolución del sistema 
    por X años en REBOUND, validando si todos los planetas se mantienen estables durante 
    esta integración según la variación de diferentes parámetros orbitales dentro de sus 
    rangos de incertidumbre.

    Criterios de inestabilidad (sólo debe cumplirse uno):
    - Distancia entre los planetas menor a el radio de Hill.
    - Excentricidad mayor a 1 (scattering).

    Cómo funciona:
    - Obtenemos los parámetros orbitales desde el NASA Exoplanet Archive.
    - Enlistamos los parámetros orbitales de cada planeta junto con su error en un diccionario.
    - Luego tenemos 2 opciones:
    * Variamos solamente un parámetro orbital (ya sea para uno de los cuerpos, o todos los 
    planetas a la vez) dentro de su rango de error.
    * Variamos todos los parámetros orbitales de cada planeta de forma aleatoria y uniforme 
    dentro de sus rangos de error.
    - Después se integra el sistema por X años utilizando REBOUND. Durante la integración
    se valida si el sistema se mantiene estable o no según los criterios de inestabilidad.
    - Se repite el proceso N veces para poder hacer estadística sobre la estabilidad según
    diferentes parámetros orbitales.

    Notas:
    - Cuando un parámetro se mantenga fijo, se usará su valor nominal reportado.
    - X será 100 años por defecto; esto se justifica en el informe.
    - Parámetros de los planetas obtenidos de Agol et al. (2021) en el NASA Exoplanet Archive.
    - Para la excentricidad, que no se reportó en Agol et al. (2021), se utilizará el valor 
    nominal reportado en Grimm et al. 2018 cuando se mantenga fija. Cuando se varíe, 
    se hará aleatoriamente entre 0 y 0.9999. Esperamos que este análisis permita limitar los 
    posibles valores de la excentricidad.
    - Según la documentación de REBOUND, siempre se usan por defecto los elementos orbitales
    de Jacobi.
    - Hay algunos parámetros que no usamos directamente, pero los dejamos por completitud.
    Por ejemplo, el semieje mayor (calculado de la tercera Ley de Kepler), o la excentricidad
    (calculada a partir de ecosomega y esinomega).
    - No consideramos las correlaciones entre parámetros, ya que la matriz de covarianza no
    se reporta explícitamente, y no pudimos descargar los datos de las cadenas de Markov
    para calcularla por nuestra cuenta. Si logramos conseguir las MCMC para las masas,
    periodos, y t0, podríamos hacer un análisis más completo considerando estas correlaciones.
    (Los resultados de MCMC por Agol et al. 2021 son los parámetros de tránsito físicos).
    - Hay muchas conversiones a float(). Esto es para evitar que los datos se guarden en
    formato Numpy dentro de los archivos de resultados.
"""

# Incluimos las cotas superior e inferior de cada error (e_upper, e_lower) por consistencia, 
# ya que un par de parámetros tienen errores asimétricos (aunque la mayoría no).
# Creo que sólo el Radio tiene error asimétrico, y ni lo usamos en las simulaciones

planets = {
    'b': {
        'period': 1.510826,  # días (Agol et al. 2021)
        'e_period_upper': 0.000006,
        'e_period_lower': 0.000006,
        'e_period-theo_upper': 1.510826*2,  # valor exploratorio
        'e_period-theo_lower': 0.0151082*20, # valor exploratorio; 10 veces el dt de WHFast
        'a': 0.01154,  # AU (Agol et al. 2021)
        'e_a_upper': 0.00010,
        'e_a_lower': 0.00010,
        'radius': 1.116,  # radios terrestres (Agol et al. 2021)
        'e_radius_upper': 0.014,
        'e_radius_lower': 0.012,
        'mass': 1.374,  # masas terrestres (Agol et al. 2021)
        'e_mass_upper': 0.069,
        'e_mass_lower': 0.069,
        'e_mass-theo_upper': 1.374*2,
        'e_mass-theo_lower': 0,
        'eccentricity': 0.00622,  # valor nominal (Grimm et al. 2018)
        'e_eccentricity-theo_upper': 0.9999,  # valor teórico máximo
        'e_eccentricity-theo_lower': 0.0,  # valor teórico mínimo
        'e_eccentricity-obs_upper': 0.00304,  # error reportado en Grimm et al. 2018
        'e_eccentricity-obs_lower': 0.00304,
        'inclination': 89.728,  # grados (Agol et al. 2021)
        'e_inclination_upper': 0.165,
        'e_inclination_lower': 0.165,
        't0': 7257.55044, # días (Agol et al. 2021)
        'e_t0_upper': 0.00015,
        'e_t0_lower': 0.00015,
        'ecosomega': -0.00215, # Agol et al. 2021
        'e_ecosomega_upper': 0.00332,
        'e_ecosomega_lower': 0.00332,
        'esinomega': 0.00217, # Agol et al. 2021
        'e_esinomega_upper': 0.00244,
        'e_esinomega_lower': 0.00244
    },
    'c': {
        'period': 2.421937,
        'e_period_upper': 0.000018,
        'e_period_lower': 0.000018,
        'e_period-theo_upper': 2.421937*2,
        'e_period-theo_lower': 0.0151082*20,
        'a': 0.01580,
        'e_a_upper': 0.00013,
        'e_a_lower': 0.00013,
        'radius': 1.097,
        'e_radius_upper': 0.014,
        'e_radius_lower': 0.012,
        'mass': 1.308,
        'e_mass_upper': 0.056,
        'e_mass_lower': 0.056,
        'e_mass-theo_upper': 1.308*2,
        'e_mass-theo_lower': 0.0,
        'eccentricity': 0.00654,
        'e_eccentricity-theo_upper': 0.9999,
        'e_eccentricity-theo_lower': 0.0,
        'e_eccentricity-obs_upper': 0.00188,
        'e_eccentricity-obs_lower': 0.00188,
        'inclination': 89.778,
        'e_inclination_upper': 0.118,
        'e_inclination_lower': 0.118,
        't0': 7258.58728,
        'e_t0_upper': 0.00027,
        'e_t0_lower': 0.00027,
        'ecosomega': 0.00055,
        'e_ecosomega_upper': 0.00232,
        'e_ecosomega_lower': 0.00232,
        'esinomega': 0.00001,
        'e_esinomega_upper': 0.00171,
        'e_esinomega_lower': 0.00171
    },
    'd': {
        'period': 4.049219,
        'e_period_upper': 0.000026,
        'e_period_lower': 0.000026,
        'e_period-theo_upper': 4.049219*2,
        'e_period-theo_lower': 0.0151082*20,
        'a': 0.02227,
        'e_a_upper': 0.00019,
        'e_a_lower': 0.00019,
        'radius': 0.788,
        'e_radius_upper': 0.011,
        'e_radius_lower': 0.010,
        'mass': 0.388,
        'e_mass_upper': 0.012,
        'e_mass_lower': 0.012,
        'e_mass-theo_upper': 0.388*2,
        'e_mass-theo_lower': 0.0,
        'eccentricity': 0.00837,
        'e_eccentricity-theo_upper': 0.9999,
        'e_eccentricity-theo_lower': 0.0,
        'e_eccentricity-obs_upper': 0.00093,
        'e_eccentricity-obs_lower': 0.00093,
        'inclination': 89.896,
        'e_inclination_upper': 0.077,
        'e_inclination_lower': 0.077,
        't0': 7257.06768,
        'e_t0_upper': 0.00067,
        'e_t0_lower': 0.00067,
        'ecosomega': -0.00496,
        'e_ecosomega_upper': 0.00186,
        'e_ecosomega_lower': 0.00186,
        'esinomega': 0.00267,
        'e_esinomega_upper': 0.00112,
        'e_esinomega_lower': 0.00112
    },
    'e': {
        'period': 6.101013,
        'e_period_upper': 0.000035,
        'e_period_lower': 0.000035,
        'e_period-theo_upper': 6.101013*2,
        'e_period-theo_lower': 0.0151082*20,
        'a': 0.02925,
        'e_a_upper': 0.00250,
        'e_a_lower': 0.00250,
        'radius': 0.920,
        'e_radius_upper': 0.013,
        'e_radius_lower': 0.012,
        'mass': 0.692,
        'e_mass_upper': 0.022,
        'e_mass_lower': 0.022,
        'e_mass-theo_upper': 0.692*2,
        'e_mass-theo_lower': 0.0,
        'eccentricity': 0.00510,
        'e_eccentricity-theo_upper': 0.9999,
        'e_eccentricity-theo_lower': 0.0,
        'e_eccentricity-obs_upper': 0.00058,
        'e_eccentricity-obs_lower': 0.00058,
        'inclination': 89.793,
        'e_inclination_upper': 0.048,
        'e_inclination_lower': 0.048,
        't0': 7257.82771,
        'e_t0_upper': 0.00041,
        'e_t0_lower': 0.00041,
        'ecosomega': 0.00433,
        'e_ecosomega_upper': 0.00149,
        'e_ecosomega_lower': 0.00149,
        'esinomega': 0.00461,
        'e_esinomega_upper': 0.00087,
        'e_esinomega_lower': 0.00087
    },
    'f': {
        'period': 9.207540,
        'e_period_upper': 0.000032,
        'e_period_lower': 0.000032,
        'e_period-theo_upper': 9.207540*2,
        'e_period-theo_lower': 0.0151082*20,
        'a': 0.03849,
        'e_a_upper': 0.00033,
        'e_a_lower': 0.00033,
        'radius': 1.045,
        'e_radius_upper': 0.013,
        'e_radius_lower': 0.012,
        'mass': 1.039,
        'e_mass_upper': 0.031,
        'e_mass_lower': 0.031,
        'e_mass-theo_upper': 1.039*2,
        'e_mass-theo_lower': 0.0,
        'eccentricity': 0.01007,
        'e_eccentricity-theo_upper': 0.9999,
        'e_eccentricity-theo_lower': 0.0,
        'e_eccentricity-obs_upper': 0.00068,
        'e_eccentricity-obs_lower': 0.00068,
        'inclination': 89.740,
        'e_inclination_upper': 0.019,
        'e_inclination_lower': 0.019,
        't0': 7257.07426,
        'e_t0_upper': 0.00085,
        'e_t0_lower': 0.00085,
        'ecosomega': -0.00840,
        'e_ecosomega_upper': 0.00130,
        'e_ecosomega_lower': 0.00130,
        'esinomega': -0.00051,
        'e_esinomega_upper': 0.00087,
        'e_esinomega_lower': 0.00087
    },
    'g': {
        'period': 12.352446,
        'e_period_upper': 0.000054,
        'e_period_lower': 0.000054,
        'e_period-theo_upper': 12.352446*2,
        'e_period-theo_lower': 0.0151082*20,
        'a': 0.04683,
        'e_a_upper': 0.00040,
        'e_a_lower': 0.00040,
        'radius': 1.129,
        'e_radius_upper': 0.015,
        'e_radius_lower': 0.013,
        'mass': 1.321,
        'e_mass_upper': 0.038,
        'e_mass_lower': 0.038,
        'e_mass-theo_upper': 1.321*2,
        'e_mass-theo_lower': 0.0,
        'eccentricity': 0.00208,
        'e_eccentricity-theo_upper': 0.9999,
        'e_eccentricity-theo_lower': 0.0,
        'e_eccentricity-obs_upper': 0.00058,
        'e_eccentricity-obs_lower': 0.00058,
        'inclination': 89.742,
        'e_inclination_upper': 0.012,
        'e_inclination_lower': 0.012,
        't0': 7257.71462,
        'e_t0_upper': 0.00103,
        'e_t0_lower': 0.00103,
        'ecosomega': 0.00380,
        'e_ecosomega_upper': 0.00112,
        'e_ecosomega_lower': 0.00112,
        'esinomega': 0.00128,
        'e_esinomega_upper': 0.00070,
        'e_esinomega_lower': 0.00070
    },
    'h': {
        'period': 18.772866,
        'e_period_upper': 0.000214,
        'e_period_lower': 0.000214,
        'e_period-theo_upper': 18.772866*2,
        'e_period-theo_lower': 0.0151082*20,
        'a': 0.06189,
        'e_a_upper': 0.00053,
        'e_a_lower': 0.00053,
        'radius': 0.755,
        'e_radius_upper': 0.014,
        'e_radius_lower': 0.014,
        'mass': 0.326,
        'e_mass_upper': 0.020,
        'e_mass_lower': 0.020,
        'e_mass-theo_upper': 0.326*2,
        'e_mass-theo_lower': 0.0,
        'eccentricity': 0.00567,
        'e_eccentricity-theo_upper': 0.9999,
        'e_eccentricity-theo_lower': 0.0,
        'e_eccentricity-obs_upper': 0.00121,
        'e_eccentricity-obs_lower': 0.00121,
        'inclination': 89.805,
        'e_inclination_upper': 0.013,
        'e_inclination_lower': 0.013,
        't0': 7249.60676,
        'e_t0_upper': 0.00272,
        'e_t0_lower': 0.00272,
        'ecosomega': -0.00365,
        'e_ecosomega_upper': 0.00077,
        'e_ecosomega_lower': 0.00077,
        'esinomega': -0.00002,
        'e_esinomega_upper': 0.00044,
        'e_esinomega_lower': 0.00044
    }
}

star = { # Piaulet-Ghorayeb et al. 2025
    'radius': 0.1192,  # radios solares
    'e_radius_upper': 0.0013,
    'e_radius_lower': 0.0013,
    'mass': 0.0898,  # masas solares
    'e_mass_upper': 0.0023,
    'e_mass_lower': 0.0023
}

import rebound
import numpy as np
from astropy.constants import M_earth, M_sun
ratio = (M_earth / M_sun).value

def transit(P, t0, e_cos, e_sin):
    # Para simular correctamente TRAPPIST-1, necesitamos la posición angular de cada planeta en el momento inicial.
    # Esta función calcula todo lo necesario.
    # La excentricidad la obtenemos directamente desde e_cos y e_sin, que son los parámetros reportados en Agol et al. 2021,
    # en lugar de usar el valor nominal reportado en Grimm et al. 2018 (última vez que se reportó la excentricidad explícita).

    e = float(np.sqrt(e_cos**2 + e_sin**2))
    if e >= 1:
        e = 0.9999
    omega = np.arctan2(e_sin, e_cos)
    
    f_t = (np.pi / 2) - omega
    cos_ft = np.cos(f_t)
    denom = 1 + e*cos_ft

    E_t = np.arctan2(np.sqrt(1 - e**2) * np.sin(f_t)/denom, (e + cos_ft)/denom)
    # Despejamos la anomalía media a partir de la ecuación de Kepler
    M_t = E_t - e * np.sin(E_t)
    
    T_peri = t0 - (M_t/(2 * np.pi)) * P
    return e, omega, T_peri

def simulation(uniform = True, vary_planets = list, vary_star_mass = False, vary_params = list, sigma = 1, integrator = "whfast", param_value=None):
    # Si uniform = False, se hará variación normal (np.random.normal) en vez de uniforme (np.random.uniform).
    # Esto permite hacer una restricción de parámetros (variación uniforme), o un análisis estadístico (variación normal).
    # Cuando hacemos variación normal, se toma el promedio de los errores como sigma. Esto es correcto para la mayoría de
    # parámetros, aunque para algunos es una aproximación.

    # La opción param_value permite fijar un parámetro específico a un valor dado (cuando sólo se varía un parámetro)
    sim = rebound.Simulation()
    # Unidades del NASA Exoplanet Archive
    sim.units = ("days", "AU", "Msun")

    # Como esta función será llamada en iteraciones, necesitamos un integrador más rápido. WHFast es lo suficientemente
    # preciso para el análisis de estabilidad, ya que no nos importa cómo evoluciona el sistema después de que se
    # vuelve inestable.
    sim.integrator = integrator
    if integrator == "whfast":
        sim.dt = 0.0151082*2 # "Because WHFast is not an adaptive integrator, the user needs to set an appropriaye timestep" - Documentación
    # Este valor de 0.0151082 días se trata del 1% del menor periodo orbital posible de TRAPPIST-1b (considerando su máximo error).
    # Lo dejamos en (factor de 2) para duplicar la velocidad de la simulación, manteniendo fiabilidad física.
    # WHFast funciona calculando las interacciones entre planetas cada dt unidades de tiempo, 

    # Como usamos los parámetros reportados en Agol et al. 2021, dejamos el tiempo inicial en su mismo valor 
    # (explícito en el paper, descripción de la Tabla 2).
    sim.t = 7257.93115525

    # Para convertir las masas de los planetas de masas terrestres a masas solares, se multiplica por el ratio M_earth / M_sun.

    # ¿Variar masa de la estrella? (True/False)
    #vary_star_mass = False

    # Lista de parámetros orbitales a variar (planetas)
    # Parámetros posibles: 'period', 'eccentricity-theo', 'eccentricity-calc', 'inclination', 'mass', 'esinomega', 'ecosomega', 't0'
    # Notar que el semieje mayor lo calcula REBOUND a partir del periodo (que se mide por tránsito).
    # Sólo se puede variar una de las dos excentricidades (teórica y calculada) a la vez.
    # Para variar la excentricidad calculada, se debe variar al menos uno de los parámetros de los que depende 
    # (ecosomega, esinomega).
    #vary_params = ['period', 'eccentricity-calc', 'inclination', 'mass', 'esinomega', 'ecosomega', 't0']
    # Lista de planetas a variar
    # Planetas posibles: 'b', 'c', 'd', 'e', 'f', 'g', 'h'
    #vary_planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']

    params_used_star = []
    params_used_planets = []

    # Añadir estrella
    random_star_mass = 0
    if vary_star_mass:
        if param_value is not None:
            random_star_mass = float(param_value)
        if uniform and param_value is None and vary_params == ['star']:
            random_star_mass = np.abs(np.random.uniform(star['mass'] - star['e_mass_lower']*sigma, star['mass'] + star['e_mass_upper']*sigma))
        else:
            random_star_mass = np.abs(np.random.normal(star['mass'], (star['e_mass_upper'] + star['e_mass_lower'])*sigma / 2))
    else:
        # Valor nominal
        random_star_mass = star['mass']
    sim.add(m=random_star_mass)
    params_used_star.append(('mass', random_star_mass))

    # Falta calcular las posiciones angulares iniciales de cada planeta. Por suerte, en Agol et al. 2o21, se reporta el t0, ecosomega, y esinomega.

    # No modificar esta lista
    planet_names = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
    period_ratios = []
    # Añadir planetas
    for planet in planet_names:
        random_mass = planets[planet]['mass'] * ratio
        random_period = planets[planet]['period']
        random_inclination = float(np.radians(planets[planet]['inclination']))
        random_t0 = planets[planet]['t0']
        random_esinomega = planets[planet]['esinomega']
        random_ecosomega = planets[planet]['ecosomega']
        random_eccentricity = float(np.sqrt(random_esinomega**2 + random_ecosomega**2))
        if planet in vary_planets:
            # Los valores dentro de np.abs() son para evitar valores negativos en caso de que la variación normal se aleje mucho del valor nominal.
            if 'period' in vary_params:
                if param_value is not None:
                    random_period = float(param_value)
                elif uniform:
                    random_period = float(np.abs(np.random.uniform(planets[planet]['period'] - planets[planet]['e_period_lower']*sigma, planets[planet]['period'] + planets[planet]['e_period_upper']*sigma)))
                else:
                    random_period = float(np.abs(np.random.normal(planets[planet]['period'], (planets[planet]['e_period_upper'] + planets[planet]['e_period_lower'])*sigma / 2)))
            
            if 't0' in vary_params:
                if uniform:
                    random_t0 = float(np.abs(np.random.uniform(planets[planet]['t0'] - planets[planet]['e_t0_lower']*sigma, planets[planet]['t0'] + planets[planet]['e_t0_upper']*sigma)))
                else:
                    random_t0 = float(np.abs(np.random.normal(planets[planet]['t0'], (planets[planet]['e_t0_upper'] + planets[planet]['e_t0_lower'])*sigma / 2)))

            if 'esinomega' in vary_params:
                if uniform:
                    random_esinomega = np.random.uniform(planets[planet]['esinomega'] - planets[planet]['e_esinomega_lower']*sigma, planets[planet]['esinomega'] + planets[planet]['e_esinomega_upper']*sigma)
                else:
                    random_esinomega = np.random.normal(planets[planet]['esinomega'], (planets[planet]['e_esinomega_upper'] + planets[planet]['e_esinomega_lower'])*sigma / 2)

            if 'ecosomega' in vary_params:
                if uniform:
                    random_ecosomega = np.random.uniform(planets[planet]['ecosomega'] - planets[planet]['e_ecosomega_lower']*sigma, planets[planet]['ecosomega'] + planets[planet]['e_ecosomega_upper']*sigma)
                else:
                    random_ecosomega = np.random.normal(planets[planet]['ecosomega'], (planets[planet]['e_ecosomega_upper'] + planets[planet]['e_ecosomega_lower'])*sigma / 2)

            if 'eccentricity-theo' in vary_params:
                if param_value is not None:
                    random_eccentricity = float(param_value)
                else:
                    random_eccentricity = np.random.uniform(planets[planet]['e_eccentricity-theo_lower'], planets[planet]['e_eccentricity-theo_upper'])

                # Cuando usamos excentricidad aleatoria independiente de los otros parámetros, es necesario recalcular esinomega y ecosomega 
                # para que sean consistentes con esta excentricidad.
                transit_angle = np.arctan2(random_esinomega, random_ecosomega)
                random_ecosomega = random_eccentricity * np.cos(transit_angle)
                random_esinomega = random_eccentricity * np.sin(transit_angle)
            elif 'eccentricity-calc' in vary_params:
                # Esto ya tiene en cuenta si es variación uniforme o normal, según los otros parámetros.
                #random_eccentricity = transit(random_period, random_t0, random_ecosomega, random_esinomega)[0]
                () # No es necesario calcularla aquí, ya que se pide más abajo de todas formas.

            if 'inclination-theo' in vary_params:
                # Esto lo añadimos para experimentar si la co-planaridad es requisito de estabilidad, o es simplemente un producto de la formación del sistema
                # Ojo que como las proyecciones sin(i) cos(i) no son lineales, variamos el seno/coseno de la inclinación, y luego obtenemos la inclinación a partir de eso.
                # Lo variamos sólo de forma uniforme.
                # No es compatible con variar 'inclination' a la vez.
                # No lo añadí a los parámetros de los planetas porque... si
                random_inclination = float(np.arccos(np.random.uniform(-1, 1)))

            if 'period-theo' in vary_params:
                # Experimento de variación semi-aleatoria del periodo
                # No es compatible con variar 'period' a la vez.
                random_period = float(np.abs(np.random.uniform(planets[planet]['e_period-theo_lower'], planets[planet]['e_period-theo_upper'])))

            if 'mass-theo' in vary_params:
                # Experimento de variación semi-aleatoria de la masa
                # No es compatible con variar 'mass' a la vez.
                random_mass = float(np.abs(np.random.uniform(planets[planet]['e_mass-theo_lower'] * ratio, planets[planet]['e_mass-theo_upper'] * ratio)))

            if 'inclination' in vary_params:
                if param_value is not None:
                    random_inclination = float(np.radians(param_value))
                elif uniform:
                    random_inclination = float(np.radians(np.random.uniform(planets[planet]['inclination'] - planets[planet]['e_inclination_lower']*sigma, planets[planet]['inclination'] + planets[planet]['e_inclination_upper']*sigma)))
                else:
                    random_inclination = float(np.radians(np.random.normal(planets[planet]['inclination'], (planets[planet]['e_inclination_upper'] + planets[planet]['e_inclination_lower'])*sigma / 2)))
            
            if 'mass' in vary_params:
                if param_value is not None:
                    random_mass = float(param_value) * ratio
                elif uniform:
                    random_mass = float(np.abs(np.random.uniform((planets[planet]['mass'] - planets[planet]['e_mass_lower']*sigma) * ratio, (planets[planet]['mass'] + planets[planet]['e_mass_upper']*sigma) * ratio)))
                else:
                    random_mass = float(np.abs(np.random.normal(planets[planet]['mass'], (planets[planet]['e_mass_upper'] + planets[planet]['e_mass_lower'])*sigma / 2))) * ratio

            # Parámetros calculados:
            e_calc, omega_calc, time_peri_calc = transit(random_period, random_t0, random_ecosomega, random_esinomega)
            sim.add(m=random_mass, e=e_calc, omega=omega_calc, T=time_peri_calc, inc=random_inclination, P=random_period, name=planet)
            params_used_planets.append([planet, float(random_mass), float(random_period), float(e_calc), float(random_inclination), float(random_esinomega), float(random_ecosomega), float(random_t0)])

        else:
            # Parámetros calculados con valores nominales:
            eccentricity_calc, omega_calc, time_peri_calc = transit(planets[planet]['period'], planets[planet]['t0'], planets[planet]['ecosomega'], planets[planet]['esinomega'])
            sim.add(m=planets[planet]['mass'] * ratio, e=eccentricity_calc, omega=omega_calc, T=time_peri_calc, inc=float(np.radians(planets[planet]['inclination'])), P=planets[planet]['period'], name=planet)
            params_used_planets.append([planet, float(planets[planet]['mass'] * ratio), float(planets[planet]['period']), float(eccentricity_calc), float(np.radians(planets[planet]['inclination'])), float(planets[planet]['esinomega']), float(planets[planet]['ecosomega']), float(planets[planet]['t0'])])
    
    sim.move_to_com()

    # Detener la simulación si se detecta inestabilidad:

    # Radio de Hill: - Lo validamos en la ejecución de la simulación, no acá
    #sim.exit_min_distance = min_hill_radius
    # 100 veces el apoastro del planeta más lejano
    sim.exit_max_distance = 100 * planets['h']['a']*(1 + random_eccentricity)
        

    return sim, params_used_planets, random_star_mass


# Loop principal (aplicando la idea de Monte Carlo):

# Disclaimer: La lógica del Loop paralelizado (sólo el uso de joblib Parallel) y la llamada de algunas funciones (para facilitar el trabajo) lo hicimos con ayuda de Claude AI (Sonnet 4.6). Todo el resto del código está hecho a mano.
from joblib import Parallel, delayed

def run_simulation(i, uniform=False, vary_planets1=list, vary_star_mass1=False, vary_params1=list, sigma1=1):
    # Ejecutar la simulación y obtener los parámetros utilizados para esta iteración
    sim, params_planets, random_star_mass = simulation(uniform, vary_planets=vary_planets1, vary_star_mass=vary_star_mass1, vary_params=vary_params1, sigma=sigma1)
    
    #sim.ri_whfast.safe_mode = 0 # prueba - Resulta ser que no ofrece mejora de rendimiento notable, a si que lo dejamos por defecto

    # Calcular radio de Hill mínimo (se calcula en pares: b-c, c-d, d-e, e-f, f-g, g-h)
    # Aprovechamos también a calcular las razones entre los periodos orbitales de planetas vecinos.
    # Ojo que hay que usar los valores aleatorios, no los nominales.
    min_hill_radius = float('inf')
    period_ratios = []
    for planet in ['b', 'c', 'd', 'e', 'f', 'g']: # No incluimos a TRAPPIST-1h porque no tiene un planeta exterior con el que comparar
        m1 = planets[planet]['mass'] * ratio
        m2 = planets[chr(ord(planet) + 1)]['mass'] * ratio
        # Semiejes desde la tercera Ley de Kepler
        p1 = planets[planet]['period']
        p2 = planets[chr(ord(planet) + 1)]['period']
        a1 = (sim.G * (random_star_mass + m1) * (p1 / (2 * np.pi))**2)**(1/3)
        a2 = (sim.G * (random_star_mass + m2) * (p2 / (2 * np.pi))**2)**(1/3)
        Mstar = random_star_mass
        if planet in vary_planets1:
            m1 = params_planets[[p[0] for p in params_planets].index(planet)][1]
            p1 = params_planets[[p[0] for p in params_planets].index(planet)][2]
            a1 = (sim.G * (random_star_mass + m1) * (p1 / (2 * np.pi))**2)**(1/3)
            if chr(ord(planet) + 1) in vary_planets1:
                m2 = params_planets[[p[0] for p in params_planets].index(chr(ord(planet) + 1))][1]
                p2 = params_planets[[p[0] for p in params_planets].index(chr(ord(planet) + 1))][2]
                a2 = (sim.G * (random_star_mass + m2) * (p2 / (2 * np.pi))**2)**(1/3)
        hill_radius = (a1 + a2)/2 * ((m1 + m2) / (3 * Mstar))**(1/3)
        if hill_radius < min_hill_radius:
            min_hill_radius = hill_radius
        # Y la razón entre los periodos orbitales:
        period_ratios.append((i, planet + '+' + chr(ord(planet) + 1), float(p2) / float(p1)))
    
    sim.exit_min_distance = min_hill_radius

    DIAS_SIMULACION = 100 * 365.25
    sobrevivio = True
    t_colision = np.nan
    inestabilidad_tipo = 'estable'
    
    try:
        sim.integrate(DIAS_SIMULACION)
    except rebound.Encounter:
        sobrevivio = False
        inestabilidad_tipo = 'colision'
        t_colision = sim.t - 7257.93115525 # Tiempo exacto en el que chocaron
    except rebound.Escape:
        sobrevivio = False
        inestabilidad_tipo = 'escape'
        t_colision = sim.t - 7257.93115525 # Tiempo exacto en el que escaparon
    
    return {
        'iter': i,
        'sigma': sigma1,
        'sobrevivio': sobrevivio,
        'inestabilidad_tipo': inestabilidad_tipo,
        't_colision_dias': t_colision, # Sólo si hubo colisión o escape
        'random_star_mass': random_star_mass,
        'period_ratios': period_ratios,
        'params_planets': params_planets,
    }

# Por razones fuera de mi comprensión, ponemos esta redundancia
if __name__ == '__main__':
    import pandas as pd
    import time
    from pathlib import Path

    def main(planetas_a_variar2=list, parametros_a_variar2=list, variar_masa_estrella2=False, n_sigmas2=1, uniform2=False, N2=10, folder_prefix=''):
        # Parámetros para variar:
        planetas_a_variar = planetas_a_variar2 #['b', 'c', 'd', 'e', 'f', 'g', 'h'] # Elegir los planetas cuyos parámetros variarán
        parametros_a_variar = parametros_a_variar2 #['period', 'eccentricity-calc', 'inclination', 'mass', 'esinomega', 'ecosomega', 't0'] # Elegir los parámetros orbitales a variar
        variar_masa_estrella = variar_masa_estrella2 # Elegir si se variará la masa de la estrella (True/False)
        n_sigmas = n_sigmas2 # Número de sigmas para considerar el error
        uniforme = uniform2 # Elegir si se usará distribución uniforme

        
        # Número de simulaciones
        N = N2

        print(f"Iniciando Monte Carlo: {N} simulaciones de 100 años...")
        
        # n_jobs=-1 usa todos los hilos
        tiempo_inicio = time.perf_counter()
        # n_jobs=12 para jugar CS mientras se ejecuta
        resultados = Parallel(n_jobs=-1, verbose=10)(
            delayed(run_simulation)(i, uniform=uniforme, vary_planets1=planetas_a_variar, vary_star_mass1=variar_masa_estrella, vary_params1=parametros_a_variar, sigma1=n_sigmas) for i in range(N)
        )
        tiempo_fin = time.perf_counter()
        print(f"Monte Carlo completado en {(tiempo_fin - tiempo_inicio) / 60:.2f} minutos.")
        
        datos_estabilidad = {}
        for resultado in resultados:
            for key, value in resultado.items():
                if key != 'params_planets' and key != 'period_ratios':
                    if key not in datos_estabilidad:
                        datos_estabilidad[key] = []
                    datos_estabilidad[key].append(value)

        datos_param_planets = {'iter': []}  # Inicializamos con la columna de iteración
        for resultado in resultados:
            datos_param_planets['iter'].append(resultado['iter'])
            if 'params_planets' in resultado:
                for param in resultado['params_planets']:
                    if param[0] not in datos_param_planets:
                        datos_param_planets[param[0]] = []
                    datos_param_planets[param[0]].append(param[1:])

        datos_period_ratios = {'iter': []}  # Inicializamos con la columna de iteración
        for resultado in resultados:
            datos_period_ratios['iter'].append(resultado['iter'])
            if 'period_ratios' in resultado:
                for ratio in resultado['period_ratios']:
                    # ratio tiene la estructura: (iter, 'par_de_planetas', valor_razon)
                    if ratio[1] not in datos_period_ratios:
                        datos_period_ratios[ratio[1]] = []
                    datos_period_ratios[ratio[1]].append(float(ratio[2]))

        carpeta = (folder_prefix + time.strftime("%Y-%m-%d %X") + '_resultados').replace(":", "-")
        Path('Results/' + carpeta).mkdir(exist_ok=True)
        df = pd.DataFrame(datos_estabilidad)
        df2 = pd.DataFrame(datos_param_planets)
        df3 = pd.DataFrame(datos_period_ratios)
        df.to_csv(f'Results/{carpeta}/estabilidad_resultados.csv', index=False)
        df2.to_csv(f'Results/{carpeta}/parametros_planetas.csv', index=False)
        df3.to_csv(f'Results/{carpeta}/period_ratios.csv', index=False)
        

        print(df)
        
        # Resumen rápido
        vivos = df["sobrevivio"].sum()
        muertos = len(df) - vivos
        print("\n=== RESUMEN MONTE CARLO ===")
        print(folder_prefix)
        print(f"Sobrevivieron : {vivos} ({vivos/len(df)*100:.1f}%)")
        print(f"Inestables    : {muertos} ({muertos/len(df)*100:.1f}%)")
        print('Inestabilidad por tipo:')
        print(df['inestabilidad_tipo'].value_counts())
        return

    # Dejar ejecutando diferentes configuraciones de Monte Carlo para poder dejarlo durante la noche
    # Sí, el código no es muy bonito (tiene muchas redundancias), pero funciona :D
    # Y perdón por algunos de los comentarios de relleno. Así programo yo.... hmm, tengo hambre

    todos_los_planetas = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
    todos_los_parametros = ['period', 'eccentricity-calc', 'inclination', 'mass', 'esinomega', 'ecosomega', 't0']

    """
        Pruebas de simulación única a largo plazo
    """
    def unica(years=100, vary_planets2=[], vary_star_mass2=False, vary_params2=[], sigma2=1, integrator2="whfast", param_value2=None):
        sim, params_planets, random_star_mass = simulation(uniform=False, vary_planets=vary_planets2, vary_star_mass=vary_star_mass2, vary_params=vary_params2, sigma=sigma2, integrator=integrator2, param_value=param_value2)

        min_hill_radius = float('inf')
        period_ratios = []
        for planet in ['b', 'c', 'd', 'e', 'f', 'g']: # No incluimos a TRAPPIST-1h porque no tiene un planeta exterior con el que comparar
            m1 = planets[planet]['mass'] * ratio
            m2 = planets[chr(ord(planet) + 1)]['mass'] * ratio
            # Semiejes desde la tercera Ley de Kepler
            p1 = planets[planet]['period']
            p2 = planets[chr(ord(planet) + 1)]['period']
            a1 = (sim.G * (random_star_mass + m1) * (p1 / (2 * np.pi))**2)**(1/3)
            a2 = (sim.G * (random_star_mass + m2) * (p2 / (2 * np.pi))**2)**(1/3)
            Mstar = random_star_mass

            hill_radius = (a1 + a2)/2 * ((m1 + m2) / (3 * Mstar))**(1/3)
            if hill_radius < min_hill_radius:
                min_hill_radius = hill_radius
            # Y la razón entre los periodos orbitales:
            period_ratios.append((planet + '+' + chr(ord(planet) + 1), float(p2) / float(p1)))
        
        sim.exit_min_distance = min_hill_radius

        DIAS_SIMULACION = years * 365.25
        sobrevivio = True
        t_colision = np.nan
        inestabilidad_tipo = 'estable'

        culpable = None
        E_inicial = np.nan
        E_final = np.nan
        

        e_finales_dict = {f"e_{p}_final": np.nan for p in ['b', 'c', 'd', 'e', 'f', 'g', 'h']}
        a_finales_dict = {f"a_{p}_final": np.nan for p in ['b', 'c', 'd', 'e', 'f', 'g', 'h']}
        try:
            E_inicial = sim.energy()
            sim.integrate(DIAS_SIMULACION)
            E_final = sim.energy()

            # Datos finales
            for idx, p_name in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h'], start=1):
                e_finales_dict[f"e_{p_name}_final"] = float(sim.particles[idx].e)
                a_finales_dict[f"a_{p_name}_final"] = float(sim.particles[idx].a)
        except rebound.Encounter:
            sobrevivio = False
            inestabilidad_tipo = 'colision'
            t_colision = sim.t - 7257.93115525 # Tiempo exacto en el que chocaron

            min_dist = float('inf')
            for idx in range(1, sim.N - 1):
                p1 = sim.particles[idx]
                p2 = sim.particles[idx+1]
                dist = np.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)
                if dist < min_dist:
                    min_dist = dist
                    culpable = f"{p1.name}+{p2.name}" # Ej: 'b+c' o 'f+g'
        except rebound.Escape:
            sobrevivio = False
            inestabilidad_tipo = 'escape'
            t_colision = sim.t - 7257.93115525 # Tiempo exacto en el que escaparon

            max_dist = 0.0
            for idx in range(1, sim.N):
                p = sim.particles[idx]
                dist = np.sqrt(p.x**2 + p.y**2 + p.z**2)
                if dist >= max_dist:
                    max_dist = dist
                    culpable = p.name
        
        return {
            'sobrevivio': sobrevivio,
            't_integración_yr' : years,
            'inestabilidad_tipo': inestabilidad_tipo,
            't_colision_dias': t_colision, # Sólo si hubo colisión o escape
            'random_star_mass': random_star_mass,
            'culpable_colision': culpable if not sobrevivio and inestabilidad_tipo == 'colision' else None,
            'culpable_escape': culpable if not sobrevivio and inestabilidad_tipo == 'escape' else None,
            'energia_inicial': E_inicial,
            'energia_final': E_final,
            'e_finales': e_finales_dict,
            'a_finales': a_finales_dict,
            'period_ratios': period_ratios,
            'params_planets': params_planets,
        }
