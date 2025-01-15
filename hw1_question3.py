import numpy as np
from astropy import units
from astropy.constants import h, c, k_B, R_sun, R_earth, au
import matplotlib.pyplot as plt


# Part a -> set up a blackbody calcuation function
# function should take wavelength or frequency and temperature and compute the flux density (spectrum)
# plot the blackbody spectrum of the Sun as observed at 10 pc using logarithmic x and y scale

def bb_func(wavelength, temperature, radius, distance):
    h_cgs = h.to(units.erg * units.s)
    c_cgs = c.to(units.cm / units.s)
    kB_cgs = k_B.to(units.erg / units.K)

    constants = 2 * h_cgs * (c_cgs ** 2)
    frac = radius**2/distance**2
    return constants * frac * (1 / wavelength ** 5) \
        * (1/(np.exp((h_cgs*c_cgs)/(wavelength*temperature*kB_cgs)) - 1))

def part_a(wavelength, temperature, radius, distance, title):
    '''

    :param wavelength:
    :param temperature:
    :return:
    '''
    # Equation 1.52 from our textbook
    bb_parta = bb_func(wavelength, temperature, radius, distance)
    print('bb_func: ', bb_parta)

    # Plotting
    fig, ax = plt.subplots()
    ax.plot(wavelength, bb_parta)
    ax.invert_xaxis()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Wavelength (cm)')
    ax.set_ylabel('Flux Density (erg/(s*cm^3)')
    ax.set_title(title)
    plt.show()


def part_b_one(wavelength, temperatures, radius, distance, titles):
    '''

    :param wavelengths:
    :param temperatures:
    :param titles:
    :return:
    '''
    # bullet point one T = 5800K, M5 dwarf companion T = 2800K
    # bullet point two T = 2800K and from homework T = 10000K WD companion
    # bullet point three T = 4290K and from homework neutron star = 1000000K

    loop1 = temperatures.shape[0]
    loop2 = temperatures.shape[1]
    print("Loop sizes: ", loop1, loop2)

    for ii in range(loop1):
        print(ii)
        radius1 = radius[ii, 0]*units.pc
        radius2 = radius[ii, 1]*units.pc
        fig, ax = plt.subplots()
        temp1 = temperatures[ii, 0] * units.K
        temp2 = temperatures[ii, 1] * units.K

        bb_temp1 = bb_func(wavelength, temp1, radius1, distance)
        bb_temp2 = bb_func(wavelength, temp2, radius2, distance)

        sum_bb_spectrum = bb_temp1 + bb_temp2

        ax.plot(wavelength, sum_bb_spectrum)
        ax.plot(wavelength, bb_temp1, 'r', ls='dashed')
        ax.plot(wavelength, bb_temp2, 'g', ls='dashed')
        ax.invert_xaxis()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Wavelength (cm)')
        ax.set_ylabel('Flux Density (erg/(s*cm^3)')
        ax.set_title(titles[ii])
        fig.savefig(f"{titles[ii]}")
        fig.clear()

    #plt.show()

def part_c(wavelength, temperature, radius, distance, title):
    #print(wavelength, temperature, radius, distance)
    bb_partc = bb_func(wavelength, temperature, radius, distance)
    c_cgs = c.to(units.cm / units.s)
    tau = c_cgs/((3*10**14*units.Hz)*wavelength)
    print(tau)
    pure_absorption = bb_partc*np.exp(-tau)


wavelength = np.arange(0.1,100, 0.1) * units.um
wavelength_cm = wavelength.to(units.cm)
temperature_a = 5800 *units.K   #citation needed
radius_sun = R_sun.to(units.pc)
radius_wd = R_earth.to(units.pc)
radius_arcturus = 25*R_sun.to(units.pc)
radius_ns_org = 5*units.km
radius_ns = radius_ns_org.to(units.pc)
distance = 10*units.pc
#part_a(wavelength_cm, temperature_a, radius_sun, distance, title="Part a")

temperature_b = np.array([[5800, 2800],
                 [2800, 10000],
                 [4290, 1000000]])
radius_b = np.array([[radius_sun.value, 0.1*radius_sun.value],
                     [0.1*radius_sun.value, radius_wd.value],
                     [radius_arcturus.value, radius_ns.value]])
#part_b_one(wavelength_cm, temperature_b, radius_b, distance, titles=["Bullet Point 1", "Bullet Point 2", "Bullet Point 3"])

temp_c = 3000*units.K
radius_c = 5*radius_sun
distance_c_au = 100*units.au
distance_c = distance_c_au.to(units.pc)

part_c(wavelength_cm, temp_c, radius_c, distance_c, title="Part c")