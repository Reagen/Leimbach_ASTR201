import numpy as np
from astropy import units
from astropy.constants import h, c, k_B, R_sun, R_earth, au, sigma_sb
import matplotlib.pyplot as plt
from scipy.integrate import quad


# Part a -> set up a blackbody calcuation function
# function should take wavelength or frequency and temperature and compute the flux density (spectrum)
# plot the blackbody spectrum of the Sun as observed at 10 pc using logarithmic x and y scale

def bb_func(wavelength, temperature, radius, distance):
    '''
    For an array of wavelengths, a given temperature, the radius of the object, and the distance you are away
    return the blackbody function
    :param wavelength: array in units cm
    :param temperature: number in units of Kelvin
    :param radius: number must have same units as distance
    :param distance: number must have same units as radius
    :return: blackbody function Equation 1.52 from our textbook
    '''
    # Constants (found in textbook)
    h_cgs = h.to(units.erg * units.s)
    c_cgs = c.to(units.cm / units.s)
    kB_cgs = k_B.to(units.erg / units.K)

    constants = 2 * h_cgs * (c_cgs ** 2)
    frac = radius**2/distance**2     # for testing *(units.pc)**2
    print(constants * frac * (1 / wavelength ** 5) \
        * (1/(np.exp((h_cgs*c_cgs)/(wavelength*temperature*kB_cgs)) - 1)))
    return constants * frac * (1 / wavelength ** 5) \
        * (1/(np.exp((h_cgs*c_cgs)/(wavelength*temperature*kB_cgs)) - 1))

def bb_func_part4(wavelength, temperature, radius, distance, absorption):
    h_cgs = h.to(units.erg * units.s)
    c_cgs = c.to(units.cm / units.s)
    kB_cgs = k_B.to(units.erg / units.K)

    constants = 2 * h_cgs.value * (c_cgs.value ** 2)
    print('Wavelength', wavelength)
    frac = radius ** 2 / distance ** 2  # for testing *(units.pc)**2
    bbfunction = (constants * frac.value * (1 / wavelength ** 5) * (1 / (np.exp((h_cgs.value * c_cgs.value) / (wavelength * temperature.value * kB_cgs.value)) - 1)))
    return bbfunction*absorption

def part_a(wavelength, temperature, radius, distance, title):
    '''
    Plotting blackbody spectrum of the Sun as observed at 10pc.
    :param wavelength: array in units cm
    :param temperature: number in units of Kelvin
    :param radius: number must have same units as distance
    :param distance: number must have same units as radius
    :param title: str title of the plot
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

def part_b(wavelength, temperatures, radius, distance, titles):
    '''
    Plotting the individual and combined blackbody spectra for combinations of binaries, observed at 10ppc
    # bullet point one T = 5800K, M5 dwarf companion T = 2800K
    # bullet point two T = 2800K and from homework T = 10000K WD companion
    # bullet point three T = 4290K and from homework neutron star = 1000000K

    :param wavelength: array in units cm
    :param temperature: number in units of Kelvin
    :param radius: number must have same units as distance
    :param distance: number must have same units as radius
    :param title: str title of the plot
    :return:
    '''

    loop1 = temperatures.shape[0]
    #loop2 = temperatures.shape[1]
    #print("Loop sizes: ", loop1, loop2)

    #Looping over the three different combinations
    for ii in range(loop1):
        print(ii)
        radius1 = radius[ii, 0]*units.pc
        radius2 = radius[ii, 1]*units.pc
        fig, ax = plt.subplots()
        temp1 = temperatures[ii, 0] * units.K
        temp2 = temperatures[ii, 1] * units.K

        #Calculate the individual Blackbody functions
        bb_temp1 = bb_func(wavelength, temp1, radius1, distance)
        bb_temp2 = bb_func(wavelength, temp2, radius2, distance)

        #Sum the blackbody functions to see the combined function
        sum_bb_spectrum = bb_temp1 + bb_temp2

        # Plotting
        ax.plot(wavelength, sum_bb_spectrum, label = "Sum of Blackbody functions")
        ax.plot(wavelength, bb_temp1, 'r', ls='dashed', label="First blackbody function")
        ax.plot(wavelength, bb_temp2, 'g', ls='dashed', label="Second blackbody function")
        ax.invert_xaxis()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Wavelength (cm)')
        ax.set_ylabel('Flux Density (erg/(s*cm^3)')
        ax.set_title(titles[ii])
        ax.legend()
        fig.savefig(f"{titles[ii]}")
        fig.clear()

    #plt.show()

def part_cd(wavelength, temperature, radius, distance, title):
    '''

    :param wavelength:
    :param temperature:
    :param radius:
    :param distance:
    :param title:
    :return:
    '''
    #print(wavelength, temperature, radius, distance)
    bb_partc = bb_func(wavelength, temperature, radius, distance)
    c_cgs = c.to(units.cm / units.s)
    sb_const = sigma_sb.to(units.erg/units.cm**2/units.s/units.K**4)
    #print(sb_const)
    tau = c_cgs/((3*10**14*units.Hz)*wavelength)
    #print('tau', tau)
    absorption_part = (1-np.exp(-tau))

    integrate = quad(bb_func_part4, 0, np.inf, args=(temperature.value, radius.value, distance.value, absorption_part))

    #flux = np.trapz(integrand)
    print(integrate)
    temp = (integrate/(2*sb_const))**(-1/4)
    print(temp)


# Constants and defined variables
wavelength = np.arange(0.1,100, 0.1) * units.um
wavelength_cm = wavelength.to(units.cm)
temperature_a = 5800 *units.K   #citation (https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
# The radius of the sun and the earth are given constants in astropy and I used the citation link to check
radius_sun = R_sun.to(units.pc)  #citation (https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
# We learned in our Survey of Astronomy class that a White Dwarf is about the size of the Earth
radius_wd = R_earth.to(units.pc) # citation (https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
radius_arcturus = 25*R_sun.to(units.pc) #citation (https://skyandtelescope.org/astronomy-news/meet-arcturus-guardian/)
radius_ns_org = 5*units.km # We learned in our Survey of Astronomy class that a Neutron star is ~10km in diameter
radius_ns = radius_ns_org.to(units.pc)
distance = 10*units.pc
temp_c = 3000*units.K
radius_c = 5*radius_sun
distance_c_au = 100*units.au
distance_c = distance_c_au.to(units.pc)
#Citations for unknown temperatures given in write-up
temperature_b = np.array([[5800, 2800],
                 [2800, 10000],
                 [4290, 1000000]])
#Citations for unknown radi given in write-up
radius_b = np.array([[radius_sun.value, 0.1*radius_sun.value],
                     [0.1*radius_sun.value, radius_wd.value],
                     [radius_arcturus.value, radius_ns.value]])


#part_a(wavelength_cm, temperature_a, radius_sun, distance, title="Part a")
#part_b(wavelength_cm, temperature_b, radius_b, distance, titles=["Bullet Point 1", "Bullet Point 2", "Bullet Point 3"])
part_cd(wavelength_cm, temp_c, radius_c, distance_c, title="Part c")


'''wavelength2 = np.arange(0.0001,100000, 0.1) * units.um
wavelength_cm2 = wavelength.to(units.cm)

part_a(wavelength_cm2, 1000000*units.K, radius_ns.value, distance, title="test")'''