""" A fantastic python code to determine the quenched SFH parameters of galaxies using emcee (http://dan.iel.fm/emcee/current/). This file contains all the functions needed to determine the mean SFH parameters of a population.
    
    N.B. The data files .ised_ASCII contain the extracted bc03 models and have a 0 in the origin at [0,0]. The first row contains the model ages (from the second column) - data[0,1:]. The first column contains the model lambda values (from the second row) - data[1:,0]. The remaining data[1:,1:] are the flux values at each of the ages (columns, x) and lambda (rows, y) values 
    """

import numpy as N
import scipy as S
import pylab as P
#import pyfits as F
import astropy.io.fits as F
from scipy.io.idl import readsav
#import pyfits as F
import emcee
#import triangle  # the new version is below
import corner
import time
import os
import matplotlib.image as mpimg
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import kde
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp2d
from itertools import product
import sys

cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)

font = {'family':'serif', 'size':16}
P.rc('font', **font)
P.rc('xtick', labelsize='medium')
P.rc('ytick', labelsize='medium')
P.rc('axes', labelsize='medium')




'''
# defaults
tq = N.linspace(0.003, 13.8, 100)
tau = N.linspace(0.003, 4, 100)
ages = N.linspace(10.88861228, 13.67023409, 50)
use_table = False

paramfile = 'posterior_params.in'
try:
    with open(paramfile) as f:
        for line in f:
            arg = line.rstrip('\n').lower().strip(' ').split('=')

            if (arg[0] in ['lookup', 'lookups', 'lookuptable', 'lookuptables', 'lu', 'lut']):
                use_table = True
                tables = arg[1].split(',')

                if tables[0] in ['default']:
                    col1_file  = 'nuv_look_up_ssfr.npy'
                    col2_file  = 'ur_look_up_ssfr.npy'
                    lu_default = True
                else:
                    if len(tables) < 2:
                        tables = arg[1].split(' ')
                        if len(tables) < 2:
                            print('Error: if not "default", 2 lookup tables needed, filenames separated by "," or " "')
                            exit(-1)
                    col1_file  = tables[0]
                    col2_file  = tables[1]
                    lu_default = False

            elif (arg[0] in ['tq', 't_q', 'tquench', 't_quench', 'quench_time', 'quenching_time']):
                valstr = arg[1].split(',')
                if len(valstr) != 3:
                    print('Error: if specifying quenching times tq in %s, must be formatted "tq = tstart, tend, n_vals"' % paramfile)
                    exit(-1)
                    
                tq = N.linspace(float(valstr[0]), float(valstr[1]), int(valstr[2]))

            elif (arg[0] in ['tau', 'exptau', 'quenchingrate', 'quenching_rate']):
                valstr = arg[1].split(',')
                if len(valstr) != 3:
                    print('Error: if specifying quenching rates tau in %s, must be formatted "tau = taustart, tauend, n_vals"' % paramfile)
                    exit(-1)
                    
                tau = N.linspace(float(valstr[0]), float(valstr[1]), int(valstr[2]))

            elif (arg[0] in ['ages', 'age', 't_obs', 'age_obs']):
                valstr = arg[1].split(',')
                if len(valstr) != 3:
                    print('Error: if specifying ages in %s, must be formatted "age = agestart, ageend, n_vals"' % paramfile)
                    exit(-1)
                    
                ages = N.linspace(float(valstr[0]), float(valstr[1]), int(valstr[2]))

            else:
                if (line.strip(' ').startswith("#")) | (len(line.rstrip('\n').strip(' ').strip('\t')) < 1):
                    pass
                else:
                    print("WARNING: unable to parse line in %s:\n   %s" % paramfile, line)

        # end loop through file
    # end with open(paramfile)
    
    grid = N.array(list(product(ages, tau, tq)))
    
    nuv_pred = N.load(col1_file)
    ur_pred  = N.load(col2_file)
    lu = N.append(nuv_pred.reshape(-1,1), ur_pred.reshape(-1,1), axis=1)

except (FileNotFoundError, IOError):
    print("Input file %s not found, trying inputs from STDIN..." % paramfile)
    

    method = raw_input('Do you wish to use a look-up table? (yes/no) :')
    if method == 'yes' or method =='y':
        use_table = True
        prov = raw_input('Do you wish to use the provided u-r and NUV-u look up tables? (yes/no) :')
        if prov == 'yes' or prov =='y':
            print 'gridding...'

            tq = N.linspace(0.003, 13.8, 100)
            tau = N.linspace(0.003, 4, 100)
            ages = N.linspace(10.88861228, 13.67023409, 50)
            grid = N.array(list(product(ages, tau, tq)))
            print 'loading...'
            nuv_pred = N.load('nuv_look_up_ssfr.npy')
            ur_pred = N.load('ur_look_up_ssfr.npy')
            lu = N.append(nuv_pred.reshape(-1,1), ur_pred.reshape(-1,1), axis=1)
        elif prov=='no' or prov=='n':
            col1 = str(raw_input('Location of your NUV-u colour look up table :'))
            col2 = str(raw_input('Location of your u-r colour look up table :'))
            one = N.array(input('Define first axis values (ages) of look up table start, stop, len(axis1); e.g. 10, 13.8, 50 :'))
            ages = N.linspace(float(one[0]), float(one[1]), float(one[2]))
            two = N.array(input('Define second axis values (tau) of look up table start, stop, len(axis1); e.g. 0, 4, 100 : '))
            tau = N.linspace(float(two[0]), float(two[1]), float(two[2]))
            three = N.array(input('Define third axis values (tq) of look up table start, stop, len(axis1); e.g. 0, 13.8, 100 : '))
            tq = N.linspace(float(three[0]), float(three[1]), float(three[2]))
            grid = N.array(list(product(ages, tau, tq)))
            print 'loading...'
            nuv_pred = N.load(col1)
            ur_pred = N.load(col2)
            lu = N.append(nuv_pred.reshape(-1,1), ur_pred.reshape(-1,1), axis=1)
        else:
            sys.exit("You didn't give a valid answer (yes/no). Try running again.")



if use_table:

    def lnlike_one(theta, ur, sigma_ur, nuvu, sigma_nuvu, age):
        """ Function for determining the likelihood of ONE quenching model described by theta = [tq, tau] for all the galaxies in the sample. Simple chi squared likelihood between predicted and observed colours of the galaxies. 
        
        :theta:
        An array of size (1,2) containing the values [tq, tau] in Gyr.
        
        :tq:
        The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time.
        
        :tau:
        The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5.
        
        :ur:
        Observed u-r colour of a galaxy; k-corrected.
        
        :sigma_ur:
        Error on the observed u-r colour of a galaxy
        
        :nuvu:
        Observed nuv-u colour of a galaxy; k-corrected.
        
        :sigma_nuvu:
        Error on the observed nuv-u colour of a galaxy
        
        :age:
        Observed age of a galaxy, often calculated from the redshift i.e. at z=0.1 the age ~ 12.5. Must be in units of Gyr.
        
        RETURNS:
        Array of same shape as :age: containing the likelihood for each galaxy at the given :theta:
            """
        tq, tau = theta
        pred_nuvu, pred_ur = lookup_col_one(theta, age)
        return -0.5*N.log(2*N.pi*sigma_ur**2)-0.5*((ur-pred_ur)**2/sigma_ur**2)-0.5*N.log10(2*N.pi*sigma_nuvu**2)-0.5*((nuvu-pred_nuvu)**2/sigma_nuvu**2)
    
elif method == 'no' or method =='n':
    """We first define the directory in which we will find the BC03 model, extracted from the original files downloaded from the BC03 website into a usable format. Here we implement a solar metallicity model with a Chabrier IMF."""
    model = str(raw_input('Location of the extracted (.ised_ASCII) SPS model to use to predict the u-r and NUV-u colours, e.g. ~/extracted_bc2003_lr_m62_chab_ssp.ised_ASCII :'))
    data = N.loadtxt(model)
    import fluxes 
    def lnlike_one(theta, ur, sigma_ur, nuvu, sigma_nuvu, age):
        """ Function for determining the likelihood of ONE quenching model described by theta = [tq, tau] for all the galaxies in the sample. Simple chi squared likelihood between predicted and observed colours of the galaxies. 
        
            :theta:
            An array of size (1,2) containing the values [tq, tau] in Gyr.
        
            :tq:
            The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time.
            
            :tau:
            The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5.
        
            :ur:
            Observed u-r colour of a galaxy; k-corrected.
        
            :sigma_ur:
            Error on the observed u-r colour of a galaxy
        
            :nuvu:
            Observed nuv-u colour of a galaxy; k-corrected.
        
            :sigma_nuvu:
            Error on the observed nuv-u colour of a galaxy
        
            :age:
            Observed age of a galaxy, often calculated from the redshift i.e. at z=0.1 the age ~ 12.5. Must be in units of Gyr.
        
            RETURNS:
            Array of same shape as :age: containing the likelihood for each galaxy at the given :theta:
            """
        tq, tau = theta
        pred_nuvu, pred_ur = predict_c_one(theta, age)
        return -0.5*N.log(2*N.pi*sigma_ur**2)-0.5*((ur-pred_ur)**2/sigma_ur**2)-0.5*N.log10(2*N.pi*sigma_nuvu**2)-0.5*((nuvu-pred_nuvu)**2/sigma_nuvu**2)

else:
    sys.exit("You didn't give a valid answer (yes/no). Try running again.")

'''


def lnlike_one(theta, ur, sigma_ur, nuvu, sigma_nuvu, age, get_c_one):
    """ Function for determining the likelihood of ONE quenching model described by theta = [tq, tau] for all the galaxies in the sample. Simple chi squared likelihood between predicted and observed colours of the galaxies. 

        :theta:
        An array of size (1,2) containing the values [tq, tau] in Gyr.

        :tq:
        The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time.

        :tau:
        The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5.

        :ur:
        Observed u-r colour of a galaxy; k-corrected.

        :sigma_ur:
        Error on the observed u-r colour of a galaxy

        :nuvu:
        Observed nuv-u colour of a galaxy; k-corrected.

        :sigma_nuvu:
        Error on the observed nuv-u colour of a galaxy

        :age:
        Observed age of a galaxy, often calculated from the redshift i.e. at z=0.1 the age ~ 12.5. Must be in units of Gyr.

        :get_c_one:
        Function to return one colour prediction, maybe predict_c_one if not using lookup table, and lookup_col_one if using a lookup table

        RETURNS:
        Array of same shape as :age: containing the likelihood for each galaxy at the given :theta:
        """
    # why is this next line here, we don't use it at all in this function?
    tq, tau = theta
    pred_nuvu, pred_ur = get_c_one(theta, age)
    return -0.5*N.log(2*N.pi*sigma_ur**2)-0.5*((ur-pred_ur)**2/sigma_ur**2)-0.5*N.log10(2*N.pi*sigma_nuvu**2)-0.5*((nuvu-pred_nuvu)**2/sigma_nuvu**2)



    
n=0



# had to move this to above predict_c_one() and get_colours()

""" Load the magnitude bandpass filters using idl save """
filters = readsav('ugriz.sav')
fuvwave= filters.ugriz.fuvwave[0]
fuvtrans = filters.ugriz.fuvtrans[0]
nuvwave= filters.ugriz.nuvwave[0]
nuvtrans = filters.ugriz.nuvtrans[0]
uwave= filters.ugriz.uwave[0]
utrans = filters.ugriz.utrans[0]
gwave= filters.ugriz.gwave[0]
gtrans = filters.ugriz.gtrans[0]
rwave= filters.ugriz.rwave[0]
rtrans = filters.ugriz.rtrans[0]
iwave= filters.ugriz.iwave[0]
itrans = filters.ugriz.itrans[0]
zwave= filters.ugriz.zwave[0]
ztrans = filters.ugriz.ztrans[0]
vwave= filters.ugriz.vwave[0]
vtrans = filters.ugriz.vtrans[0]
jwave= filters.ugriz.jwave[0]
jtrans = filters.ugriz.jtrans[0]
hwave= filters.ugriz.hwave[0]
htrans = filters.ugriz.htrans[0]
kwave= filters.ugriz.kwave[0]
ktrans = filters.ugriz.ktrans[0]

""" and the HST bandpass filters using numpy save """
""" filter transmission curves are from the SVO filter service,
    e.g. http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=HST/WFC3_IR.F160W 
    
    Roughly speaking, the bandpasses 
    (F435W, F606W, F775W, F814W, F850LP, F105W, F125W, F140W, F160W)
    correspond to
    (B, V, i, I, z, Y, J, JH, H)
    but these aren't exact copies of those filters already read in
    and e.g. it's easy to confuse i and I anyway
    so to be explicit, don't abbreviate to single-letter filters
"""
hst_filters = N.load('HST_filters.npy').flat[0]
f435wave  = hst_filters['HST_ACS_WFC.F435W_77']['wave']
f435trans = hst_filters['HST_ACS_WFC.F435W_77']['throughput']
f606wave  = hst_filters['HST_ACS_WFC.F606W_77']['wave']
f606trans = hst_filters['HST_ACS_WFC.F606W_77']['throughput']
f775wave  = hst_filters['HST_ACS_WFC.F775W_77']['wave']
f775trans = hst_filters['HST_ACS_WFC.F775W_77']['throughput']
f814wave  = hst_filters['HST_ACS_WFC.F814W_77']['wave']
f814trans = hst_filters['HST_ACS_WFC.F814W_77']['throughput']
f850wave  = hst_filters['HST_ACS_WFC.F850LP_77']['wave']
f850trans = hst_filters['HST_ACS_WFC.F850LP_77']['throughput']
f105wave  = hst_filters['HST_WFC3_IR.F105W']['wave']
f105trans = hst_filters['HST_WFC3_IR.F105W']['throughput']
f125wave  = hst_filters['HST_WFC3_IR.F125W']['wave']
f125trans = hst_filters['HST_WFC3_IR.F125W']['throughput']
f140wave  = hst_filters['HST_WFC3_IR.F140W']['wave']
f140trans = hst_filters['HST_WFC3_IR.F140W']['throughput']
f160wave  = hst_filters['HST_WFC3_IR.F160W']['wave']
f160trans = hst_filters['HST_WFC3_IR.F160W']['throughput']


# these are just temporary, don't really use them for science
hst_filters_z08 = N.load('HST_filters_z0.8.npy').flat[0]
f435wave_z08  = hst_filters_z08['HST_ACS_WFC.F435W_77']['wave']
f435trans_z08 = hst_filters_z08['HST_ACS_WFC.F435W_77']['throughput']
f606wave_z08  = hst_filters_z08['HST_ACS_WFC.F606W_77']['wave']
f606trans_z08 = hst_filters_z08['HST_ACS_WFC.F606W_77']['throughput']
f775wave_z08  = hst_filters_z08['HST_ACS_WFC.F775W_77']['wave']
f775trans_z08 = hst_filters_z08['HST_ACS_WFC.F775W_77']['throughput']
f814wave_z08  = hst_filters_z08['HST_ACS_WFC.F814W_77']['wave']
f814trans_z08 = hst_filters_z08['HST_ACS_WFC.F814W_77']['throughput']
f850wave_z08  = hst_filters_z08['HST_ACS_WFC.F850LP_77']['wave']
f850trans_z08 = hst_filters_z08['HST_ACS_WFC.F850LP_77']['throughput']
f105wave_z08  = hst_filters_z08['HST_WFC3_IR.F105W']['wave']
f105trans_z08 = hst_filters_z08['HST_WFC3_IR.F105W']['throughput']
f125wave_z08  = hst_filters_z08['HST_WFC3_IR.F125W']['wave']
f125trans_z08 = hst_filters_z08['HST_WFC3_IR.F125W']['throughput']
f140wave_z08  = hst_filters_z08['HST_WFC3_IR.F140W']['wave']
f140trans_z08 = hst_filters_z08['HST_WFC3_IR.F140W']['throughput']
f160wave_z08  = hst_filters_z08['HST_WFC3_IR.F160W']['wave']
f160trans_z08 = hst_filters_z08['HST_WFC3_IR.F160W']['throughput']


## TEMPORARY FOR HST PROPOSAL PURPOSES
nuvwave  = f435wave_z08
nuvtrans = f435trans_z08
uwave    = f606wave_z08
utrans   = f606trans_z08
rwave    = f850wave_z08
rtrans   = f850trans_z08



def expsfh(tq, tau, time):
    """ This function when given a single combination of [tq, tau] values will calcualte the SFR at all times. First calculate the sSFR at all times as defined by Peng et al. (2010) - then the SFR at the specified time of quenching, tq and set the SFR at this value  at all times before tq. Beyond this time the SFR is an exponentially declining function with timescale tau. 
        
        INPUT:
        :tau:
        The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5. 
        
        :tq: 
        The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time.
        :time:
        An array of time values at which the SFR is calcualted at each step. 
        
        RETURNS:
        :sfr:
        Array of the same dimensions of time containing the sfr at each timestep.
        """
    ssfr = 2.5*(((10**10.27)/1E10)**(-0.1))*(time/3.5)**(-2.2)
    c = time.searchsorted(3.0)
    ssfr[:c] = N.interp(3.0, time, ssfr)
    c_sfr = N.interp(tq, time, ssfr)*(1E10)/(1E9)
    ### definition is for 10^10 M_solar galaxies and per gyr - convert to M_solar/year ###
    a = time.searchsorted(tq)
    sfr = N.ones(len(time))*c_sfr
    sfr[a:] = c_sfr*N.exp(-(time[a:]-tq)/tau)
    return sfr

def expsfh_mass(ur, Mr, age, tq, tau, time):
    """Calculate exponential decline star formation rates at each time step input by matching to the mass of the observed galaxy at the observed time. This is calculated from the mass-to-light ratio that is a function of one color band u-r as in Bladry et al. (2006; see Figure 5) who fit to data from Glazebrrok et al (2004) and Kauffmann et al (2003).
       
       INPUT:

        :ur:
        u-r optical colour, needed to calculate the mass of the observed galaxy

        :Mr:
        Absolute r-band magnitude, needed to calculate the mass of the observed galaxy 

        :age:
        Observed age of a galaxy, often calculated from the redshift i.e. at z=0.1 the age ~ 12.5. Must be in units of Gyr.

        :tq: 
        The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time.

        :tau:
        The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5. 

        :time:
        An array of time values at which the SFR is calcualted at each step. 
        
        RETURNS:
        :sfr:
        Array of the same dimensions of time containing the sfr at each timestep.
        """ 

    t_end = age # time at which to integrate under the exponential curve until to gain the final mass 
    if ur <=2.1:
        log_m_l = -0.95 + 0.56 * ur
    else:
        log_m_l = -0.16 + 0.18 * ur
    m_msun = 10**(((4.62 - Mr)/2.5) + log_m_l)
    print 'Mass [M_solar]', m_msun
    c_sfr = (m_msun/(tq + tau*(1 - N.exp((tq - t_end)/tau)))) / 1E9
    a = time.searchsorted(tq)
    sfr = N.ones(len(time))*c_sfr
    sfr[a:] = c_sfr*N.exp(-(time[a:]-tq)/tau)
    return sfr 






''' BDS modified predict_c_one and get_colors to take more general versions of color prediction requests
    can take any list of filter pairs and return all the colors
    default is still to do it the old way though
'''




def predict_c_one(theta, age, nuv=[nuvwave, nuvtrans], u=[uwave, utrans], r=[rwave, rtrans], filter_pairs=None):
    """ This function predicts the u-r and nuv-u colours of a galaxy with a SFH defined by [tq, tau], according to the BC03 model at a given "age" i.e. observation time. It calculates the colours at all times then interpolates for the observed age - it has to do this in order to work out the cumulative mass across the SFH to determine how much each population of stars contributes to the flux at each time step. 
        
        :theta:
        An array of size (1,2) containing the values [tq, tau] in Gyr.
        
        :tq:
        The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time.
        
        :tau:
        The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5.
        
        :age:
        Observed age of a galaxy, often calculated from the redshift i.e. at z=0.1 the age ~ 12.5. Must be in units of Gyr. 

        If you want colors of the format (a-b), (b-c), you can specify that as:
        nuv=[a_wave, a_trans], u=[b_wave, b_trans], r=[c_wave, c_trans]

        OR for an arbitrary set of colors, instead use:
        filter_pairs = [[a, b], [c, d], ... , [y, z]]
        where each of the filters a, b, c, ... has the format [a_wave, a_trans] etc. as above.
        Note that:
        filter_pairs = [[a, b], [b, c]]
        will have the same result as assigning a, b, and c to nuv, u, and r. However, the nuv=, u=, r= method is a bit faster. 

        
        RETURNS:

        if not specifying filter_pairs:

           :nuv_u_age:
           Array the same shape as :age: with the nuv-u colour values at each given age for the specified :theta: values
        
           :u_r_age:
           Array the same shape as :age: with the u-r colour values at each given age for the specified :theta: values

        otherwise:

           :colours:
           a list of N colours (or colour arrays, if :age: is an array) where N == len(filter_pairs), ordered as filter_pairs is.



        modified 11/7/2018 by BDS to allow user to specify different colours from nuv-u and u-r (those are still default)
           nuv still corresponds to the shortest-wavelength filter, r to the longest

        modified 20/7/2018 by BDS to allow an arbitrary list of pairs of filters
        """
    ti = N.arange(0, 0.01, 0.003)
    t = N.linspace(0,14.0,100)
    t = N.append(ti, t[1:])
    tq, tau = theta
    sfr = expsfh(tq, tau, t)
    ### Work out total flux at each time given the sfh model of tau and tq (calls fluxes function) ###
    total_flux = fluxes.assign_total_flux(data[0,1:], data[1:,0], data[1:,1:], t*1E9, sfr)
    ### Calculate fluxes from the flux at all times then interpolate to get one colour for the age you are observing the galaxy at - if many galaxies are being observed, this also works with an array of ages to give back an array of colours ###


    if filter_pairs is None:
        nuv_u, u_r = get_colours(t*1E9, total_flux, data, nuv=nuv, u=u, r=r)
        nuv_u_age = N.interp(age, t, nuv_u)
        u_r_age = N.interp(age, t, u_r)
        return nuv_u_age, u_r_age
    else:
        colours = get_colours(t*1E9, total_flux, data, filter_pairs)
        colours_age = [N.interp(age, t, thecol) for thecol in colours]
        return colours_age

        '''
        # loop through pairs to get colors
        colours  = []
        
        for thepair in filter_pairs:
            thecolour = get_colour(t*1.0e9, total_flux, data, thepair)
            colours.append(N.interp(age, t, thecolour))

        return colours
        '''
    

def get_colour(time, flux, data, filter_pair):
    """" Calculates the colour of a given sfh fluxes across time given the BC03 models from the magnitudes of the SED.
        
        :time:
        Array of times at which the colours should be calculated. In units of Gyrs. 
        
        :flux:
        SED of fluxes describing the calculated SFH. Returned from the assign_total_flux function in fluxes.py
        
        :data:
        BC03 model values for wavelengths, time steps and fluxes. The wavelengths are needed to calculate the magnitudes. 

        :filter_pair:
        the two filter transmission curves you wish to use to make colours with. 
        format: [[ bluer_wavelengths,  bluer_transmission], 
                 [redder_wavelengths, redder_transmission]]
        
        RETURNS:
        :colour:
        Array the same shape as :time: with the predicted colour

        """
    bluer  = filter_pair[0]
    redder = filter_pair[1]
    i_wave  = 0
    i_trans = 1

    bluer_mag  = fluxes.calculate_AB_mag(time, data[1:,0], flux,  bluer[i_wave],  bluer[i_trans])
    redder_mag = fluxes.calculate_AB_mag(time, data[1:,0], flux, redder[i_wave], redder[i_trans])
    return bluer_mag - redder_mag




def get_colours(time, flux, data, nuv=[nuvwave, nuvtrans], u=[uwave, utrans], r=[rwave, rtrans], filter_pairs=None):
    """" Calculates the colours of a given sfh fluxes across time given the BC03 models from the magnitudes of the SED.
        
        :time:
        Array of times at which the colours should be calculated. In units of Gyrs. 
        
        :flux:
        SED of fluxes describing the calculated SFH. Returned from the assign_total_flux function in fluxes.py
        
        :data:
        BC03 model values for wavelengths, time steps and fluxes. The wavelengths are needed to calculate the magnitudes. 
        
        RETURNS:
        :nuv_u: :u_r:
        Arrays the same shape as :time: with the predicted nuv-u and u-r colours


        modified 11/7/2018 by BDS to allow user to specify different colours from nuv-u and u-r (those are still default)
          note: it's still important that nuv is the bluest wavelength and r is the reddest wavelength
        """

    if filter_pairs is None:
        nuvwave_local  = nuv[0]
        nuvtrans_local = nuv[1]
        uwave_local    = u[0]
        utrans_local   = u[1]
        rwave_local    = r[0]
        rtrans_local   = r[1]
        
        nuvmag = fluxes.calculate_AB_mag(time, data[1:,0], flux, nuvwave_local, nuvtrans_local)
        umag = fluxes.calculate_AB_mag(time, data[1:,0], flux, uwave_local, utrans_local)
        rmag = fluxes.calculate_AB_mag(time, data[1:,0], flux, rwave_local, rtrans_local)
        nuv_u = nuvmag - umag
        u_r = umag - rmag
        return nuv_u, u_r
    else:
        colours  = []
        
        for thepair in filter_pairs:
            bluermag  = fluxes.calculate_AB_mag(time, data[1:,0], flux, thepair[0][0], thepair[0][1])
            reddermag = fluxes.calculate_AB_mag(time, data[1:,0], flux, thepair[1][0], thepair[1][1])
            
            colours.append(bluermag - reddermag)

        return colours
      















def lookup_col_one(theta, age):
    ur_pred = u(theta[0], theta[1])
    nuv_pred = v(theta[0], theta[1])
    return nuv_pred, ur_pred

# Prior likelihood on theta values given the inital w values assumed for the mean and stdev
def lnprior(theta):
    """ Function to calcualted the prior likelihood on theta values given the inital w values assumed for the mean and standard deviation of the tq and tau parameters. Defined ranges are specified - outside these ranges the function returns -N.inf and does not calculate the posterior probability. 
        
        :theta: 
        An array of size (1,4) containing the values [tq, tau] for both smooth and disc galaxies in Gyr.
        
        :tq:
        The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time. Can be either for smooth or disc galaxies.
        
        :tau:
        The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5. Can be either for smooth or disc galaxies.
        
        RETURNS:
        Value of the prior at the specified :theta: value.
        """
    tq, tau = theta
    if 0.003 <= tq <= 13.807108309208775 and 0.003 <= tau <= 4.0:
        return 0.0
    else:
        return -N.inf

# Overall likelihood function combining prior and model
def lnprob(theta, ur, sigma_ur, nuvu, sigma_nuvu, age, get_c_one):
    """Overall posterior function combiningin the prior and calculating the likelihood. Also prints out the progress through the code with the use of n. 
        
        :theta:
        An array of size (1,4) containing the values [tq, tau] for both smooth and disc galaxies in Gyr.
        
        :tq:
        The time at which the onset of quenching begins in Gyr. Allowed ranges from the beginning to the end of known cosmic time. Can be either for smooth or disc galaxies.
        
        :tau:
        The exponential timescale decay rate of the star formation history in Gyr. Allowed range from the rest of the functions is 0 < tau [Gyr] < 5. Can be either for smooth or disc galaxies.
        
        :ur:
        Observed u-r colour of a galaxy; k-corrected. An array of shape (N,1) or (N,).
        
        :sigma_ur:
        Error on the observed u-r colour of a galaxy. An array of shape (N,1) or (N,).
        
        :nuvu:
        Observed nuv-u colour of a galaxy; k-corrected. An array of shape (N,1) or (N,).
        
        :sigma_nuvu:
        Error on the observed nuv-u colour of a galaxy. An array of shape (N,1) or (N,).
        
        :age:
        Observed age of a galaxy, often calculated from the redshift i.e. at z=0.1 the age ~ 12.5. Must be in units of Gyr. An array of shape (N,1) or (N,).
        
        RETURNS:
        Value of the posterior function for the given :theta: value.
        
        """
    global n
    n+=1
    if n %100 == 0:
        print 'step number', n/100
    lp = lnprior(theta)
    if not N.isfinite(lp):
        return -N.inf
    return lp + lnlike_one(theta, ur, sigma_ur, nuvu, sigma_nuvu, age, get_c_one)

def sample(ndim, nwalkers, nsteps, burnin, start, ur, sigma_ur, nuvu, sigma_nuvu, age, id, ra, dec, get_c_one, use_table, thegrid, lu=None):
    """ Function to implement the emcee EnsembleSampler function for the sample of galaxies input. Burn in is run and calcualted fir the length specified before the sampler is reset and then run for the length of steps specified. 
        
        :ndim:
        The number of parameters in the model that emcee must find. In this case it always 2 with tq, tau.
        
        :nwalkers:
        The number of walkers that step around the parameter space. Must be an even integer number larger than ndim. 
        
        :nsteps:
        The number of steps to take in the final run of the MCMC sampler. Integer.
        
        :burnin:
        The number of steps to take in the inital burn-in run of the MCMC sampler. Integer. 
        
        :start:
        The positions in the tq and tau parameter space to start for both disc and smooth parameters. An array of shape (1,4).
        
        
        :ur:
        Observed u-r colour of a galaxy; k-corrected. An array of shape (N,1) or (N,).
        
        :sigma_ur:
        Error on the observed u-r colour of a galaxy. An array of shape (N,1) or (N,).
        
        :nuvu:
        Observed nuv-u colour of a galaxy; k-corrected. An array of shape (N,1) or (N,).
        
        :sigma_nuvu:
        Error on the observed nuv-u colour of a galaxy. An array of shape (N,1) or (N,).
        
        :age:
        Observed age of a galaxy, often calculated from the redshift i.e. at z=0.1 the age ~ 12.5. Must be in units of Gyr. An array of shape (N,1) or (N,).
        
        :id:
        ID number to specify which galaxy this run is for.
        
        :ra:
        right ascension of source, used for identification purposes
        
        :dec:
        declination of source, used for identification purposes
        
        RETURNS:
        :samples:
        Array of shape (nsteps*nwalkers, 4) containing the positions of the walkers at all steps for all 4 parameters.
        :samples_save:
        Location at which the :samples: array was saved to. 
        
        """
    tq, tau, ages = thegrid
    grid = N.array(list(product(ages, tau, tq)))
    
    if use_table:
        global u
        global v
        a = N.searchsorted(ages, age)
        b = N.array([a-1, a])
        print 'interpolating function, bear with...'
        g = grid[N.where(N.logical_or(grid[:,0]==ages[b[0]], grid[:,0]==ages[b[1]]))]
        values = lu[N.where(N.logical_or(grid[:,0]==ages[b[0]], grid[:,0]==ages[b[1]]))]
        f = LinearNDInterpolator(g, values, fill_value=(-N.inf))
        look = f(age, grid[:10000, 1], grid[:10000, 2])
        lunuv = look[:,0].reshape(100,100)
        v = interp2d(tq, tau, lunuv)
        luur = look[:,1].reshape(100,100)
        u = interp2d(tq, tau, luur)
    else:
        pass
    print 'emcee running...'
    p0 = [start + 1e-4*N.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=2, args=(ur, sigma_ur, nuvu, sigma_nuvu, age, get_c_one))
    """ Burn in run here..."""
    pos, prob, state = sampler.run_mcmc(p0, burnin)
    lnp = sampler.flatlnprobability
    N.save('lnprob_burnin_'+str(int(id))+'_'+str(ra)+'_'+str(dec)+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.npy', lnp)
    samples = sampler.chain[:,:,:].reshape((-1,ndim))
    samples_save = 'samples_burn_in_'+str(int(id))+'_'+str(ra)+'_'+str(dec)+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.npy'
    N.save(samples_save, samples)
    sampler.reset()
    print 'Burn in complete...'
    """ Main sampler run here..."""
    sampler.run_mcmc(pos, nsteps)
    lnpr = sampler.flatlnprobability
    N.save('lnprob_run_'+str(int(id))+'_'+str(ra)+'_'+str(dec)+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.npy', lnpr)
    samples = sampler.chain[:,:,:].reshape((-1,ndim))
    samples_save = 'samples_'+str(int(id))+'_'+str(ra)+'_'+str(dec)+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.npy'
    N.save(samples_save, samples)
    print 'Main emcee run completed.'
    return samples, samples_save


#Define function to plot the walker positions as a function of the step
def walker_plot(samples, nwalkers, limit, id):
    """ Plotting function to visualise the steps of the walkers in each parameter dimension for smooth and disc theta values. 
        
        :samples:
        Array of shape (nsteps*nwalkers, 4) produced by the emcee EnsembleSampler in the sample function.
        
        :nwalkers:
        The number of walkers that step around the parameter space used to produce the samples by the sample function. Must be an even integer number larger than ndim.
        
        :limit:
        Integer value less than nsteps to plot the walker steps to. 
        
        :id:
        ID number to specify which galaxy this plot is for.
        
        RETURNS:
        :fig:
        The figure object
        """
    s = samples.reshape(nwalkers, -1, 2)
    s = s[:,:limit, :]
    fig = P.figure(figsize=(8,5))
    ax1 = P.subplot(2,1,1)
    ax2 = P.subplot(2,1,2)
    for n in range(len(s)):
        ax1.plot(s[n,:,0], 'k')
        ax2.plot(s[n,:,1], 'k')
    ax1.tick_params(axis='x', labelbottom='off')
    ax2.set_xlabel(r'step number')
    ax1.set_ylabel(r'$t_{quench}$')
    ax2.set_ylabel(r'$\tau$')
    P.subplots_adjust(hspace=0.1)
    save_fig = 'walkers_steps_'+str(int(id))+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.pdf'
    fig.savefig(save_fig)
    return fig

def corner_plot(s, labels, extents, bf, id):
    """ Plotting function to visualise the gaussian peaks found by the sampler function. 2D contour plots of tq against tau are plotted along with kernelly smooth histograms for each parameter.
        
        :s:
         Array of shape (#, 2) for either the smooth or disc produced by the emcee EnsembleSampler in the sample function of length determined by the number of walkers which resulted at the specified peak.
        
        :labels:
        List of x and y axes labels i.e. disc or smooth parameters
        
        :extents:
        Range over which to plot the samples, list shape [[xmin, xmax], [ymin, ymax]]
        
        :bf:
        Best fit values for the distribution peaks in both tq and tau found from mapping the samples. List shape [(tq, poserrtq, negerrtq), (tau, poserrtau, negerrtau)]
        
        :id:
        ID number to specify which galaxy this plot is for. 
        
        RETURNS:
        :fig:
        The figure object
        """
    x, y = s[:,0], s[:,1]
    fig = P.figure(figsize=(6.25,6.25))
    ax2 = P.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
    ax2.set_xlabel(labels[0])
    ax2.set_ylabel(labels[1])
    # newer versions of matplotlib will give a "contours must be increasing" error,
    # which is fixed in the next-generation version of triangle (i.e. corner)
    #triangle.hist2d(x, y, ax=ax2, bins=100, extent=extents, plot_contours=True)
    corner.hist2d(x, y, ax=ax2, bins=100, extent=extents, plot_contours=True)
    ax2.axvline(x=bf[0][0], linewidth=1)
    ax2.axhline(y=bf[1][0], linewidth=1)
    [l.set_rotation(45) for l in ax2.get_xticklabels()]
    [j.set_rotation(45) for j in ax2.get_yticklabels()]
    ax2.tick_params(axis='x', labeltop='off')
    ax1 = P.subplot2grid((3,3), (0,0),colspan=2)
    den = kde.gaussian_kde(x[N.logical_and(x>=extents[0][0], x<=extents[0][1])])
    pos = N.linspace(extents[0][0], extents[0][1], 750)
    ax1.plot(pos, den(pos), 'k-', linewidth=1)
    ax1.axvline(x=bf[0][0], linewidth=1)
    ax1.axvline(x=bf[0][0]+bf[0][1], c='b', linestyle='--')
    ax1.axvline(x=bf[0][0]-bf[0][2], c='b', linestyle='--')
    ax1.set_xlim(extents[0][0], extents[0][1])
    ax12 = ax1.twiny()
    ax12.set_xlim(extents[0][0], extents[0][1])
    ax12.set_xticks(N.array([1.87, 3.40, 6.03, 8.77, 10.9, 12.5]))
    ax12.set_xticklabels(N.array([3.5, 2.0 , 1.0, 0.5, 0.25, 0.1]))
    [l.set_rotation(45) for l in ax12.get_xticklabels()]
    ax12.tick_params(axis='x', labelbottom='off')
    ax12.set_xlabel(r'$z$')
    ax1.tick_params(axis='x', labelbottom='off', labeltop='off')
    ax1.tick_params(axis='y', labelleft='off')
    ax3 = P.subplot2grid((3,3), (1,2), rowspan=2)
    ax3.tick_params(axis='x', labelbottom='off')
    ax3.tick_params(axis='y', labelleft='off')
    den = kde.gaussian_kde(y[N.logical_and(y>=extents[1][0], y<=extents[1][1])])
    pos = N.linspace(extents[1][0], extents[1][1], 750)
    ax3.plot(den(pos), pos, 'k-', linewidth=1)
    ax3.axhline(y=bf[1][0], linewidth=1)
    ax3.axhline(y=bf[1][0]+bf[1][1], c='b', linestyle='--')
    ax3.axhline(y=bf[1][0]-bf[1][2], c='b', linestyle='--')
    ax3.set_ylim(extents[1][0], extents[1][1])
    if os.path.exists(str(int(id))+'.jpeg') == True:
        ax4 = P.subplot2grid((3,3), (0,2), rowspan=1, colspan=1)
        img = mpimg.imread(str(int(id))+'.jpeg')
        ax4.imshow(img)
        ax4.tick_params(axis='x', labelbottom='off', labeltop='off')
        ax4.tick_params(axis='y', labelleft='off', labelright='off')
    P.tight_layout()
    P.subplots_adjust(wspace=0.0)
    P.subplots_adjust(hspace=0.0)
    return fig



