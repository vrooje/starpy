from posterior import *
from astropy.cosmology import FlatLambdaCDM
import numpy as N
import sys, os, time
from scipy.stats import kde
from scipy import interpolate
from scipy.integrate import simps
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp2d
from itertools import product


# Use sys to assign arguments for the galaxy data from the command line
try:
    # this is the default, running starpy once on one source
    u_r, err_u_r, nuv_u, err_nuv_u, z, dr8, ra, dec = sys.argv[1:]
    rows = [[u_r, err_u_r, nuv_u, err_nuv_u, z, dr8, ra, dec]]
    many_sources = False
except:
    # if the above doesn't work assume the first input points to a file with a list of colors for many sources
    # the inputs should have the same structure as the above, with spaces between parameters
    objlist_file = sys.argv[1]
    many_sources = True
    #lists = [[] for i in range(8)]
    #u_r, err_u_r, nuv_u, err_nuv_u, z,  dr8,  ra,  dec = lists
    rows = []

    with open(objlist_file) as fobj:
        for i_l, line in enumerate(fobj):
            arg = line.strip('\n').strip(' ').split(' ')
            if not(len(arg) == 8):
                print("Something wrong at line %d in file %s, got %d values instead of 8" % (i_l, objlist_file, len(arg)))
                exit(-1)

            rows.append(arg)

    print(" Read %d objects from file %s." % (len(rows), objlist_file))

    
# Use astropy to calculate the age from the redshift in the data 
cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
#age = N.array(cosmo.age(float(z)))

'''

 26/07/2018 - edited by BDS to move input parameters into this file so they don't have to be 
  read in separately in posterior.py and fluxes.py

'''

# defaults
tq = N.linspace(0.003, 13.8, 100)
tau = N.linspace(0.003, 4, 100)
ages = N.linspace(10.88861228, 13.67023409, 50)
col1_file  = 'nuv_look_up_ssfr.npy'
col2_file  = 'ur_look_up_ssfr.npy'
model      = 'models/Padova1994/chabrier/ASCII/extracted_bc2003_lr_m62_chab_ssp.ised_ASCII'
use_table  = False
paramfile  = ''
write_params = False

plotdir = "./"
savedir = "./"

paramfile = 'posterior_params.in'
try:
    with open(paramfile) as f:
        for line in f:
            arg = line.rstrip('\n').strip(' ').split('=')
            #print(arg)
            #print("------\n")

            if (arg[0].lower().strip() in ['lookup', 'lookups', 'lookuptable', 'lookuptables', 'lu', 'lut']):
                use_table = True
                tables = arg[1].split(',')

                if tables[0].strip() in ['default']:
                    lu_default = True
                else:
                    if len(tables) < 2:
                        tables = arg[1].split(' ')
                        if len(tables) < 2:
                            print('Error: if not "default", 2 lookup tables needed, filenames separated by "," or " "')
                            exit(-1)
                    col1_file  = tables[0].strip()
                    col2_file  = tables[1].strip()
                    lu_default = False

            elif (arg[0].lower().strip() in ['tq', 't_q', 'tquench', 't_quench', 'quench_time', 'quenching_time']):
                valstr = arg[1].split(',')
                if len(valstr) != 3:
                    print('Error: if specifying quenching times tq in %s, must be formatted "tq = tstart, tend, n_vals"' % paramfile)
                    exit(-1)
                    
                tq = N.linspace(float(valstr[0]), float(valstr[1]), int(valstr[2]))

            elif (arg[0].lower().strip() in ['tau', 'exptau', 'quenchingrate', 'quenching_rate']):
                valstr = arg[1].split(',')
                if len(valstr) != 3:
                    print('Error: if specifying quenching rates tau in %s, must be formatted "tau = taustart, tauend, n_vals"' % paramfile)
                    exit(-1)
                    
                tau = N.linspace(float(valstr[0]), float(valstr[1]), int(valstr[2]))

            elif (arg[0].lower().strip() in ['ages', 'age', 't_obs', 'age_obs']):
                valstr = arg[1].split(',')
                if len(valstr) != 3:
                    print('Error: if specifying ages in %s, must be formatted "age = agestart, ageend, n_vals"' % paramfile)
                    exit(-1)
                    
                ages = N.linspace(float(valstr[0]), float(valstr[1]), int(valstr[2]))

            elif (arg[0].lower().strip() in ['model', 'models']):
                if (arg[1].strip() in ['default']):
                    # we've already defined the default above
                    pass
                else:
                    model = arg[1].strip()

            elif (arg[0].lower().strip() in ['save', 'savedir', 'save_dir']):
                if (arg[1].strip() in ['default']):
                    # we've already defined the default above
                    pass
                else:
                    savedir = arg[1].strip()
                    if not savedir.endswith("/"):
                        savedir += "/"

            elif (arg[0].lower().strip() in ['plot', 'plotdir', 'plot_dir']):
                if (arg[1].strip() in ['default']):
                    # we've already defined the default above
                    pass
                else:
                    plotdir = arg[1].strip()
                    if not plotdir.endswith("/"):
                        plotdir += "/"

            elif (arg[0].lower().strip() in ['params_out', 'paramfile', 'paramfile_out']):
                paramfile = arg[1].strip()
                write_params = True
                fparams = open(paramfile, "w")


            else:
                if (line.strip(' ').startswith("#")) | (len(line.rstrip('\n').strip(' ').strip('\t')) < 1):
                    pass
                else:
                    print("WARNING: unable to parse line in %s:\n%s" % (paramfile, line))
                    #print(arg)

        # end loop through file
    # end with open(paramfile)
    
    grid = N.array(list(product(ages, tau, tq)))
    
    nuv_pred = N.load(col1_file)
    ur_pred  = N.load(col2_file)
    lu = N.append(nuv_pred.reshape(-1,1), ur_pred.reshape(-1,1), axis=1)

except IOError as e:
    print("Oops!\n\n")
    print(e)
    print("\n")
    print("Input file %s not found or there was an error reading in a file within it, trying inputs from STDIN..." % paramfile)
    
    model = str(raw_input('Tell me the location of the extracted (.ised_ASCII) SPS model to use to predict the u-r and NUV-u colours, e.g. ~/extracted_bc2003_lr_m62_chab_ssp.ised_ASCII :'))
    
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
            col1_file = str(raw_input('Location of your NUV-u colour look up table :'))
            col2_file = str(raw_input('Location of your u-r colour look up table :'))
            one = N.array(input('Define first axis values (ages) of look up table start, stop, len(axis1); e.g. 10, 13.8, 50 :'))
            ages = N.linspace(float(one[0]), float(one[1]), float(one[2]))
            two = N.array(input('Define second axis values (tau) of look up table start, stop, len(axis1); e.g. 0, 4, 100 : '))
            tau = N.linspace(float(two[0]), float(two[1]), float(two[2]))
            three = N.array(input('Define third axis values (tq) of look up table start, stop, len(axis1); e.g. 0, 13.8, 100 : '))
            tq = N.linspace(float(three[0]), float(three[1]), float(three[2]))
            grid = N.array(list(product(ages, tau, tq)))
            print 'loading...'
            nuv_pred = N.load(col1_file)
            ur_pred = N.load(col2_file)
            lu = N.append(nuv_pred.reshape(-1,1), ur_pred.reshape(-1,1), axis=1)
        else:
            sys.exit("You didn't give a valid answer (yes/no). Try running again.")



print("Parameters and models used:")
print("Model file: %s" % model)
if use_table:
    print("Lookup files used: \n   bluer colour: %s\n   redder colour: %s" % (col1_file, col2_file))
else:
    print("Not using lookup table, predicting colours from model directly (this is VERY SLOW).")
    print(".... seriously, if you are running this a lot you should make a lookup table first!")

print("Saving plots to %s" % plotdir)
print("Saving .npy files to %s" % savedir)
    
print("Grid used:\n")
print("   quenching time tq  varies from %.4f to %.4f Gyr, in %d steps" % (min(tq), max(tq), len(tq)))
print("   quenching rate tau varies from %.4f to %.4f,     in %d steps" % (min(tau), max(tau), len(tau)))
print("   pop ages covered   varies from %.4f to %.4f Gyr, in %d steps" % (min(ages), max(ages), len(ages)))

if many_sources:
    print("Beginning computations for %s sources..." % len(rows))
    
# this bit was previously in fluxes.py

data = N.loadtxt(model)
model_ages = data[0,1:]
model_lambda = data[1:,0]
model_fluxes = data[1:,1:]
time_flux = N.arange(0, 0.01, 0.003)
t_flux = N.linspace(0,14.0,100)
time_steps_flux = N.append(time_flux, t_flux[1:])*1E9
#First mask the ages of the very young stars hidden in birth clouds
mask = model_ages[model_ages<4E6]
model_fluxes[:,0:len(mask)] = 0.0
# Calculate the fluxes at the ages specified by the time steps rather than in the models using numpy/scipy array manipulations rather than a for loop
f = interpolate.interp2d(model_ages, model_lambda, model_fluxes)
interp_fluxes_sim = f(time_steps_flux, model_lambda)









# Define parameters needed for emcee 
nwalkers = 100 # number of monte carlo chains
nsteps= 800 # number of steps in the monte carlo chain
start = [7.5, 1.5] # starting place of all the chains
burnin = 400 # number of steps in the burn in phase of the monte carlo chain

#The rest calls the emcee module which is initialised in the sample function of the posterior file.
if use_table:
    the_c_function = lookup_col_one
    lookup = lu #N.append(nuv_pred.reshape(-1,1), ur_pred.reshape(-1,1), axis=1)
else:
    the_c_function = predict_c_one
    lookup = None


for i_row in range(len(rows)):
    u_r, err_u_r, nuv_u, err_nuv_u, z, dr8, ra, dec = rows[i_row]

    if many_sources:
        print("======= Beginning run %d =======")

    age = N.array(cosmo.age(float(z)))

    print("Input colors are:\n   bluer = %s +/- %s\b   redder = %s +/- %s" % (nuv_u, err_nuv_u, u_r, err_u_r))
    print("for source %s at redshift z = %s, i.e. age = %.2f Gyr,\n   and (RA, Dec) = (%s, %s)" % (dr8, z, age, ra, dec))

    it_worked = False
    
    try:
        samples, samples_save = sample(2, nwalkers, nsteps, burnin, start, float(u_r), float(err_u_r), float(nuv_u), float(err_nuv_u), age, dr8, ra, dec, the_c_function, use_table, (tq, tau, ages), lu=lookup, savedir=savedir)
        it_worked = True
        
    except Exception as e:
        print("******************* WHOOPS -- SOMETHING WENT WRONG FOR ID %s *******************")
        print(e)
        print("\n                    We shall skip this one.... onwards!\n")
        print("********************************************************************************\n\n")
        
        
    if it_worked:
        tq_mcmc, tau_mcmc,  = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0], v[4]-v[1], v[1]-v[3]), zip(*N.percentile(samples, [16,50,84,2.5,97.5],axis=0)))
        print 'Best fit [t, tau] values found by starpy for input parameters are : [', tq_mcmc[0], tau_mcmc[0], ']'
        fig = corner_plot(samples, labels = [r'$ t_{quench}$', r'$ \tau$'], extents=[[N.min(samples[:,0]), N.max(samples[:,0])],[N.min(samples[:,1]),N.max(samples[:,1])]], bf=[tq_mcmc, tau_mcmc], id=dr8)
        fig.savefig(plotdir+'starpy_output_'+str(dr8)+'_'+str(ra)+'_'+str(dec)+'.pdf')

        if write_params:
            fparams.write("%f %f %f %f %f %f %f %f %f %f" % (tq_mcmc[0], tau_mcmc[0], tq_mcmc[1], tau_mcmc[1], tq_mcmc[2], tau_mcmc[2], tq_mcmc[3], tau_mcmc[3], tq_mcmc[4], tau_mcmc[4]))



        
if write_params:
    fparams.close()
