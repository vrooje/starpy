from posterior import *
from itertools import product
import numpy as N

""" Function to generate a look up table of u-r and NUV-u colours using the predict_c_one function in StarfPy. Defaults are for a 50 x 100 x 100 look up table in age, tau and t. 
    """
# defaults
# tq = N.linspace(0.003, 13.8, 100)
tau = N.linspace(0.003, 4, 100)
# age = N.linspace(10.88861228, 13.67023409, 50)  #this is approx. 0.25 > z > 0

# don't have quenching time go all the way to z=0 when we aren't going to go z < 0.2
tq = N.linspace(0.003, 11.5, 100)
# higher-redshift (times) for HST images - approx 4 > z > 0.2 
# in equal dt as the above. which means a LOT more bins
age = N.linspace(1.571, 11.258, 175)

"""
# temporary, just testing things out
tq = N.linspace(0.003, 11.0, 10)
tau = N.linspace(0.003, 4, 10)
age = N.linspace(1.571, 11.258, 10)
"""

print 'making product list of inputs...'
grid = N.array(list(product(age, tau, tq)))

ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))

# savename1 = str(raw_input('What should I save the first lookup table as? e.g. "~/col1_look_up.npy" : '))
# savename2 = str(raw_input('What should I save the second lookup table as? e.g. "~/col2_look_up.npy" : '))

savenames = {'f435-f606': 'f435mf606_look_up.npy',
             'f606-f775': 'f606mf775_look_up.npy',
             'f775-f850': 'f775mf850_look_up.npy',
             'f850-f105': 'f850mf105_look_up.npy',
             'f850-f125': 'f850mf125_look_up.npy',
             'f606-f814': 'f606mf814_look_up.npy',
             'f814-f105': 'f814mf105_look_up.npy',
             'f814-f125': 'f814mf125_look_up.npy',
             'f105-f125': 'f105mf125_look_up.npy',
             'f125-f160': 'f125mf160_look_up.npy',
             'f606-f850': 'f606mf850_look_up.npy',
             'f435-f606 z0.8': 'f435mf606_z08_look_up.npy',
             'f606-f850 z0.8': 'f606mf850_z08_look_up.npy'}

'calculating colours...'

'''
  writing to the savefile takes ~50-100 ms each time for a file ~10 MB size, so as the iterations
  increase, the time spent writing increases considerably as well.
  For a grid of 50 x 50 x 175 and writing time of 0.1 ms, writing each iteration adds up to 48 hours
  of CPU time. So let's not write *every* iteration, but instead just every so often. That way if
  we crash we won't have that much to re-do and can pick up again at close to where we left off.
  Note: adding an if statement to check if we're at the right n value costs slightly extra, but the
  time saved more than makes up for it.
'''
# scaled integer values between 1 and 500
i_savefile = max(min(int(len(grid)/3500), 500), 1)


colstr =  'F435-F606, F606-F850 at z = 0.8'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    # nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0])
    # N.save(savename1, nuv)
    # N.save(savename2, ur)
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f435wave_z08, f435trans_z08], u=[f606wave_z08, f606trans_z08], r=[f850wave_z08, f850trans_z08])
    if n % i_savefile == 0:
        N.save(savenames['f435-f606 z0.8'], nuv)
        N.save(savenames['f606-f850 z0.8'], ur)

# write the final array (needed if we're not writing to the savefile with each iteration)
N.save(savenames['f435-f606 z0.8'], nuv)
N.save(savenames['f606-f850 z0.8'], ur)
print("Done %s\n" % colstr)
    
'''
colstr =  'F435-F606, F606-F775'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    # nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0])
    # N.save(savename1, nuv)
    # N.save(savename2, ur)
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f435wave, f435trans], u=[f606wave, f606trans], r=[f775wave, f775trans])
    if n % i_savefile == 0:
        N.save(savenames['f435-f606'], nuv)
        N.save(savenames['f606-f775'], ur)

# write the final array (needed if we're not writing to the savefile with each iteration)
N.save(savenames['f435-f606'], nuv)
N.save(savenames['f606-f775'], ur)
print("Done %s\n" % colstr)
    

ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F775-F850, F850-F105' 
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f775wave, f775trans], u=[f850wave, f850trans], r=[f105wave, f105trans])
    if n % i_savefile == 0:
        N.save(savenames['f775-f850'], nuv)
        N.save(savenames['f850-f105'], ur)

# write the final array (needed if we're not writing to the savefile with each iteration)
N.save(savenames['f775-f850'], nuv)
N.save(savenames['f850-f105'], ur)
print("Done %s\n" % colstr)
    

ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F850-F125, F125-F160'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f850wave, f850trans], u=[f125wave, f125trans], r=[f160wave, f160trans])
    if n % i_savefile == 0:
        N.save(savenames['f850-f125'], nuv)
        N.save(savenames['f125-f160'], ur)

# write the final array (needed if we're not writing to the savefile with each iteration)
N.save(savenames['f850-f125'], nuv)
N.save(savenames['f125-f160'], ur)
print("Done %s\n" % colstr)
    

ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F606-F814, F814-F125'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f606wave, f606trans], u=[f814wave, f814trans], r=[f125wave, f125trans])
    if n % i_savefile == 0:
        N.save(savenames['f606-f814'], nuv)
        N.save(savenames['f814-f125'], ur)

# write the final array (needed if we're not writing to the savefile with each iteration)
N.save(savenames['f606-f814'], nuv)
N.save(savenames['f814-f125'], ur)
print("Done %s\n" % colstr)
    


ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F814-F105, F105-F125'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f814wave, f814trans], u=[f105wave, f105trans], r=[f125wave, f125trans])
    if n % i_savefile == 0:
        N.save(savenames['f814-f105'], nuv)
        N.save(savenames['f105-f125'], ur)

# write the final array (needed if we're not writing to the savefile with each iteration)
N.save(savenames['f814-f105'], nuv)
N.save(savenames['f105-f125'], ur)
print("Done %s\n" % colstr)
    
'''

#bye 
