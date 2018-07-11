from posterior import *
from itertools import product
import numpy as N

""" Function to generate a look up table of u-r and NUV-u colours using the predict_c_one function in StarfPy. Defaults are for a 50 x 100 x 100 look up table in age, tau and t. 
    """
# defaults
tq = N.linspace(0.003, 13.8, 100)
tau = N.linspace(0.003, 4, 100)
# age = N.linspace(10.88861228, 13.67023409, 50)  #this is approx. 0.25 > z > 0

# higher-redshift (times) for HST images - approx 4 > z > 0.2
# in equal dt as the above. which means a LOT more bins
age = N.linspace(1.571, 11.258, 175)

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
             'f125-f160': 'f125mf160_look_up.npy'}

'calculating colours...'


colstr =  'F435-F606, F606-F775'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    # nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0])
    # N.save(savename1, nuv)
    # N.save(savename2, ur)
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f435wave, f435trans], u=[f606wave, f606trans], r=[f775wave, f775trans])
    N.save(savenames['f435-f606'], nuv)
    N.save(savenames['f606-f775'], ur)


ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F775-F850, F850-F105' 
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f775wave, f775trans], u=[f850wave, f850trans], r=[f105wave, f105trans])
    N.save(savename['f775-f850'], nuv)
    N.save(savename['f850-f105'], ur)


ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F850-F125, F125-F160'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f850wave, f850trans], u=[f125wave, f125trans], r=[f160wave, f160trans])
    N.save(savename['f850-f125'], nuv)
    N.save(savename['f125-f160'], ur)


ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F606-F814, F814-F125'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f606wave, f606trans], u=[f814wave, f814trans], r=[f125wave, f125trans])
    N.save(savename['f606-f814'], nuv)
    N.save(savename['f814-f125'], ur)



ur = N.zeros(len(grid))
nuv = N.zeros(len(grid))
colstr =  'F814-F105, F105-F125'
for n in range(len(grid)):
    if n%10000 == 0:
        print('%d percent complete for %s' % ((float(n)/len(grid))*100, colstr))
    nuv[n], ur[n] = predict_c_one([grid[n,2], grid[n,1]], grid[n,0], nuv=[f814wave, f814trans], u=[f105wave, f105trans], r=[f125wave, f125trans])
    N.save(savename['f814-f105'], nuv)
    N.save(savename['f105-f125'], ur)


