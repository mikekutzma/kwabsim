import cavitySolverMain as csm


LAM = 1064*10**(-7) #cm
M2 = 1
n=1

cfile = '/home/tarkan/kwabsim/datatest.dat'

cavity = csm.getCavity(csm.getStringCav(cfile))

print(csm.check_stab(csm.get_RTM(cavity)))

actions = csm.check_actions(cavity)
#print(actions)

if ((actions[0] + actions[1]) > 0):
    cavities = csm.split_cavity(cavity)
    csm.plot_xy(cavities,LAM,100)
    
else:
    [Z,Y] = csm.cavityRun(cavity,LAM,100)
    csm.Plot(Z[0],Y[0])




