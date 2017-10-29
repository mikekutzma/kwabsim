import cavitySolverMain as csm
import numpy as np

LAM = 1064*10**(-7) #cm
M2 = 1
n=1

cfile = '/home/tarkan/kwabsim/datatest.dat'

cavity = csm.getCavity(csm.getStringCav(cfile))

qx1 = 38.78+16.55j
qy1 = 33.04+16.642j

Qx =-38.598+8.8216j
Qy = -33.564 + 10.495j

minx = 99999
miny = 99999
start = 1
end = 30

for fl in np.linspace(start,end,(end-start+1)):
    print('fl = ',fl)
    cavity = [['L',fl]]

    qx2 = csm.prop_q(qx1,csm.get_RTM(cavity))
    qy2 = csm.prop_q(qy1,csm.get_RTM(cavity))
    
    deltax = qx2 - Qx
    deltay = qy2 - Qy

    if(abs(deltax) < minx):
        minx = abs(deltax)
        minx_fl = fl
    if(abs(deltay)<miny):
        miny = abs(deltay)
        miny_fl = fl
    
#end loop

print('minx_fl = ',minx_fl,'\t deltax = ',deltax,)
print('miny_fl = ',miny_fl,'\t deltay = ',deltay,)



