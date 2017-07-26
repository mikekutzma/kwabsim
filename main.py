import numpy as np
import math
import matplotlib.pyplot as plt

#Ray Trace Matricies
#================================================

#Free Space
def m_freeSpace(dist):
    #|1 d|
    #|0 1|
    return np.array([[1,dist],[0,1]])

#Lens
def m_lens(f):
    #|1    0|
    #|-1/f 1|
    if(f==0):
        return np.array([[1,0],[0,1]])
    return np.array([[1,0],[-1/f,1]])

#Mirror
def m_mirror(roc):
    #|1     0|
    #|-2/f  1|
    return m_lens(roc/2)

#===================================================

def get_RTM(cavity):
    rt = np.array([[1,0],[0,1]])
    for optic in list(reversed(cavity)):
        m = optics[optic[0]](optic[1])
        rt = np.matmul(rt,m)
    return rt

def unfold_Cav(cavity):
    return cavity[0:-1]+list(reversed(cavity))[0:-1]

def check_stab(rtm):
    x = (rtm[0][0]+rtm[1][1]+2)/4.0
    if((x>0) & (x<1)):
        return True
    return False

def get_q(rtm,lam):
    a,b,c,d =[rtm[0][0],rtm[0][1],rtm[1][0],rtm[1][1]]
    rover1 = (d-a)/(2*b)
    w2 = ((lam/math.pi)*abs(b)/(math.sqrt(1 - ((a+d)/2)**2)))
    print(w2*math.pi/(1064*10**(-7)))
    q = 1/(rover1 - (lam/(math.pi*w2))*1j)
    q = -q.real + q.imag
    return q

def prop_q(q,rtm):
    a,b,c,d =[rtm[0][0],rtm[0][1],rtm[1][0],rtm[1][1]]
    qo = (q*a + b)/(q*c + d)
    return qo

def get_hyperbs(q,cavity): 
    inds = []
    for i in range(len(cavity)):
        if(cavity[i][0]=='D'):
            inds.append(i)
    ds = [cavity[i] for i in inds]
    ops = []
    for m in range(len(inds)-1):
        i,j = [inds[m],inds[m+1]]
        small_cav = cavity[i+1:j]
        ops.append(small_cav)
    qs = [q]
    for m in range(len(ds)-1):
        small_cav = [ds[m]]+ops[m]
        qs.append(prop_q(qs[m],get_RTM(small_cav)))
    return [qs,[d[1] for d in ds]]

optics = {
        'D':m_freeSpace,
        'L':m_lens,
        'eM':m_mirror,
        'iM':m_mirror
        }
LAM = 1064*10**(-7)
cavity = [['eM',1],['D',0.99999999],['eM',1]]
rtm = get_RTM(unfold_Cav(cavity))
print(check_stab(rtm))
q = get_q(rtm,LAM)
z0 = q.real
zr = q.imag

x = np.array(range(0,101))/100
w2 = (4*LAM/math.pi)*(zr + ((x-z0)**2)/zr)
w = np.sqrt(w2)
plt.plot(x,w)
plt.show()
