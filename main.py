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
        rt = np.dot(rt,m)
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
    #print(w2*math.pi/(1064*10**(-7)))
    q = 1/(rover1 - (lam/(math.pi*w2))*1j)
    q = -q.real + q.imag*1j
    return q

def prop_q(q,rtm):
    a,b,c,d =[rtm[0][0],rtm[0][1],rtm[1][0],rtm[1][1]]
    qo = (q*a + b)/(q*c + d)
    return qo

def get_hyperbs(q,x): 
    #Input list of domain and initial q, return list of y's

    z0 = q.real
    zr = q.imag
    ys = np.sqrt((4*LAM/math.pi)*(zr + ((x-z0)**2)/zr))
    
    return ys

def getDs(cavity):
    #Compile list of where each free space matrix begins

    dtotal = 0 #running total of distance
    Ds = [] #running list of initial D locations

    for opt in cavity:
    #scanning each indivual matrix of cavity
        
        if opt[0]=='D':
            #print('opt: ',opt)
            Ds.append([dtotal,opt[1]]) #(D start point, length of D)
            dtotal = dtotal + opt[1]
    
    return Ds

def getXs(Ds,res):
    #Return a list of arrays for local hyperbolas (x=0 for each 
    #hyperbola will be the beginning of each free space)

    Xs = []
    dtotal = Ds[-1][0]+ Ds[-1][1]
    #print('Ds[-1] :',Ds[-1])
    #print('dtotal = ', dtotal)

    for D in Ds:
        #scan through each free space
        Xs.append(np.linspace(0,D[1],num=round(res/dtotal*D[1])))

    return Xs

def getY(q,X):
    #return list of Ys for an individual X block

    z0 = -q.real
    zr = q.imag
    
    #print('z0,zr = ',z0,zr)

    Ys = np.sqrt((4*LAM/math.pi)*(zr + ((X-z0)**2)/zr))

    return Ys

def getQs(q,cavity):
    #create list of Qs for the beginning of each free space
    Qs = [q]
    firstD = False #This will tell us if we've passed the first 'D'
    qrun = q #This is our running Q paramter
    Dstart = 1 #Tells to start cavity search with first element 
# or not

    if(cavity[0][0]=='D'):
        Dstart = 0

    for opt in cavity[Dstart:]:
        #we will keep propogating q through each cavity part. When
        # we reach a 'D' matrix, we'll stop a record our progress
        
        #We must ignore this for the first 'D' though because we 
        # already have the q

        #print('opt: ',opt)
        if(opt[0]=='D'):
            #Our current optic is free space


            
            if(firstD==False):
                firstD = True
                qrun = prop_q(qrun,optics[opt[0]](opt[1]))
                #print('Converted firstD check: ',firstD)
                continue
                #We acknowledge that we already have Q for the first
                #free space and skip it
                
            Qs.append(qrun)
            qrun = prop_q(qrun,optics[opt[0]](opt[1]))
            #We record qrun into Qs, and then propogate it

        else:
            #We aren't on a 'D'. Who gives a shit then, propogate
            qrun = prop_q(qrun,optics[opt[0]](opt[1]))

    return Qs

def combine_graphs(Xs,Ds,Qs):
    #We will now connect blocks of graphs
    newX = [] #This will hold all D block domains
    newY = [] #This will hold all the Y values 
    
    for ind in range(0,len(Xs)):
        #ind denotes index of block number

        tempY =  getY(Qs[ind],Xs[ind])
        for y in tempY:
            newY = newY + [y]
        #We place ind'th Y block into Y
        
        for x in Xs[ind]:
            #x is now floating type
            newX = newX + [x + Ds[ind][0]]
           
    return [newX,newY]

def getRs(Qs):
    #Convert a list of Qs to get the radius of curvature
    
    Rs = []

    for q in Qs:
        Rs.append((((1/q).real))**(-1))

    return Rs


optics = {
        'D':m_freeSpace,
        'L':m_lens,
        'M':m_mirror
         }
LAM = 1064*10**(-7)
#LAM = 1000*10**(-5)

#cavity = [['eM',0],['D',60],['L',200],['D',50],['eM',0]]
#cavity = [['eM',1],['D',.4],['D',.5],['eM',1]]
#cavity = [['M',0],['D',40],['M',5],['D',5.2],['M',5],['D',40],['M',0]]
cavity = [['M',200],['D',34],['M',0],['D',3],['L',18],['D',3],['M',-40],['D',61],['L',30],['D',105],['L',30],['D',41],['M',0]]

rtm = get_RTM(unfold_Cav(cavity))


if(check_stab(rtm)==True):
    print('This cavity is stable')

    res = 100
    
    q = get_q(rtm,LAM)
    Ds = getDs(cavity)
    Xs = getXs(Ds,res)
    Qs = getQs(q,cavity)

    Rs = getRs(Qs)
    print('Rs: \n',Rs)

    X,Y = combine_graphs(Xs,Ds,Qs)

    plt.plot(X,Y)
    plt.show()

else:
    print('Cavity not stable')
