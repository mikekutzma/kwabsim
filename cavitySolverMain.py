import numpy as np
import math
import matplotlib.pyplot as plt


LAM = 1064*10**(-7)
M2 = 1
n=1
Display_Parameters = True

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

optics = {
    'D':m_freeSpace,
    'L':m_lens,
    'M':m_mirror,
    'Cx':m_lens,
    'Cy':m_lens
}

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

def getZs(Ds,res):
    #Return a list of arrays for local hyperbolas (x=0 for each 
    #hyperbola will be the beginning of each free space)

    Zs = []
    dtotal = Ds[-1][0]+ Ds[-1][1]
    #print('Ds[-1] :',Ds[-1])
    #print('dtotal = ', dtotal)

    for D in Ds:
        #scan through each free space
        Zs.append(np.linspace(0,D[1],num=round(res/dtotal*D[1])))

    return Zs

def getY(q,Z,LAM):
    #return list of Ys for an individual Z block

    z0 = -q.real
    zr = q.imag
    
    #print('z0,zr = ',z0,zr)

    Ys = np.sqrt((4*LAM/math.pi)*(zr + ((Z-z0)**2)/zr))

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

def combine_graphs(Zs,Ds,Qs,LAM):
    #We will now connect blocks of graphs
    newZ = [] #This will hold all D block domains
    newY = [] #This will hold all the Y values 
    
    for ind in range(0,len(Zs)):
        #ind denotes index of block number

        tempY =  getY(Qs[ind],Zs[ind],LAM)
        for y in tempY:
            newY = newY + [y]
        #We place ind'th Y block into Y
        
        for z in Zs[ind]:
            #x is now floating type
            newZ = newZ + [z + Ds[ind][0]]
           
    return [[newZ],[newY]]

def cavityABCD(cavity,LAM):
    optics = {
        'D':m_freeSpace,
        'L':m_lens,
        'M':m_mirror
    }
    
    rtm = get_RTM(unfold_Cav(cavity))

    return rtm


def cavityRun(cavity,LAM,res):
    #input cavity, lambda, and resolution of graph desired
    
    rtm = cavityABCD(cavity,LAM)
    if(check_stab(rtm)==True):
        print('This cavity is stable')

        #True if you want parameters to be shown at each optic

        q = get_q(rtm,LAM)
        #print('q at flat: ',q)
        Ds = getDs(cavity)
        Zs = getZs(Ds,res)
        Qs = getQs(q,cavity)
        

        if(Display_Parameters==True):
            printPar(Qs,cavity)


        #print('Qs: ',Qs)

        #Rs = getRs(Qs)
        #print('Rs: \n',Rs)

        Z,Y = combine_graphs(Zs,Ds,Qs,LAM)
        return [Z,Y]

    else:
        print('Cavity not stable')
        return [[],[]]
        
def Plot(Z,Y):
    plt.plot(Z,Y)
    plt.xlabel('Distance (cm)')
    plt.ylabel('Width (cm)')
    plt.title('Cavity Propogation provided by AOMoney and tntBizzle')
    plt.show()
    return 

#=========================================
#Interpretting cavity files

def getStringCav(sfile):
    cfile = open(sfile,'r+')
    scavity = cfile.readline()

    if(scavity[-1:]=='\n'):
        scavity = scavity[:-1]
        #Removes line break in case there is one
    while(scavity[-1]==' '):
        scavity = scavity[:-1]
        #Removes spaces at the end. Having one can cause us to 
        #never add the last optic/paramter pair

    return scavity
# get the first line of a file. Presumabley we're not idiots and 
# load the file with the cavity inside


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
#We'll need this to determine is a string represents a paramter or 
# not


def getCavity(scavity):
    cavity = []
    subcav = []
    # initialize cavity. add components with .append()

    skip = 0
    end = False #Indicates if we've reached the end
    
    for i in np.array(range(0,len(scavity))):
        #scanning individual characters

        if(scavity[i] == ' '):
            continue
            #ignore spaces

        if(skip > 0):
            #we wish to skip over elements already accounted for
            skip = skip -1 
            continue

        length = 0 #this will be how many characters it take to
        #reach a space
        while(scavity[i+length] != ' '):
            #Increase length until i+length gets a space
            #We must also check to see if we're at the end of the
            #line
            length = length + 1 
            if(i+length == len(scavity)):
                end = True
                break
                
        #Now scavity[i:i+length] is a block of strings between 2 ' '
        if(is_number(scavity[i:i+length])):
            subcav.append(float(scavity[i:i+length]))
            #submit a block of numbers to subcavity.

            if(end==True):
                cavity.append(subcav)
                #If we're at the last parameter, this will trigger
                #the last dump

        else:
            #This block contains letter(s)
            #We need to dump everything subcav has to cavity before 
            #appending subcav and then reset it

            cavity.append(subcav)
            subcav = []
            subcav.append(scavity[i:i+length])
            
        skip = length - 1 
        # set skip to 
        # continue over the rest of the characters
        # Example: We don't want 1.03 to give [1.03,.03,03,3]
        # skip = 3 will assure we read only the first 
        
    return cavity[1:]


def check_actions(cavity):
    # checks if cavity requires at least two graphs
    
    action = [0,0,0]
    
    for optic in cavity:
        
        if ((optic[0]=='Cx') or (optic[0]=='Cy')):
            action[0] = 1
            
        elif (optic[0]=='B'):
            action[1] = 1
        
        elif (optic[0] == 'X'):
            action[2] = 1
            
    return action
# possibly combine check_actions and split_cavity into one function?????    
def split_cavity(cavity):
    
    cavities = [cavity,cavity] #(x,y)
    
    i=0
    off = [0,0]

    while(i<len(cavity)):
    #for i in range(0,len(cavity)):

        #print('\n i:',cavity[i][0])
        #print('matrix i:',optics[cavity[i][0]](cavity[i][1]))

        if cavity[i][0] == 'Cx':        
            cavities[1] = remove_optic(i+off[1],cavities[1])
            #print('entered Cx and attempted remove_optic')
            off[1] = off[1] -1
        
        elif cavity[i][0] == 'Cy':
            cavities[0] = remove_optic(i+off[0],cavities[0])
            off[0] = off[0] -1

        elif cavity[i][0] == 'B':
            cavities[1] = remove_optic(i+off[1],cavities[1])
            off[1] = off[1] -1
            #print('Cavities[0]:',cavities[0],'/n Cavities[1]:',cavities[1])
        #print('off = ',off[0],' , ',off[1])

        i = i+1

    return cavities
    
def plot_xy(cavities,LAM,res):
    Zx,X = cavityRun(cavities[0],LAM,res)
    Zy,Y = cavityRun(cavities[1],LAM,res)
    
    plt.plot(Zx[0],X[0],label='X Axis')
    plt.plot(Zy[0],Y[0],label='Y Axis')

    #print('Zx:',Zx[0],'\n X:',X[0])

    plt.xlabel('Distance (cm)')
    plt.ylabel('Width (cm)')
    plt.title('Cavity Propogation provided by AOMoney and tntBizzle')

    plt.legend()
    
    plt.show()
    return           

def insert_optic(optic,par,pos,cavity):
    #Inserts ['optic',par] as the pos'th position in [cavity]
    cav1 = cavity[:pos]
    cav2 = cavity[pos:]
    return cav1 + [[optic,par]]+ cav2

def multiply_optic(pos,start,stop,step,cavity):
    # Add ['X',start,stop,step] after the pos'th optic in [cavity]
    return cavity[:pos+1] + [['X',start,stop,step]] + cavity[pos+1:]

def remove_optic(pos,cavity):
    #Removes pos'th optic
    return cavity[:pos]+cavity[pos+1:]




#=======================================
#Display parameters and all things numerical

def header():
    #Displays a header to organize the data
    print('\n=================================================================')
    print('Q \t |Width(mm) \t |Divergence(mrad) \t |Rayleigh Range')
    print('=================================================================\n')
    return

def getZr(q):
    #get Rayleigh Range
    return q.imag

def getWaist(qZ):
    #get Waist Width
    #You can input q or Zr, we'll figure it out
    global LAM
    global M2
    global n

    if(qZ.imag != 0):
        Zr = getZr(qZ)

    else:
        Zr = qZ

    #Now we are definitely dealing with Zr
    return(np.sqrt(4*LAM*M2/math.pi/n*Zr))

def getFullDiv(q):
    #get Full Angle of Divergence from q
    Zr = q.imag
    W0 = getWaist(Zr)
    return(W0/Zr)

def dispQPar(q,name):
    #print a line showing q paramters
    print(name, '\t',getSpot(q)*10,'\t',getDiv(q)*10,'\t',getZr(q))

def getR(q):
    #get radius of curvature
    return(1/((1/q).real))

def getSpot(q):
    #get spot size
    qover1 = 1/q
    return((np.sqrt(-LAM/np.pi/qover1.imag))*2)

def getDiv(q):
    #get divergence at q
    #We solve for the derivative of w, wp, and apply arctan to it

    Theta = getFullDiv(q)
    z0 = q.real
    w = getSpot(q)
    wp = Theta**2 * z0/w
    return(np.arctan(wp))

def getFFR(qx,qy):
    #get far-field roundness, the ratio of far-field divergence
    divx = getFullDiv(qx)
    divy = getFullDiv(qy)
    return(min(dicx,divy)/max(divx,divy))

def getAsy(qx,qy):
    Wx = getWaist(qx)
    Wy = getWaist(qy)
    return(max(Wx,Wy)/min(Wx,Wy))

def nameL(optic):
    return('Lens w/ FL='+str(optic[1]))

def nameM(optic):
    return('Mirror w/ ROC='+str(optic[1]))

def printPar(Qs,cavity):
    #Organizes the parameters and conducts the printing
    
    header()

    names = {
        'L':nameL,
        'M':nameM,
        'Cx':nameL,
        'Cy':nameL
    }
    Qind = 0
    for i in range(len(cavity)):
        #We'll only display parameters at boundaries of free space

        #print('cavity[',i,'][0] = ', cavity[i][0])
        
        if(cavity[i][0] == 'D'):
            #We now are looking exclusively at free space
            
            if(cavity[i-1][0] != 'D'):
                #Before free space is an optic
                
                if(cavity[i-1][0] == 'M'):
                    dispQPar(Qs[Qind],'After mirror w/ ROC='+str(cavity[i-1][1]))
                else:
                    dispQPar(Qs[Qind],'After lens  w/ FL='+str(cavity[i-1][1]))
                print('\n')

            if(cavity[i+1][0] != 'D'):
                q = prop_q(Qs[Qind],m_freeSpace(cavity[i][1]))
                if(cavity[i+1][0] == 'M'):
                    dispQPar(q,'Before  mirror  w/ ROC='+str(cavity[i+1][1]))
                else:
                    dispQPar(q,'Before  lens  w/ FL='+str(cavity[i+1][1]))
                   
            Qind = Qind +1
            #print('Qind = ',Qind)
        #End looking at 'D' components
    #End for loop
    
    return()

            
    
    


    
