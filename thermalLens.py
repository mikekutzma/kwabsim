import cavitySolver as cs
import numpy as np
import matplotlib.pyplot as plt


#More polishing: 
#Element before 'X' shouldn't need a paramter
#Should deny cavities with more than one 'X'


def xDecipher(cavity,LAM,res):
    #input cavity, lambda, and resolution. This function is 
    #cavityRun but can interpret X as a multiplier. 

    # ['X',start,end,step]
    
    xind = 0

    while (cavity[xind][0] != 'X' ):
        # sweep, looking for 'X'
        
        xind = xind +1

        if((xind + 1 == len(cavity)) & (cavity[xind][0] != 'X')):
            #We are on the last element and still haven't found 
            #a fucking 'X'. We aint gun

            return cs.cavityRun(cavity,LAM,res)

            #now cavity[xind][0] == 'X'
    
    start,end,step = [cavity[xind][1],cavity[xind][2],cavity[xind][3]]

    if(xind == 0):
        #X is the first element. No target indicated. User is a fag
        print('X multiplier listed as first element. X element ignored.')
        return cs.cavityRun(cavity[1:],LAM,res)

    else: 
        before_x = cavity[:xind]
        after_x = cavity[xind+1:]
        cavities = before_x + after_x
        #print('cavities = ',cavities)
        #print('xind = ',xind)
        Xs = []
        Ys = []

        for i in np.linspace(start,end,( abs((end-start))/step) +1):
            # i is the paramter for the element before
            
            cavities[xind-1][-1]= i # replace paramter in X's target
            #with what X demands
        
            if(cs.check_stab(cs.get_RTM(cavities)) != False):
                X,Y = cs.cavityRun(cavities,LAM,res)
                Ys = Ys + Y
                Xs = Xs + X
            else:
                print(cavities[xind-1][0],' = ',cavities[xind-1][-1],'provides an unstable cavity')
                
    return [Xs,Ys]

def PlotXs(Xs,Ys):
    for i in np.array(range(0,len(Ys))):
        #print('len(Ys[i])',len(Ys[i]))
        plt.plot(Xs[i],Ys[i])
    plt.show()
    return
