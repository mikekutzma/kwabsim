import thermalLens as tl

#X,Ys = tl.xDecipher([['M',1],['D',.49],['L',1],['X',0,1,.25],['D',.49],['M',1]],1064*10**-5,200)
X,Y = tl.xDecipher([['M',0],['D',25.6],['M',-40],['D',2.5],['L',18],['D',2.5],['M',0],['D',36.6],['M',200]],1064*10**(-5),100)



#print('len(X): ',len(X),'len(Y): ',len(Y))
#print('len(Y[0]) = \n',len(Y[0]))


tl.PlotXs(X,Y)

#negative curvatures in X provide no results. Check linspace line?
