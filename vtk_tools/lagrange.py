import sympy as sy 

def LagrangPoly(x,order,i,xi=None):
    if xi==None:
        xi=sy.symbols('x:%d'%(order+1))
    index = range(order+1)    
    return sy.prod([(x-xi[j])/(xi[i]-xi[j]) for j in index if j != i])
