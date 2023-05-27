def softrepulsion(r,rmin,rmax,epsilon,rcut,ron=None):
    if ron is None: 
        ron = 0.95*rcut

    if(r<ron):
        S=1
        dS=0
    elif(ron<=r and r<=rcut):
        S=((rcut**2-r**2)**2)*(rcut**2 + 2*r**2 - 3*ron**2)/((rcut**2-ron**2)**3)
        dS=((12*(r**5))-(12*(rcut**2)*(r**3))-(12*(ron**2)*(r**3))+(12*(rcut**2)*(ron**2)*r))/((rcut**2-ron**2)**3)  
    elif(r>rcut):
        S=0
        dS=0

    if(r<rcut):
        U=epsilon*(1-(r/rcut)**4)
        dU=epsilon*(-4*(r**3)/(rcut**4))
    else:
        U=0
        dU=0

    shift=0

    if(ron<rcut):
        V=S*U
        F=-(S*dU + U*dS)
    elif(ron>=rcut):
        V=U-shift
        F=-dU

    return(V,F)
