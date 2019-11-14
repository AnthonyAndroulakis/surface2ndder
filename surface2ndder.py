#surface2ndder.py, Anthony Androulakis, 2019
#uses 2nd derivative test to find the local min, local max, and saddle points of a surface given a function
#input: function
#output: a dictionary containing results for each critical point
#example usage: python3 surface2ndder.py x**3-12*x*y+8*y**3

import sympy as sym
import sys

function=sys.argv[1]

x,y=sym.symbols('x y', real=True)

#find 1st derivative
dx=eval('sym.diff('+function+',x)')
dy=eval('sym.diff('+function+',y)')

#find critical points
allcrtpts=list(sym.nonlinsolve([dx,dy],[x,y]))
crtpts=[]
for i in range(len(allcrtpts)):
    if 'I' not in str(allcrtpts[i]):
        crtpts.append(allcrtpts[i]) #only real solutions, critical points found here from setting 1st derivative to 0


#plugging in critical points to second partial derivatives, dxx, dyy, and dxy
dxx=eval('sym.diff('+str(dx)+',x)')
dxxvals=[]
for j in range(len(crtpts)):
    dxxvals.append(dxx.evalf(subs={x: crtpts[j][0],y: crtpts[j][1]}))

dyy=eval('sym.diff('+str(dy)+',y)')
dyyvals=[]
for k in range(len(crtpts)):
    dyyvals.append(dyy.evalf(subs={x: crtpts[k][0],y: crtpts[k][1]}))

dxy=eval('sym.diff('+str(dx)+',y)')
dxyvals=[]
for l in range(len(crtpts)):
    dxyvals.append(dxy.evalf(subs={x: crtpts[l][0],y: crtpts[l][1]}))

#finding determinant values (dxx*dyy-dxy^2)
determinants=[]
for m in range(len(crtpts)):
    determinants.append(dxxvals[m]*dyyvals[m]-dxyvals[m]**2)

#finding results
results={}
for n in range(len(crtpts)):
    if determinants[n]>0: #positive, so yay
        if dxxvals[n]>0:
            results[str(crtpts[n])]='local min' #flipped b/c 2nd derivative shows rate of change of slope
        elif dxxvals[n]<0:
            results[str(crtpts[n])]='local max'
    elif determinants[n]<0: #negative, so saddle point
        results[str(crtpts[n])]='saddle point'
    elif determinants[n]==0: #just 0, so no
        results[str(crtpts[n])]='inconclusive using 2nd derivative test'

print(results)
