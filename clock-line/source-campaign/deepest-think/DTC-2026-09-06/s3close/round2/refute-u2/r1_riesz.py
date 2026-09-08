import mpmath as mp
mp.mp.dps = 50
S3 = 2*mp.pi**2
def Jser(t, N=200):          # small t series: J = sum_{k>=0} c_k t^{2k}
    # (1+t^2)L-2t = 2*sum_{m>=1} t^{2m+1}*(1/(2m+1)+1/(2m-1)); J = that/(2t^3)
    s = mp.mpf(0)
    for m in range(1, N):
        s += t**(2*m-2)*(mp.mpf(1)/(2*m+1)+mp.mpf(1)/(2*m-1))
    return s
def J(t):
    if t < mp.mpf('0.3'):  return Jser(t)
    if t > mp.mpf('3'):    return Jser(1/t)/t**4
    return ((1+t*t)*mp.log((1+t)/abs(1-t)) - 2*t)/(2*t**3)
f = lambda u: S3*(mp.e**u)**3*J(mp.e**u)
tot = mp.mpf(0)
pts = [-80,-40,-20,-10,-5,-2,-1,-0.5,-0.2,-0.05,0,0.05,0.2,0.5,1,2,5,10,20,40,80]
for a,b in zip(pts[:-1],pts[1:]):
    tot += mp.quad(f,[a,b],maxdegree=12)
print('I(|x|=1) =', mp.nstr(tot,20))
print('pi^4/2   =', mp.nstr(mp.pi**4/2,20))
print('rel      =', mp.nstr(abs(tot-mp.pi**4/2)/(mp.pi**4/2),6))
# check J series vs closed form in the overlap
for t in ['0.29','0.3','3.0','3.1']:
    t=mp.mpf(t)
    cf=((1+t*t)*mp.log((1+t)/abs(1-t)) - 2*t)/(2*t**3)
    sr= Jser(t) if t<1 else Jser(1/t)/t**4
    print('overlap t=',t, mp.nstr(cf,15), mp.nstr(sr,15))
