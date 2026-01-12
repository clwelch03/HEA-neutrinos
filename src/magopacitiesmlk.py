from math import sqrt, exp, log, log10, erf, cos, tanh, sinh, cosh, pi, floor
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.special import iv as besseli
import numpy as np

#constants and units, all calculations done in fm, e^2 / 4pi = alpha
SINTW = sqrt(0.231)
HBARC = 197.3269718
NEUTRON_MASS = 939.5653 / HBARC
PROTON_MASS = 938.272 / HBARC
MN = 2 * NEUTRON_MASS * PROTON_MASS / (NEUTRON_MASS + PROTON_MASS)
MSPLIT = PROTON_MASS - NEUTRON_MASS
MUON_MASS = 105.7 / HBARC
ELECTRON_MASS = 0.511 / HBARC
GA = 1.267
MEVFM3_TO_GCM3 = 1.78e12
FM4_TO_DYNECM2 = 3.16e35
G15_TO_FM = 19.5 * sqrt(4 * pi / 137) / HBARC**2
GP = 5.5858
GN = -3.8263
COSTC = sqrt(0.95)
GF = (HBARC / 292800)**2

##helper functions
#functions for chemical potentials
def nfd(e, mu, t): return 1 / (np.exp((e - mu) / t) + 1)

def n_of_mue(mue, eb, t):
    if not type(mue) == float:
        mue = mue[0]
    maxe = max([mue, 0]) + 10 * t
    nmax = int(floor(maxe**2 / (2 * eb))) + 1
    n = quad(lambda kz: 2 * nfd(np.sqrt(kz**2 + ELECTRON_MASS**2), mue, t) - 2 * nfd(np.sqrt(kz**2 + ELECTRON_MASS**2), -mue, t), 0, maxe)[0]
    for ne in range(1, nmax):
        n += quad(lambda kz: 4 * nfd(np.sqrt(kz**2 + 2 * ne * eb + ELECTRON_MASS**2), mue, t)
            - 4 * nfd(np.sqrt(kz**2 + 2 * ne * eb + ELECTRON_MASS**2), -mue, t), 0, sqrt(maxe**2 - 2 * ne * eb))[0]
    return n * eb / (4 * pi**2)

def mue_of_n(n, eb, t):
    return fsolve(lambda mue: n_of_mue(mue, eb, t) - n, t)[0]

#lepton kinematics
def thetapm(x, sp, sn, eb, ui, pm):
    if MSPLIT + pm * ELECTRON_MASS + ui - eb / (4 * MN) * (GN * sn - (GP - 2) * sp) - x > 0:
        return 1
    else:
        return 0

def e0pm(knu, sp, sn, eb, ui, pm):
    delta = MSPLIT + ui - eb / (4 * MN) * (GN * sn - (GP - 2) * sp)
    if knu + pm * delta > ELECTRON_MASS:
        return sqrt((knu + pm * delta)**2 - ELECTRON_MASS**2)
    else:
        return 0

#for pauli blocking in degenerate CC
def ufunc(eb, t): return MN * t / eb * (1 - exp(-eb / (MN * t)))
def f1(k0, u): return sqrt(pi) / (2 * u * sqrt(1 - u)) * (sqrt( 1 - u) * erf(k0)
    - exp(-k0**2 * u) * erf(k0 * sqrt(1 - u)))
def f2(k0, u): return 1 / (4 * u**2 * sqrt(1 - u)) * (sqrt(pi) * sqrt(1- u) * (2 + u) * erf(k0)
    - 2 * sqrt(pi) * exp(-k0**2 * u) * (1 + k0**2 * u) * erf(k0 * sqrt(1 - u))
    - 2 * k0 * u * exp(-k0**2) * sqrt(1 - u))
def f3(k0, u): return 1 / (4 * u * (1 - u)**(3 / 2)) * (sqrt(pi) * (1 - u)**(3 / 2) * erf(k0)
    - exp(-k0**2 * u) * sqrt(pi) * erf(k0 * sqrt(1 - u)) + 2 * u * sqrt(1 - u) * k0 * exp(-k0**2))
def g1(k0, u): return pi**2 * exp(-k0**2) / (12 * k0 * sqrt(1 - u)) * (sqrt(1 - u)
    - k0 * u * sqrt(pi) * erf(k0 * sqrt(1 - u)))
def g2(k0, u): return pi**2 * exp(-k0**2) / (12 * sqrt(1 - u)) * (k0 * sqrt(1 - u) 
    - exp(-k0**2 * (u - 1)) * (k0**2 * u - 1) * sqrt(pi) * erf(k0 * sqrt(1 - u)))
def g3(k0, u): return pi**2 * exp(-k0**2) / (24 * (1 - u)**(3 / 2)) * (2 * k0 * sqrt(1 - u)
    - exp(-k0**2 * (u - 1)) * u * sqrt(pi) * erf(k0 * sqrt(1 - u)))

def bigb(kfn, u, w, eb, t, knu, cost):
    kfnt = kfn / sqrt(2 * MN * t)
    knut = knu / sqrt(2 * MN * t)
    return 2 * MN * t / (eb * sqrt(pi)) * exp(eb / (2 * MN * t)) * (f1(kfnt, u) + w * g1(kfnt, u) 
        + knut**2 * ((f1(kfnt, u) + w * g1(kfnt, u)) * (u + cost**2 - u * cost**2)
        + (f2(kfnt, u) + w * g2(kfnt, u)) * u**2 * (cost**2 - 1)
        + (f3(kfnt, u) + w * g3(kfnt, u)) * (u**2 + 2 * cost**2 + u**2 * cost**2)))

#for nc stuff
def delta_sspr(eb, g, s, spr): return - g * (s - spr) * eb / (4 * NUCMASS)

def k0_func(eb, g, s, spr, q0, q): return abs(MN / q * (q0 - delta_sspr(eb, g, s, spr)))

##charged current matrix elements
def mred_cc(sp, sn, nezero, cost):
    if nezero:
        match (sn, sp):
            case (1, 1):
                return 1 / 2 * (1 + GA)**2 * (1 + cost)
            case (1, -1):
                return 0
            case (-1, 1):
                return 2 * GA**2 * (1 - cost)
            case (-1, -1):
                return 1 / 2 * (1 - GA)**2 * (1 + cost)
    else:
        match (sn, sp):
            case (1, 1):
                1 / 4 * ((1 + GA)**2 * (1 + cost) + (1 - GA)**2 * (1 - cost))
            case (1, -1):
                return GA**2 * (1 + cost)
            case (-1, 1):
                return GA**2 * (1 - cost)
            case (-1, -1):
                return 1 / 4 * ((1 + GA)**2 * (1 - cost) + (1 - GA)**2 * (1 + cost))

def vt_cc(x, sp, sn, cost):
    if x < 0:
        return 0
    elif x < 1 / 2:
        return mred_cc(sp, sn, True, cost)
    else:
        return mred_cc(sp, sn, True, cost) + (2 * x - 1) * mred_cc(sp, sn, False, cost)

##neutral current matrix elements
#protons
def mred_ncp(s, spr, cost, costpr, phi):
    sint = sqrt(1 - cost**2)
    sintpr = sqrt(1 - costpr**2)
    if s == spr:
        return 1 / 2 * ((1 - 4 * SINTW**2)**2 * (1 + cost * costpr + sint * sintpr * cos(phi))
            + GA**2 * (1 + cost * costpr - sint * sintpr * cos(phi)) 
            - 2 * GA * s * (1 - 4 * SINTW**2) * (cost + costpr))
    else:
        return GA**2 * (1 - cost * costpr)

#integrate over outgoing neutrino angles
def mred_ncp_int(s, spr, cost):
    if s == spr:
        return 1 / 2 * ((1 - 4 * SINTW**2)**2 + GA**2 - 2 * (1 - 4 * SINTW**2) * GA * s * cost)
    else:
        return GA**2

#neutrons
def mred_ncn(s, spr, cost, costpr, phi):
    sint = sqrt(1 - cost**2)
    sintpr = sqrt(1 - costpr**2)
    if s == spr:
        return 1 / 2 * (1 + cost * costpr + sint * sintpr * cos(phi) 
            + GA**2 * (1 + cost * costpr - sint * sintpr * cos(phi))
            - 2 * GA * s * (cost + costpr))
    else:
        return GA**2 * (1 - cost * costpr)

#integrate over outgoing neutrino angles
def mred_ncn_int(s, spr, cost):
    if s == spr:
        return 1 / 2 * (1 + GA**2 - 2 * GA * s * cost)
    else:
        return GA**2

#integrate over outgoing neutrino angles with exp(-q^2/4MT)
def mred_ncn_intexp(s, spr, cost, knu, t):
    u = knu**2 / (2 * MN * t)
    iu = 1 / tanh(u) - 1 / u
    prefactor = 1 / (4 * u) * ( 1 - exp(-u))
    if s == spr:
        return prefactor * (1 - iu + GA**2 + (1 + iu * (2 * cost**2 - 1))
            - 2 * GA * s * cost * (1 + iu))
    else:
        return prefactor * GA**2 * (2 - 2 * cost * iu)

#electrons
def mred_nce(h, hpr, cost, costpr, isx):
    if isx:
        xsign = -1
    else:
        xsign = 1

    if h != hpr:
        return 0
    elif h == 1:
        return 4 * SINTW**4 * (1 - cost) * (1 - costpr)
    else:
        return (1 + xsign * 2 * SINTW**2)**2 * (1 + cost) * (1 + costpr)

#just the helicity weighting
def mred_nce_helweights(h, hpr, isx):
    if isx:
        xsign = -1
    else:
        xsign = 1

    if h != hpr:
        return 0
    elif h == 1:
        return 4 * SINTW**4
    else:
        return (1 + xsign * 2 * SINTW**2)**2

## below here are all of the analytical approximations for opacities

#charged current nu + n -> e + p
def kappan(eb, t, mue, nb, yp, knu, cost, ui):
    kz = knu * cost
    kperp = knu * sqrt(1 - cost**2)
    ebmt = eb / (MN * t)
    prefactor = GF**2 * COSTC**2 * nb * (1 - yp) * eb / (pi * cosh(GN * ebmt / 4))
    #isospin chemical potential
    muhat = t * log(yp / (1 - yp) * (2 * cosh(GN * ebmt / 4) * sinh(ebmt / 2) / (ebmt * cosh(GP * ebmt / 4))))
    #this has weak magnetism and stimulated absorption
    corr = (1 + 1.1 * knu / MN) / (1 - nfd(knu, mue + muhat, t))

    spin_sum = 0
    for sp in [-1, 1]:
        for sn in [-1, 1]:
            ee = e0pm(knu, sp, sn, eb, ui, 1)
            elec_kin = (1 - nfd(ee, mue, t)) * thetapm(-knu, sp, sn, eb, ui, -1) * vt_cc(ee**2 / (2 * eb), sp, sn, cost)
            blocking = exp(GN * sn * ebmt / 4) * (1 - nb * yp * exp((GP - 2) * sp * ebmt / 4) / cosh(GP * ebmt / 4) 
                * (pi / (MN * t))**(3 / 2) * (1 - exp(-ebmt)) / (ebmt + 1 - exp(-ebmt)) * cosh(kz * ee / (2 * MN * t))
                * exp(-kperp**2 / (2 * MN * t) * (1 - exp(-ebmt)) / (ebmt + 1 - exp(-ebmt)) - (kz**2 + ee**2) / (4 * MN * t)))
            spin_sum += elec_kin * blocking
    return prefactor * spin_sum * corr

def kappan_degen(eb, t, mue, mun, nb, yp, knu, cost, ui):
    ebmt = eb / (MN * t)
    prefactor = GF**2 * COSTC**2 * nb * (1 - yp) * eb / pi 
    #isospin chemical potential
    muhat = t * log(nb * yp * sinh(ebmt / 2) * (2 * pi / (MN * t))**(3 / 2) / (ebmt * cosh(GP * ebmt / 4))) - mun
    #weak magnetism and stim absorption
    corr = (1 + 1.1 * knu / MN) / (1 - nfd(knu, mue + muhat, t))

    spin_sum = 0
    for sp in [-1, 1]:
        for sn in [-1, 1]:
            ee = e0pm(knu, sp, sn, eb, ui, 1)
            bb = bigb(sqrt(mun / t + GN * sn * ebmt / 4), (1 - exp(-ebmt) / ebmt), exp(GN * sn * ebmt / 8) / cosh(GN * ebmt / 8), eb, t, knu, cost)
            elec_kin = (1 - nfd(ee, mue, t)) * thetapm(-knu, sp, sn, eb, ui, -1) * vt_cc(ee**2 / (2 * eb), sp, sn, cost)
            blocking = 1 - yp / (1 - yp) * (1 - exp(-ebmt)) * exp((GP - 2) * sp * ebmt / 4) / cosh(GP * ebmt / 4) * bigb
            spin_sum += elec_kin * blocking
    return prefactor * spin_sum * corr

#charged current nubar + p -> e+ + n
def kappap(eb, t, mue, nb, yp, knu, cost, ui):
    kz = knu * cost
    kperp = knu * sqrt(1 - cost**2)
    ebmt = eb / (MN * t)
    prefactor = GF**2 * COSTC**2 * nb * yp * eb * exp(ebmt / 2) / (pi * cosh(GP * ebmt / 4))
    #isospin chemical potential
    muhat = t * log(yp / (1 - yp) * (2 * cosh(GN * ebmt / 4) * sinh(ebmt / 2) / (ebmt * cosh(GP * ebmt / 4))))
    #this has weak magnetism and stimulated absorption
    corr = (1 - 7.1 * knu / MN) / (1 - nfd(knu, -mue - muhat, t))

    spin_sum = 0
    for sp in [-1, 1]:
        for sn in [-1, 1]:
            ee = e0pm(knu, sp, sn, eb, ui, -1)
            elec_kin = (1 - nfd(ee, -mue, t)) * (1 - thetapm(knu, sp, sn, eb, ui, 1)) * vt_cc(ee**2 / (2 * eb), sp, sn, cost)
            blocking = exp((GP - 2) * sp * ebmt / 4) * (1 - nb * (1 - yp) * exp(GN * sn * ebmt / 4) / cosh(GN * ebmt / 4) 
                * (pi / (MN * t))**(3 / 2) * (1 - exp(-ebmt)) / (ebmt + 1 - exp(-ebmt)) * cosh(kz * ee / (2 * MN * t))
                * exp(-kperp**2 / (2 * MN * t) * (1 - exp(-ebmt)) / (ebmt + 1 - exp(-ebmt)) - (kz**2 + ee**2) / (4 * MN * t)))
            spin_sum += elec_kin * blocking
    return prefactor * spin_sum * corr

def kappap_degen(eb, t, mue, mun, nb, yp, knu, cost, ui):
    ebmt = eb / (MN * t)
    prefactor = GF**2 * COSTC**2 * nb * yp * eb / pi 
    #isospin chemical potential
    muhat = t * log(nb * yp * sinh(ebmt / 2) * (2 * pi / (MN * t))**(3 / 2) / (ebmt * cosh(GP * ebmt / 4))) - mun
    #weak magnetism and stim absorption
    corr = (1 - 7.1 * knu / MN) / (1 - nfd(knu, -mue - muhat, t))

    spin_sum = 0
    for sp in [-1, 1]:
        for sn in [-1, 1]:
            ee = e0pm(knu, sp, sn, eb, ui, -1)
            bb = bigb(sqrt(mun / t + GN * sn * ebmt / 4), (1 - exp(-ebmt) / ebmt), exp(GN * sn * ebmt / 8) / cosh(GN * ebmt / 8), eb, t, knu, cost)
            elec_kin = (1 - nfd(ee, -mue, t)) * (1 - thetapm(knu, sp, sn, eb, ui, 1)) * vt_cc(ee**2 / (2 * eb), sp, sn, cost)
            blocking = 1 - (1 - yp) / yp * (1 - exp(-ebmt)) * bigb
            spin_sum += elec_kin * blocking
    return prefactor * spin_sum * corr

#neutral current nu + n -> nu + n
def kappan_nc_diff(eb, t, n, knu, cost, knupr, costpr, phi):
    cos_kkpr = cost * costpr + sqrt((1 - cost**2) * (1 - costpr**2)) * cos(phi)
    q = sqrt(knu**2 + knupr**2 - 2 * knu * knupr * cos_kkpr)
    q0 = knupr - knu
    ebmt = eb / (MN * t)

    prefactor = GF**2 * knupr**2 * n / (4 * (2 * pi)**(5 / 2) * cosh(GN * ebmt / 4))

    spin_sum = 0
    for s in [-1, 1]:
        for spr in [-1, 1]:
            k0 = k0_func(eb, GN, s, spr, q0, q)
            factor = mred_ncn(s, spr, cost, costpr, phi) * exp(GN * s * ebmt / 4)
            blocking = exp(-k0**2 / (2 * MN * t)) - n * exp(GN * spr * ebmt / 4) / (sqrt(2) * cosh(GN * ebmt / 4)) \
                * pi**2 / (MN * t * q) * exp(-q**2 / (4 * MN * t)) \
                * (erf((q + 2 * k0) / (2 * sqrt(MN * t))) + erf((q - 2 * k0) / (2 * sqrt(MN * t))))
            spin_sum += factor * blocking
    return prefactor * spin_sum

def kappan_nc(eb, t, n, knu, cost):
    ebmt = eb / (MN * t)
    prefactor = GF**2 * n / (4 * pi * cosh(GN * ebmt / 4))

    spin_sum = 0
    for s in [-1, 1]:
        for spr in [-1, 1]:
            factor = (knu + delta_sspr(eb, GN, s, spr))**2 * exp(GN * s * ebmt / 4)
            blocking = mred_ncn_int(s, spr, cost) - n * exp(GN * spr * ebmt / 4) / cosh(GN * ebmt / 4) \
                * (pi / (MN * t))**(3 / 2) * mred_ncn_intexp(s, spr, cost, knu, t)
            spin_sum += factor * blocking
    return prefactor * spin_sum

def kappan_nc_diffdegen(eb, t, mu, knu, cost, costpr, phi):
    cos_kkpr = cost * costpr + sqrt((1 - cost**2) * (1 - costpr**2)) * cos(phi)
    q = sqrt(knu**2 + knupr**2 - 2 * knu * knupr * cos_kkpr)
    q0 = knupr - knu
    ebmt = eb / (MN * t)

    prefactor = GF**2 * knupr**2 * MN**2 / (64 * pi**5)
    spin_sum = 0
    for s in [-1, 1]:
        for spr in [-1, 1]:
            k0 = k0_func(eb, GN, s, spr, q0, q)
            e0 = k0**2 / (2 * MN * t) - GN * s * ebmt / 4
            spin_sum += mred_ncn(s, spr, cost, costpr, phi) * (q0 + t * log((1 + exp(e0 - (q0 + mu) / t)) / (1 + exp(e0 - mu / t)))) / (exp(q0 / t - 1))
    return prefactor * spin_sum

def kappan_nc_degen(eb, t, mu, knu, cost):
    ebmt = eb / (MN * t)

    prefactor = GF**2 * MN * t / (4 * pi**3) * sqrt(2 * mu * MN)

    spin_sum = 0
    for s in [-1, 1]:
        for spr in [-1, 1]:
            if knu + delta_sspr(eb, GN, s, spr) > 0:
                spin_sum += mred_ncn_int(s, spr, cost) * (knu + delta_sspr(eb, GN, s, spr))
    return prefactor * spin_sum

def kappap_nc_diff(eb, t, n, knu, cost, knupr, costpr, phi):
    cos_kkpr = cost * costpr + sqrt((1 - cost**2) * (1 - costpr**2)) * cos(phi)
    q = sqrt(knu**2 + knupr**2 - 2 * knu * knupr * cos_kkpr)
    q0 = knupr - knu
    ebmt = eb / (MN * t)
    qz = knupr * costpr - knu * cost
    cosq = qz / q
    qperp = sqrt(q**2 - qz**2)
    texp = exp(-ebmt)
    z = qperp**2 * sqrt(texp) / (eb * (1 - texp))
    zpr = qperp**2 * texp / (eb * (1 - texp**2))

    prefactor = GF**2 * knupr**2 * n / (4 * (2 * pi)**(5 / 2) * cosh(GP * ebmt / 4)) * sqrt(MN / t)
    
    spin_sum = 0 
    for s in [-1, 1]:
        for spr in [-1, 1]:
            dss = delta_sspr(eb, GP - 2, s, spr)
            alpha = int(round(MN * abs(q0 - dss) / eb))
            exponent = - (MN * (q0 - dss))**2 / (2 * qz**2 * MN * t)
            exppol = - (alpha * eb - MN * (q0 - dss))**2 / (2 * qz**2 * MN * t)
            
            polar_avg = exp(exponent) / sqrt(abs(qz) * q) + (1 - cosq**2) * besseli(alpha, z) / abs(qz) \
                * exp(exppol - qperp**2 / (2 * eb) * (1 + texp) / (1 - texp))
            polar_avg_blk = exp(2 * exponent) / sqrt(abs(qz) * q) + (1 - cosq**2) * besseli(alpha, zpr) / abs(qz) \
                * exp(exppol - qperp**2 / (2 * eb) * (1 + texp**2) / (1 - texp**2))
            blocking_pref = n * sinh(ebmt / 2) / cosh(GP * ebmt / 4) * 2 * pi**(3 / 2) / (eb * sqrt(MN * t))

            spin_sum += mred_ncp(s, spr, cost, costpr, phi) * exp((GP - 2) * s * ebmt / 4) * polar_avg \
                - exp((GP - 2) * spr * ebmt / 4) * blocking_pref * polar_avg_blk
    return prefactor * spin_sum

def kappap_nc(eb, t, n, knu, cost):
    ebmt = eb / (MN * t)
    texp = exp(-ebmt)

    prefactor = GF**2 * n / (4 * pi * cosh(GP * ebmt / 4) * sqrt(texp))

    spin_sum = 0
    for s in [-1, 1]:
        for spr in [-1, 1]:
            dss = delta_sspr(eb, GP - 2, s, spr)
            if knu + dss > 0:
                sum_pref = mred_ncp_int(s, spr, cost) * (knu + dss)**2 * exp((GP - 2) * ebmt / 4)
                blocking = 1 - n * sinh(ebmt / 2) / cosh(GP * ebmt / 4) * (2 * pi)**(3 / 2) / (eb * sqrt(mt)) \
                    * exp((GP - 2) * spr * ebmt / 4)
                spin_sum += sum_pref * blocking
    return prefactor * spin_sum

#functions for electron integrals
def ifunc(qbar, mubar, tbar):
    if qbar / tbar > 100: return 0
    elif (qbar - mubar) / tbar > 50:
        return log(1 + exp(mubar / tbar)) / (exp(qbar / tbar) - 1)
    elif qbar > 0:
        return log((1 + exp(mubar / tbar)) / (1 + exp((mubar - qbar) / tbar))) / (exp(qbar / tbar) - 1)
    else:
        return log((1 + exp((mubar + 2 * qbar) / tbar)) / (1 + exp((mubar + qbar) / tbar))) / (exp(qbar / tbar) - 1)

def jint(kperp, kbar, mubar, tbar, hc):
    q0 = (kperp**2 - kbar**2 * (1 - hc**2)) / (2 * kbar * (1 + hc))
    return kperp * (1 - hc) * (1 - (hc * kbar - q0) / sqrt(kperp**2 + (hc * kbar - q0)**2)) \
        * exp(-kperp**2 / 2) * besseli(0, kperp * kbar * sqrt(1 - hc**2)) * ifunc(q0, mubar, tbar)

def jfunc(knu, mue, t, hc, eb):
    kbar = knu / sqrt(eb)
    mubar = mue / sqrt(eb)
    tbar = t / sqrt(eb)
    maxe = max([0, mubar]) + 5 * tbar

    return quad(lambda kp: jint(kp, kbar, mubar, tbar, hc), 0, maxe)[0]

def kappae_nc(eb, t, mue, knu, cost, isx):
    if isx:
        xsign = -1
    else:
        xsign = 1
    
    prefactor = eb**2 * t / (8 * pi**3) * exp(-knu**2 * (1 - cost**2) / (2 * eb))
    hel_sum = 4 * SINTW**4 * jfunc(knu, mue, t, cost, eb) + (1 + 2 * xsign * SINTW**2)**2 * jfunc(knu, mue, t, -cost, eb)
    return prefactor * hel_sum