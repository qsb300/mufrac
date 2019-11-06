#!/usr/bin/python

#formated with python yapf --style=pep8
from math import log
from math import exp
from math import pi
'''
Calculate excess chemical potential and allowed fraction analytically.
'''


class Beta(object):
    "Obtain beta=1/T from different Temprature."

    def __init__(self, val, mode='T'):
        "Beta from mode= 'T'(reduced), 'K'(kelvin), 'C'(celsius)"
        self.beta = None
        if mode == 'T':
            self.fromT(val)
        elif mode == 'K':
            self.fromK(val)
        elif mode == "C":
            self.fromK(val + 273.15)
        else:
            raise NotImplementedError("mode=%s is unknown" % (mode))

    def fromT(self, T=1.0):
        "From reduced temperature"
        self.beta = 1.0 / T

    def fromK(self, tempK=298.15):
        "From kelvin temperature"
        #Avogadro constant N_A = 6.022045000e+23;
        #http://en.wikipedia.org/wiki/Boltzmann_constant 3.2976230(30)e-24 cal/K
        self.beta = tempK * 6.022045000e+23 * 3.297623030e-24 / 1000.0


class Particle(object):
    "Abstract base class."

    def __init__(self, phi=None):
        self._bmu = None
        self._phi = None
        if phi is not None:
            self._phi = phi
            self.update()
            self.check()

    def update(self):
        "Child shoudl update self._bmu via func(self._phi,...)"
        raise NotImplementedError("Must override update.")

    def check(self):
        "Check which member is not initialized."
        members = [
            attr for attr in dir(self)
            if not callable(getattr(self, attr)) and not attr.startswith("__")
        ]
        err = []
        #print members
        for m in members:
            if getattr(self, m) is None:
                err.append("%s not initialized in class: %s" %
                           (m, type(self).__name__))
        if len(err) != 0:
            raise RuntimeError("\n".join(err))
        #print err

    def bmu(self, phi=None):
        "Obtain excess chemical potential."
        if phi is not None:
            self._phi = phi
            self.update()
        if self._bmu is None:
            self.check()
        else:
            return self._bmu

    def mu(self, phi=None, beta=1.0):
        "Obtain excess chemical potential at beta=1/T."
        if phi is not None:
            self._phi = phi
            self.update()
        if self._bmu is None:
            self.check()
        else:
            return self._bmu / beta

    def frac(self, phi=None):
        "Obtain allowed fraction."
        if phi is not None:
            self._phi = phi
            self.update()
        if self._bmu is None:
            self.check()
        else:
            #print self._bmu
            return exp(-self._bmu)


class Spt(Particle):
    '''excess chemical potential calculated based on scaled particle theory.
       See DOI: 10.1002/prot.22425 eq. 3 for detail.
    '''

    def __init__(self, phi=None, zeta=1.0):
        self._zeta = zeta
        super(Spt, self).__init__(phi)

    def update(self):
        phi = self._phi
        zeta = self._zeta
        sqrzeta = pow(zeta, 2)
        cubzeta = pow(zeta, 3)
        phirate = phi / (1.0 - phi)
        sqrphirate = pow(phirate, 2)
        cubphirate = pow(phirate, 3)
        self._bmu=-log(1.0-phi) \
                  + (3.0*zeta+3.0*sqrzeta+cubzeta)*phirate \
                  + (9.0*sqrzeta/2.0+3.0*cubzeta)*sqrphirate \
                  + 3.0*cubzeta*cubphirate


class Pyv(Particle):
    '''Percus-Yevick formula via virial route.
        See DOI: 10.1063/1.4968039 Table II for detail.
	    (1+\eta+\eta^2)/(1-\eta)^3, ?lower bound
    '''

    def update(self):
        phi = self._phi
        self._bmu=2.0*log(1.0-phi) \
                  + 2.0*phi*(5.0-2.0*phi)/pow((1.0-phi),2)


class Pyc(Particle):
    '''Percus-Yevick formula via compressibility route.
        See DOI: 10.1063/1.4968039 Table II for detail.
	    See DOI: 10.1007/BF01009789 "first derived from" scaled partical theory.
        :math:`\\alpha`.
	    (1+2\eta+3\eta^2)/(1-\eta)^2, ?upper bound
        .. math::
            (1+2\eta+3\eta^2)/(1-\eta)^2
    '''

    def update(self):
        phi = self._phi
        self._bmu=-log(1.0-phi) \
                  + phi*(14.0-13.0*phi+5.0*pow(phi,2))/(2.0*pow((1.0-phi),3))


class Pyu(Particle):
    '''Percus-Yevick formula via chemical-potential.
        See DOI: 10.1063/1.4968039 Table II for detail.
    '''

    def update(self):
        phi = self._phi
        self._bmu=-log(1.0-phi) \
                  + phi*(14.0+phi)/(2.0*pow((1.0-phi),2))


class Cs(Particle):
    '''Carnahan-Starling formula.
	    See DOI: 10.1063/1.4968039 Table II for detail.
	'''

    def update(self):
        phi = self._phi
        self._bmu = phi * (8.0 - 9.0 * phi + 3.0 * pow(phi, 2)) \
                    / pow((1.0 - phi), 3)


class Csk(Particle):
    '''Carnahan-Starling-Kolafa formula.
        See DOI: 10.1063/1.4968039 Table II for detail.
	'''

    def update(self):
        phi = self._phi
        self._bmu=5.0/3.0*log(1.0-phi) \
                  + phi*(58.0-79.0*phi+39.0*pow(phi,2)-8.0*pow(phi,3)) \
                  / (6.0*pow((1-phi),3))


class Virial(Particle):
    '''Virial coeffcient.
	    See DOI: 10.1063/1.1733209 eq. 3 for detail.
	'''
    #www.sklogwiki.org/SklogWiki/index.php/Carnahan-Starling_equation_of_state
    #From Clisby and McCoy, related to number density (N/V) by (pi/6)^(n-1)
    vtCM = [
        4.0, 10.0, 18.3647684, 28.224512, 39.8151475, 53.3444198, 68.5375488,
        85.8128384, 105.775104
    ]

    #Bn=n^2+n-2 (?CS)
    #vt=[4,10,18,28,42,54,74,88,108]

    def __init__(self, phi=None, uplim=10, vt=None):
        self._uplim = uplim
        if vt is None:
            self.vt = Virial.vtCM
        else:
            if vt == "CS":
                self.vt = [i**2 + i - 2 for i in range(2, uplim + 2)]
            elif vt == "CM":
                self.vt = Virial.vtCM
            else:
                self.vt = vt
        super(Virial, self).__init__(phi)

    def update(self):
        phi = self._phi
        self._bmu = 0.0
        for i in range(0, 9):
            if i + 2 <= self._uplim:
                self._bmu += (i + 2.0) \
                             / ( i + 1.0) * self.vt[i] * pow(phi, i + 1)
                #print self._bmu
            else:
                break


class GfmtPyc(Particle):
    '''Generalized fundamental measure theory. 
       See DOI:10.1103/PhysRevE.81.031919 eq. 1-5 for detail.
	'''
    pass


class GfmtCs(Particle):
    '''Generalized fundamental measure theory. 
       See DOI:10.1103/PhysRevE.81.031919 eq. 1-5 for detail.
    '''
    pass


#phi concetration, zeta Rp/RC
def SPT(phi, zeta):
    sqrzeta = zeta * zeta
    cubzeta = zeta * zeta * zeta
    phirate = phi / (1 - phi)
    sqrphirate = phirate * phirate
    cubphirate = phirate * phirate * phirate
    bmu=-log(1-phi)\
         + (3*zeta+3*sqrzeta+cubzeta)*phirate\
         + (9*sqrzeta/2.0+3*cubzeta)*sqrphirate\
         + 3*cubzeta*cubphirate
    return bmu


def GFMTSph(phi, lp, sp, vp, Rc):
    '''generalized fundamental measure theory with spherical crowder. 
       See DOI:10.1103/PhysRevE.81.031919 eq. 1-5 for detail.
	'''


def GFMTPY(phi, lp, sp, vp, lc, sc, vc):
    "generalized fundamental measure theory with PY."


def GFMTCS(phi, lp, sp, vp, lc, sc, vc):
    "generalized fundamental measure theory with CS."


def CS(phi):
    '''Carnahan-Starling formula.
	    See DOI: 10.1063/1.4968039 Table II for detail.
	'''
    bmu = phi * (8.0 - 9.0 * phi + 3.0 * pow(phi, 2)) / pow((1.0 - phi), 3)
    return bmu


def CSK(phi):
    '''Carnahan-Starling-Kolafa formula.
        See DOI: 10.1063/1.4968039 Table II for detail.
	'''
    bmu=5.0/3.0*log(1.0-phi) \
        + phi*(58.0-79.0*phi+39.0*pow(phi,2)-8.0*pow(phi,3)) \
        / (6.0*pow((1-phi),3))
    return bmu


def PhiScale(method, ratio, **kwargs):
    dct = {}
    for key, value in kwargs.iteritems():
        dct[key] = value
    print dct
    if len(dct) == 1:
        return method(dct['phi'] * ratio)
    elif len(dct) == 2:
        return method(dct['phi'] * ratio, dct['zeta'])
    '''
	if kwargs.has_keys('phi'):
		if len(kwargs)==1:
			return method(kwargs['phi'])
		else:
			if kwargs.has_keys('zeta'):
				return method(kwargs['phi'],kwargs['zeta'])
	'''


def testPhi(models, phi):
    print "SPT"
    bm = PhiScale(SPT, 1.0, phi=phi, zeta=1.0)
    print exp(-bm)
    print "CS"
    bm = PhiScale(CS, 1.0, phi=phi)
    print exp(-bm)

    p = Particle()
    #raise NotImplementedError
    #p.frac(phi)
    s = Spt()
    #raise RuntimeError
    #print s.frac()

    for m in models:
        a = m(phi)
        b = m()
        print m.__name__
        print b.frac(phi), exp(-a.bmu()), a.bmu(), a.mu(beta=1.2) * 1.2

    print "Virial CM 10"
    v10 = Virial(phi, 10)
    print v10.frac(), v10.bmu()
    v10CM = Virial(phi, 10, "CM")
    print v10CM.frac(), v10CM.bmu()
    print "Virial CM 2"
    v2 = Virial(phi, 2)
    print v2.frac(), v2.bmu()
    print "Virial CM 3"
    v3 = Virial(phi, 3)
    print v3.frac(), v3.bmu()
    print "Virial CM 4"
    v4 = Virial(phi, 4)
    print v4.frac(), v4.bmu()
    print "Virial CS vt"
    vtcs = [i**2 + i - 2 for i in range(2, 22)]
    v100 = Virial(phi, 20, vtcs)
    print v100.frac(), v100.bmu()
    print "Virial CS"
    vscs = Virial(phi, 20, "CS")
    print vscs.frac(), vscs.bmu()
    print "Virial CS 30"
    vscs10 = Virial(phi, 30, "CS")
    print vscs10.frac(), vscs10.bmu()


def testModel(models, start, end, dx):
    for i in range(start, end):
        phi = i * dx
        out = [phi]
        for m in models:
            a = m(phi)
            out.append(a.bmu())
        print ' '.join(["%8.4f" % u for u in out])


def testVirial(start, end, dx):
    for i in range(start, end):
        phi = i * dx
        out = [phi]
        uplims = [10, 2, 3, 4]
        for u in uplims:
            a = Virial(phi, u)
            out.append(a.bmu())
        print ' '.join(["%8.4f" % u for u in out])

def testViriaVts(start, end, dx, vts):
    for i in range(start, end):
        phi = i * dx
	out = [phi]
	nvt=len(vts)+1
	a=Virial(phi, nvt, vts)
	out.append(a.bmu())
    	print ' '.join(["%8.4f" % u for u in out])


if __name__ == "__main__":
    import sys
    models = [Spt, Virial, Cs, Csk, Pyv, Pyc, Pyu]
    if len(sys.argv) == 2:
        phi = float(sys.argv[1])
        testPhi(models, phi)
    elif len(sys.argv) == 4:
        start = int(sys.argv[1])
        end = int(sys.argv[2])
        dx = float(sys.argv[3])
        testModel(models, start, end, dx)
    elif len(sys.argv) == 5:
        start = int(sys.argv[1])
        end = int(sys.argv[2])
        dx = float(sys.argv[3])
        testVirial(start, end, dx)
    elif len(sys.argv) > 5:
        start = int(sys.argv[1])
	end = int(sys.argv[2])
	dx = float(sys.argv[3])
	vts=map(float,sys.argv[4:])
	testViriaVts(start, end, dx, vts)
