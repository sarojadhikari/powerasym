from posterior import PosteriorfNL
import numpy as np
from scipy import integrate
from scipy.stats import norm

pf=PosteriorfNL()

class PosteriorSDfNL(PosteriorfNL):
    """
    scale dependent fNL=fNL0(l/60)^nfNL
    """
    def __init__(self, ns=0.965, N=50.0, l0=60.0):
        self.ns=ns
        self.efolds=N
        self.l0=l0
        self.A0const=7.94E-10
        self.A1const=0.0258/500. # for ns=0.965
    
    def fNLl(self, l, fNL60=50., nfNL=-1.0):
        return fNL60*(l/self.l0)**nfNL
    
    def pdfl(self, l, A, nA, fNL60, nfNL):
        self.sigmaG=2.84/(2.*l+1)
        return self.pdf(self.fNLl(l, A, nA), self.fNLl(l, fNL60, nfNL))
    
    def pdf_multipoles(self, lmin, lmax, A60, nA, f60, nf):
        pdfm=0.1
        for l in np.arange(lmin, lmax+1):
            pdfm=pdfm*self.pdfl(l, A60, nA, f60, nf)
        
        return pdfm
    
    def post_multipoles(self, lmin, lmax, A60, nA):
        """
        use the combined likelihood of different multipoles given a measurement of
        (A, nA) for fNL0 and nfNL, where 
            fNL(l) = fNL0(l/60.)^nfNL, and
            A(l) = A60(l/60.)^nA
        
        to get a 2D posteriro pdf in the plane (fNL, nfNL)
        """
        lmodes=lmax-lmin+1
        Ndim=100
        f60list=np.arange(0.0, 500., 500./Ndim)
        nfNLlist=np.arange(-2.0, 2.0, 4./Ndim)
        posterior2d=np.ones((len(f60list), len(nfNLlist)))
        print np.shape(posterior2d)
        
        post1 = lambda fNL,nfNL: self.pdf_multipoles(lmin, lmax, A60, nA, fNL, nfNL)
        postcombined = lambda fNL,nfNL: self.pdf_multipoles(lmin, lmax, A60, nA, fNL, nfNL) * self.scale_dep_fNL(50., 50., -0.3, 0.3, fNL, nfNL)
        
        #norm2d = integrate.nquad(postcombined, [[0.0, 1000.],[-5.0, 5.0]])[0]
        norm2d=postcombined(50., 0.0)
        print norm2d
        
        for f in range(len(f60list)):
            for n in range(len(nfNLlist)):
                posterior2d[f][n]=postcombined(f60list[f], nfNLlist[n])

        return f60list, nfNLlist, posterior2d/norm2d
        
    def scale_dep_fNL(self, fNL60, sfNL60, nfNL, snfNL, fNLvalue, nfNLvalue):
        """
        return unnormalized 2d pdf distribution for fNL, nfNL (assuming Gaussian) and errors given
        also take fNL distribution to be folded normal for the purpose of combing with power
        asymmetry constraints
        """
        post = lambda fnl,nf: self.pdf_fold(fnl, mean=fNL60, sigma=sfNL60) * norm.pdf(nf, loc=nfNL, scale=snfNL)
        return post(fNLvalue, nfNLvalue)
        