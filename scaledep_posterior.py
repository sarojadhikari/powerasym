from posterior import PosteriorfNL

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
        return fNL0*(l/self.l0)**nfNL
    
    def pdfl(self, l, A, nA, fNL60, nfNL):
        self.sigmaG=2.84/(2.*l+1)
        return self.pdf(self.fNLl(l, A, nA), self.fNLl(l, fNL0, nfNL))
    
    def pdf_multipoles(self, lmin, lmax, A60, nA):
        """
        use the combined likelihood of different multipoles given a measurement of
        (A, nA) for fNL0 and nfNL, where 
            fNL(l) = fNL0(l/60.)^nfNL, and
            A(l) = A60(l/60.)^nA
        
        to get a 2D posteriro pdf in the plane (fNL, nfNL)
        """
        
        