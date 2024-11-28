import sage.all
from sage.rings.integer_ring import Z as ZZ
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.modules.free_module_element import free_module_element as vector
from curve import WeierstrassCurve
from curve import PointWeierstrass

class  Lollipop489201(WeierstrassCurve):

       def __init__(self, p, a, b, r, cofactor,L):
            super().__init__(p, a,b, r, cofactor)
            self.D = -547
            self.L = L
            self.cofactor = cofactor
            self.r = r
            M = Matrix([[-L,1], [r,0]])
            self.N = M.LLL()
            self.N_inv = self.N**-1

class Lollipop489201Point(PointWeierstrass):
      
        def fast_scalar_mul(self, n):
            psiP = self.psi()
            beta = vector([n,0]) * self.curve.N_inv
            b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
            k1 = n-b[0]
            k2 = -b[1]
            return self.multi_scalar_mul(k1, psiP, k2)