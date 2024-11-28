import sage.all
from sage.rings.integer_ring import Z as ZZ
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.modules.free_module_element import free_module_element as vector
from curve import WeierstrassCurve
from curve import PointWeierstrass

class  Lollipop489201(WeierstrassCurve):

       def __init__(self, p, a, b, r, cofactor):
            super().__init__(p, a,b, r, cofactor)