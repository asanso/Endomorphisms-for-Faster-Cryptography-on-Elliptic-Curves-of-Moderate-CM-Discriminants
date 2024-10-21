p = 2^255 +3225
a =  p -3
b = 28091019353058090096996979000309560759124368558014865957655842872397301267595
E = EllipticCurve(GF(p), [a,b])
t = E.trace_of_frobenius()
disc = t^2 - 4*p)
factor(disc)
