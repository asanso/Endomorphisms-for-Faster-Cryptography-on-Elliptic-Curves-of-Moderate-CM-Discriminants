p = 2163160611951109656514578155686584121369935567340260530653195705987773102423967
aboldhat=1600224790493789893768878688270410785959510786167424698577648903928123990501736
bboldhat=442173884976592923728839681693322529291546262090176625621658004289727200834915
E0 = EllipticCurve(GF(p), [aboldhat,bboldhat])
t =E0.trace_of_frobenius()
disc = t^2 - 4*p
print(factor(disc))
print()
#list of curves on the crater
#Fp = GF(p)
#FpX = Fp['X']
#X = FpX.gen()
#D = fundamental_discriminant(-547)
#H = FpX(hilbert_class_polynomial(D))
#print(H.roots())
#print()
