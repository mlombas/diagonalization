from polynomials import *

p1 = Poly.from_str("3x^2 + 4x^3 + 2x^4 + -3x^6 + 4")
p2 = Poly.from_str("2x^5 + 3x^4 + 6x^10 + 5x + 3")
print(Term(2,4))
print(p1)
print(p2)
print(p1 + p2)
print(p1 - p2)
print(p1 * p2)
print(p1.derivate())
print(p2.derivate())
print(Poly.copy(p1) == p1)
print(p1.all_terms())

print(square_free(Poly.from_str("3x^800 + -2x^400")))

p3 = Poly.from_str("3x^2 + 2x + -5")
roots = calc_roots(p2)
print(roots)
print([abs(p2.evaluate(r)) for r in roots])
