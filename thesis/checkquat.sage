import unittest


def show_quaternion(expr, **kwargs) -> None:
    """If you're using a Jupyter notebook you can use this utility function to
    pretty-print the four coefficients of the quaternion, each on its own line.

    :param expr: the quaternion to print
    :param **kwargs: passed to c.subs for each of the four coefficients c
    """
    coeffs = [c.subs(**kwargs) for c in expr.coefficient_tuple()]
    if coeffs[1] == 0 and coeffs[2] == 0 and coeffs[3] == 0:
        show(coeffs[0])
    else:
        show(1, "." * 12, coeffs[0])
        show(LatexExpr(r"\sqrt{\epsilon}"), "." * 10, coeffs[1])
        show(LatexExpr(r"\Pi"), "." * 11, coeffs[2])
        show(LatexExpr(r"\sqrt{\epsilon}\Pi"), "." * 9, coeffs[3])


def project_to_trace_zero(expr):
    return expr - expr.coefficient_tuple()[0]


def hermitian_form(x, y):
    coeffs = (x * y.conjugate()).coefficient_tuple()
    return coeffs[0] + coeffs[1] * sqrt_eps


def symmetric_form(x, y):
    return (x * y.conjugate() + y * x.conjugate()) / 2


# Variables
R = QQ["s0", "s1", "t0", "t1", "a0", "a1", "b0", "b1", "eps", "varpi", "z0", "z1"]
s0, s1, t0, t1, a0, a1, b0, b1, eps, varpi, z0, z1 = R.gens()
DD = QuaternionAlgebra(Frac(R), eps, varpi, names=("sqrt_eps", "Pi", "sqrt_eps_Pi"))
sqrt_eps, Pi, sqrt_eps_Pi = DD.gens()

# Main variables
alpha = a0 + a1 * sqrt_eps
beta = b0 + b1 * sqrt_eps
s = s0 + s1 * sqrt_eps
t = t0 + t1 * sqrt_eps
lam_inv = z0 + z1 * sqrt_eps

# Their conjugates
alphac = a0 - a1 * sqrt_eps
betac = b0 - b1 * sqrt_eps
sc = s0 - s1 * sqrt_eps
tc = t0 - t1 * sqrt_eps
lamc_inv = z0 - z1 * sqrt_eps

r = 100


# The pair (g,u)
def g(x):
    return lam_inv * x * (alpha + beta * Pi)


u = s + t * Pi


# Derived quantities
uu = hermitian_form(u, u)
x = project_to_trace_zero(u.conjugate() * sqrt_eps * u)
y = project_to_trace_zero(varpi ^ r * (alpha + beta * Pi))
xx = symmetric_form(x, x)
xy = symmetric_form(x, y)
yy = symmetric_form(y, y)


class QuaternionTestCase(unittest.TestCase):
    def assertQtrnsEqualWhen(self, expr1, expr2, **kwargs):
        coeffs1 = [c.subs(**kwargs) for c in expr1.coefficient_tuple()]
        coeffs2 = [c.subs(**kwargs) for c in expr2.coefficient_tuple()]
        for i in range(4):
            self.assertEqual(coeffs1[i], coeffs2[i])

    def assertQtrnsEqual(self, expr1, expr2):
        self.assertQtrnsEqualWhen(expr1, expr2, s0=0, s1=0)
        self.assertQtrnsEqualWhen(expr1, expr2, t0=0, t1=0)

    def assertQtrnsEqualExactly(self, expr1, expr2):
        self.assertQtrnsEqualWhen(expr1, expr2)

    def test_uu(self):
        self.assertQtrnsEqualExactly(uu, s * sc - t * tc * varpi)

    def test_gu_u(self):
        self.assertQtrnsEqualWhen(
            hermitian_form(g(u), u), lam_inv * alphac * uu, s0=0, s1=0
        )
        self.assertQtrnsEqualWhen(
            hermitian_form(g(u), u), lam_inv * alpha * uu, t0=0, t1=0
        )

    def test_x(self):
        self.assertQtrnsEqual(x, (s * sc + t * tc * varpi) * sqrt_eps)

    def test_xx(self):
        self.assertQtrnsEqual(xx, -eps * uu**2)

    def test_xy(self):
        self.assertQtrnsEqual(xy, -varpi ^ r * a1 * eps * (s * sc + t * tc * varpi))

    def test_yy(self):
        self.assertQtrnsEqual(
            yy, varpi ^ (2 * r) * (-(a1**2) * eps - beta * betac * varpi)
        )

    def test_det(self):
        self.assertQtrnsEqual(
            xx * yy - xy**2,
            eps * uu ^ 2 * varpi ^ 200 * (varpi * beta * betac),
        )


if __name__ == "__main__":
    unittest.main()
