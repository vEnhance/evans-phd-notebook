# Setup

import argparse
import unittest

q = var("q")
qs = var("q_s")  # = q^s


def irange(start, stop):
    return range(start, stop + 1)


def O(r, vb, vc, ve, vda):
    """Value of the overall orbital integral in the semi-Lie case"""
    assert vb % 2 != vc % 2
    assert r >= 0
    assert vb + vc >= 0
    S = 0
    for k in irange(-vb - r, 2 * ve + vc + r):
        n = min(
            (k + (vb + r)) // 2,
            (2 * ve + vc + r - k) // 2,
            ve,
            min((vb + vc) // 2, vda) + r,
        )
        S += (-qs) ** k * sum([q**i for i in irange(0, n)])
    if vda < ve - r and vb + vc > 2 * vda:
        for k in irange(2 * vda - vb + r, 2 * ve + vc - 2 * vda - r):
            c = min(k - (2 * vda - vb + r), 2 * ve + vc - 2 * vda - r - k, ve - vda - r)
            S += q ** (vda + r) * (-qs) ** k * c
    return S


def delO(r, vb, vc, ve, vda):
    """This is 1/(log q) times the derivative of the orbital for 1_(<= r)"""
    assert r >= 0
    assert vb + vc >= 0 and vb % 2 != vc % 2
    varkappa = ve - vda - r
    N = min(ve, floor((vb + vc) // 2 + r), vda + r)
    j = var("j")
    S = sum(q**j * (floor((2 * ve + vb + vc + 1) / 2) + r - 2 * j), j, 0, N)

    if varkappa >= 0 and vb + vc > 2 * vda:
        if varkappa % 2 == 0:
            S += q ** (vda + r) * (-varkappa / 2)
        else:
            S += q ** (vda + r) * (varkappa / 2 - (ve + (vb + vc) / 2 - 2 * vda - r))
    return (-1) ** (r + vc) * S


def delO_combo(r, vb, vc, ve, vda):
    """This is 1/(log q) times the derivative of the orbital for 1_(<= r) + 1_(<= r-1)"""

    assert r >= 0

    N = min(ve, (vb + vc) // 2 + r, vda + r)
    j = var("j")
    S = sum(q**j, j, 0, N)

    if ve >= vda + r and vb + vc > 2 * vda:
        assert N > 0
        if (r + ve + vda) % 2 == 1:
            C = (ve - r - vda) // 2
        else:
            C = (ve + vb + vc - r - 3 * vda) // 2
        S += C * q**N + (C + 1) * q ** (N - 1)

    elif 2 * vda > vb + vc and ve >= (vb + vc) // 2 + r:
        C = ve - r - (vb + vc) // 2
        S += C * q**N

    else:
        pass

    return (-1) ** (r + vc) * S


# Brute force auxiliary functions for the semi-Lie orbital
def O_brute_odd(r, vb, vc, ve, vda):
    """This is Case 5 in the situation where theta is odd"""
    assert vb % 2 != vc % 2
    assert vda >= 0
    assert r >= 0
    assert vda * 2 > vb + vc

    voffset = vda - vc
    S = 0

    for n2 in irange(0, ve):
        for n1 in irange(n2 - vb - r, n2 + vc + r):
            rho1 = max(-n1, -r - vc)
            rho2 = ceil((n2 - n1 - vc - r) / 2)
            if voffset < min(rho1, rho2):
                continue
            else:
                S += (
                    q ** (n2 - n1)
                    * q ** (-max(rho1, rho2))
                    * (-1) ** (n1 + n2)
                    * qs ** (n1 + n2)
                )
    return S


def O_brute_even(r, vb, vc, ve, vda):
    """This is Case 5/6+/6- in the situation where theta is odd"""
    assert vb % 2 != vc % 2
    assert vda >= 0
    assert r >= 0
    assert vda * 2 < vb + vc
    theta = vda * 2
    vt = theta // 2 - vc

    S = 0

    for n2 in irange(0, ve):
        # Case 5
        for m in irange(0, theta + 2 * r):
            n1 = n2 + vc + r - m
            S += (
                q ** (m - max(m - n2, ceil(m / 2)))
                * (-1) ** (n1 + n2)
                * qs ** (n1 + n2)
            )

        # Case 6+/6-
        center_plus = vb - theta // 2
        center_minus = theta // 2 - vc

        for m in irange(
            theta + 2 * r + 1, max(n2 + vc + r, 2 * vc + 2 * r + vt) + center_plus
        ):
            n1 = n2 + vc + r - m
            rho1 = max(m - n2, 0) - vc - r
            rho2 = m - 2 * vc - 2 * r - vt
            if center_plus >= min(rho1, rho2):
                S += (
                    q ** (n2 - n1)
                    * q ** (-max(rho1, rho2))
                    * (-1) ** (n1 + n2)
                    * qs ** (n1 + n2)
                )
            if center_minus >= min(rho1, rho2):
                S += (
                    q ** (n2 - n1)
                    * q ** (-max(rho1, rho2))
                    * (-1) ** (n1 + n2)
                    * qs ** (n1 + n2)
                )

    return S


## Formulas for the group AFL on S3(F)
def ARCH(a0, a1, w1, w2, k):
    assert a0 <= a1
    assert w1 + w2 <= (a1 - a0) / 2
    assert a0 <= k <= a1

    if a0 <= k <= a0 + w1:
        return k - a0
    elif a0 + w1 <= k <= a0 + w1 + w2:
        return w1 + floor((k - (a0 + w1)) / 2)
    elif a0 + w1 + w2 <= k <= a1 - (w1 + w2):
        return w1 + floor(w2 / 2)
    elif a1 - (w1 + w2) <= k <= a1 - w1:
        return w1 + floor(((a1 - w1) - k) / 2)
    elif a1 - w1 <= k <= a1:
        return a1 - k
    else:
        raise ValueError


def ARCH_sum_n(a0, a1, w1, w2):
    j = var("j")
    S = 0
    for k in irange(a0, a1):
        S += sum(q**j, j, 0, ARCH(a0, a1, w1, w2, k)) * (-qs) ** k
    return S


def ARCH_sum_c(a0, a1, w1, w2):
    S = 0
    for k in irange(a0, a1):
        S += ARCH(a0, a1, w1, w2, k) * (-qs) ** k
    return S


def O_group(r, l, delta, lam):
    if l % 2 == 1:
        assert lam == l
    else:
        assert l < lam
    if l < 0 or delta < 0:
        assert l == delta and l % 2 == 0
    else:
        assert l <= 2 * delta

    S = 0
    if l % 2 == 1:
        return ARCH_sum_n(-2 * r, 2 * delta + l + 2 * r, r, l)
    elif l % 2 == 0 and l >= 0:
        return ARCH_sum_n(-2 * r, 2 * delta + lam + 2 * r, r, l) + q ** (
            l / 2 + r
        ) * ARCH_sum_c(
            l - r, 2 * delta + lam - l + r, delta - l / 2, min(lam - l - 1, 2 * r)
        )
    else:
        vd = l // 2
        assert (avd := -vd) > 0  # avd := abs(vd)
        if r < avd:
            return 0
        return ARCH_sum_n(-2 * r, lam + 2 * (r - 2 * avd), r - avd, 0) + q ** (
            r - avd
        ) * ARCH_sum_c(-r - avd, lam + r - 3 * avd, 0, min(lam - 1, 2 * (r - avd)))


# Brute force auxiliary functions for the inhomogeneous group AFL orbital
def O_zero(r, l, delta):
    j = var("j")
    return sum(qs ** (2 * j), j, -r, delta + r)


def volDiskSingle(n, vxx, rho):
    assert n >= 1 and n >= rho
    if vxx < rho:
        return 0
    elif rho <= 0:
        return q ** (-n) * (1 - q ** (-2))
    else:
        return q ** (-(n + rho)) * (1 - q ** (-1))


def volDiskDouble(n, vxx1, vx12, rho1, rho2):
    assert rho1 >= rho2
    assert n >= 1 and n >= rho1

    if vxx1 >= rho1 and vx12 >= rho2:
        return factor(
            q ** (-(n + rho1)) * (1 - q ** (-1))
            if rho1 >= 1
            else q ** (-n) * (1 - q ** (-2))
        )
    else:
        return 0


def qs_weight(n, m):
    kappa = 1 / ((1 - q ** (-1)) * (1 - q ** (-2)))
    return (
        kappa
        * (-1) ** n
        * qs ** (2 * m - n)
        * q ** (2 * n - 2 * m)
        * q ** (2 * m)
        * (1 - q ** (-2))
    )


def O_case_1_2_brute(r, l, delta, lam=None):
    assert 0 < l < 2 * delta
    assert r >= 0
    if lam is None:
        assert l % 2 == 1
        lam = l
    else:
        assert l % 2 == 0
    INFINITY = 3 * (abs(r) + abs(l) + abs(delta) + 5)

    S = 0
    for n in irange(1, l + r):
        for m in irange(n - r, n + delta + r):
            r_n = ceil((n - r) / 2)
            r_m = m - delta - r
            if r_n >= r_m:
                S += volDiskDouble(
                    n, vxx1=min(l, delta), vx12=ceil(l / 2), rho1=r_n, rho2=r_m
                ) * qs_weight(n, m)
            else:
                S += volDiskDouble(
                    n, vxx1=lam, vx12=ceil(l / 2), rho1=r_m, rho2=r_n
                ) * qs_weight(n, m)
    return S


def O_case_3_4_brute(r, l, delta, lam):
    assert 0 <= l < 2 * delta
    assert r >= 0
    assert l % 2 == 0
    assert lam % 2 == 1
    INFINITY = 3 * (abs(r) + abs(l) + abs(delta) + abs(lam) + 5)

    S = 0
    for n in irange(l + r + 1, INFINITY):
        for m in irange(n - r, n + delta + r):
            r_n = n - l / 2 - r
            r_m = m - delta - r

            # Case 3+ and 4+
            if r_n > r_m:  # Case 3+
                S += volDiskDouble(
                    n,
                    vxx1=lam + delta - l,
                    vx12=lam - l / 2,
                    rho1=r_n,
                    rho2=r_m,
                ) * qs_weight(n, m)
            else:  # Case 4+
                S += volDiskDouble(
                    n,
                    vxx1=lam,
                    vx12=lam - l / 2,
                    rho1=r_m,
                    rho2=r_n,
                ) * qs_weight(n, m)

            # Cases 3- and 4-
            if r_n > r_m:  # Case 3-
                S += volDiskDouble(
                    n,
                    vxx1=delta,
                    vx12=l / 2,
                    rho1=r_n,
                    rho2=r_m,
                ) * qs_weight(n, m)
            else:  # Case 4-
                assert (
                    volDiskDouble(
                        n,
                        vxx1=lam,
                        vx12=l / 2,
                        rho1=r_m,
                        rho2=r_n,
                    )
                    == 0
                )

    return S


def O_ell_odd_brute(r, l, delta):
    return O_zero(r, l, delta) + O_case_1_2_brute(r, l, delta)


def O_ell_even_brute(r, l, delta, lam):
    return (
        O_zero(r, l, delta)
        + O_case_1_2_brute(r, l, delta, lam)
        + O_case_3_4_brute(r, l, delta, lam)
    )


def O_ell_neg_brute(r, vb, lam):
    assert r >= 0
    assert vb < 0
    assert lam % 2 == 1
    l = 2 * vb
    delta = 2 * vb
    INFINITY = 3 * (abs(r) + abs(l) + abs(delta) + abs(lam) + 5)

    S = 0
    for n in irange(1, INFINITY):
        for m in irange(n - r, n + delta + r):
            if n <= l + r:
                S += volDiskSingle(n, vxx=lam, rho=m - delta - r) * qs_weight(n, m)
            elif n > l + r:
                rho1 = max(n - l / 2 - r, m - delta - r)
                rho2 = min(n - l / 2 - r, m - delta - r)
                S += volDiskDouble(
                    n, vxx1=lam, vx12=lam, rho1=rho1, rho2=rho2
                ) * qs_weight(n, m)

    return S + O_zero(r, l, delta)


class SemiLieTest(unittest.TestCase):
    def get_params(self, r_min=0, r_max=20):
        params = {
            "r": randint(r_min, r_max),
            "vb": randint(-10, 20),
            "vda": randint(0, 20),
            "ve": randint(0, 20),
        }
        params["vc"] = randrange(1 - params["vb"], 21, 2)
        assert params["vb"] % 2 != params["vc"] % 2
        assert params["vb"] + params["vc"] >= 0
        return params

    def test_O_correct(self):
        params = self.get_params()
        brute_res = (
            O_brute_odd(**params)
            if params["vb"] + params["vc"] < 2 * params["vda"]
            else O_brute_even(**params)
        )

        orb = O(**params)
        self.assertEqual(orb.subs(q_s=1), 0)
        self.assertEqual(brute_res.subs(q=17, q_s=1337), orb.subs(q=17, q_s=1337))

    def test_delO_correct(self):
        params = self.get_params()
        self.assertEqual(
            delO(**params),
            derivative(O(**params), qs).subs(q_s=1),
        )

    def test_delO_combo_correct(self):
        params = self.get_params(r_min=1)
        r = params.pop("r")
        self.assertEqual(
            delO(r, **params) + delO(r - 1, **params),
            delO_combo(r, **params),
        )

    def test_matrix_for_full_rank_correct(self):
        N = 7
        params = self.get_params(r_min=0, r_max=N)
        del params["ve"]
        del params["r"]
        vb, vc, vda = params["vb"], params["vc"], params["vda"]
        theta = min(vb + vc, 2 * vda)

        # Here M0 = M, M1 = M', M2 = M''
        M0 = matrix(
            [
                vector(
                    (-1) ** (r + vc) * delO(r, vb, vc, i, vda) for r in range(0, n + 1)
                )
                for i in irange(0, n + theta // 2 + 1)
            ]
        )
        M1 = matrix([M0[0]] + [M0[i + 1] - M0[i] for i in range(M0.nrows() - 1)])
        M2 = matrix([M1[0], M1[1]] + [M1[i + 2] - M1[i] for i in range(M1.nrows() - 2)])

        for r in range(0, N + 1):
            t = r + theta // 2
            for i in irange(t + 2, theta // 2 + N + 1):
                self.assertEqual(M2[i][r], 0)
            C = (vb + vc - 1 - 2 * vda) / 2
            if theta % 2 == 1:
                self.assertEqual(M2[t + 1][r], q**t - (q ** (t - 1) if t > 0 else 0))
            elif t > 0:
                self.assertEqual(
                    M2[t + 1][r], -C * q**t - (C + 1) * (q ** (t - 1) if t > 0 else 0)
                )

    def test_kernel_large_r(self):
        params = self.get_params(r_min=2)
        r = max(params.pop("r"), params["ve"] + 2)
        self.assertEqual(
            delO(r, **params) + 2 * delO(r - 1, **params) + delO(r - 2, **params), 0
        )

    def test_kernel_full(self):
        params = self.get_params(r_min=5)
        vb, vc, vda, ve = params["vb"], params["vc"], params["vda"], params["ve"]
        r = params.pop("r")

        def delPhi(r):
            return (
                delO(r, **params)
                + delO(r - 1, **params)
                - q**2 * delO(r - 2, **params)
                - q**2 * delO(r - 3, **params)
            )

        left_bound = ve - min((vb + vc - 1) / 2, vda) + 2
        right_bound = left_bound + 2
        expr = expand(delPhi(r) + (q + 1) * delPhi(r - 1) + q * delPhi(r - 2))
        if not left_bound <= r <= right_bound:
            self.assertEqual(expr, 0)
        # in fact I think the following is the exact criteria for it to be nonzero
        # although I don't claim this in the paper
        # but we'll put it in the test case just for kicks
        elif left_bound <= r <= right_bound - 1:
            self.assertNotEqual(expr, 0)
        elif r == right_bound - 1 and vb + vc < 2 * vda:
            self.assertNotEqual(expr, 0)
        elif r == right_bound - 1 and vb + vc >= 2 * vda:
            self.assertEqual(expr, 0)


class InhomogeneousGroupTest(unittest.TestCase):
    def get_params(self):
        l = randint(-5, 20)
        if l < 0:
            l *= 2

        params = {
            "r": randint(0, 20),
            "l": l,
            "lam": l if l % 2 == 1 else randrange(max(0, l + l % 2) + 1, 22, 2),
            "delta": l if l < 0 else randint(l // 2, 20),
        }
        return params

    def test_O_correct(self):
        params = self.get_params()
        l = params["l"]
        orb = O_group(**params)
        if l < 0:
            brute_res = O_ell_neg_brute(r=params["r"], vb=l // 2, lam=params["lam"])
        elif l % 2 == 0:
            brute_res = O_ell_even_brute(**params)
        elif l % 2 == 1:
            del params["lam"]
            brute_res = O_ell_odd_brute(**params)

        self.assertEqual(brute_res.subs(q=17, q_s=1337), orb.subs(q=17, q_s=1337))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "checkthesis",
        description="Checks the formulas in Evan's thesis for consistency",
    )
    parser.add_argument(
        "--trials", default=5, type=int, help="Number of trials to run."
    )

    parser.add_argument(
        "--test",
        default="",
        type=str,
        help="The name of a specific test to run (if empty, runs all).",
    )

    args = parser.parse_args()
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    for _ in range(args.trials):
        if args.test:
            suite.addTest(SemiLieTest(args.test))
            suite.addTest(InhomogeneousGroupTest(args.test))
        else:
            suite.addTest(loader.loadTestsFromTestCase(SemiLieTest))
            suite.addTest(loader.loadTestsFromTestCase(InhomogeneousGroupTest))
    runner = unittest.TextTestRunner()
    runner.run(suite)
