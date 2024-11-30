import argparse
import unittest

q = var("q")
qs = var("q_s")  # = q^s


def irange(start, stop):
    return range(start, stop + 1)


def print_coeffs(expression) -> None:
    """
    If you're working in a Juptyer notebook, and you have a polynomial in q and qs,
    you can use this utility function to print out the coefficients of q^s each on their
    own line.

    :param expression: The polynomial to eb pretty-printed
    """
    if expand(expression) == 0:
        show(0)
    else:
        for c in expand(expression).coefficients(qs, sparse=True):
            show(c[1], "." * 10, c[0])


# Semi-Lie Orbital and its derivatives
def O(r, vb, vc, ve, vda):
    assert vb % 2 != vc % 2, (vb, vc)
    assert r >= 0, r
    assert vb + vc >= 0, (vb, vc)
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
    assert r >= 0, r
    assert vb + vc >= 0 and vb % 2 != vc % 2, (vb, vc)
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
    assert r >= 1
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


# Formulas for the group AFL on S3(F)
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


# Orbital for S3
def O_for_S3(r, l, delta, lam):
    if l % 2 == 1:
        assert lam == l, (l, lam)
    else:
        assert l < lam, (l, lam)
    if l < 0 or delta < 0:
        assert l == delta and l % 2 == 0, (l, delta)
    else:
        assert l <= 2 * delta, (l, delta)

    S = 0
    if l % 2 == 1:
        return ARCH_sum_n(-2 * r, 2 * delta + l + 2 * r, r, l)
    elif l % 2 == 0 and l >= 0:
        n_sum = ARCH_sum_n(-2 * r, 2 * delta + lam + 2 * r, r, l)
        c_sum = ARCH_sum_c(
            l - r, 2 * delta + lam - l + r, delta - l / 2, min(lam - l - 1, 2 * r)
        )
        return n_sum + q ** (r + l / 2) * c_sum
    else:
        assert (avd := -l / 2) > 0  # avd := abs(vd)
        if r < avd:
            return 0
        n_sum = ARCH_sum_n(-2 * r, lam + 2 * (r - 2 * avd), r - avd, 0)
        c_sum = ARCH_sum_c(-r - avd, lam + r - 3 * avd, 0, min(lam - 1, 2 * (r - avd)))
        return n_sum + q ** (r - avd) * c_sum


# Derivatives of the orbital integral
def ARCH_deriv_n(r, C, W, H):
    j = var("j")
    assert W > 4 * H >= 0, (W, H)
    assert W % 2 == 1, W
    return -(
        (-1) ** (r + C) * sum(((W + 1) / 2 + r - 2 * (j - r)) * q**j, j, r + 1, r + H)
        + sum((-1) ** (j + C) * ((W + 1) / 2 + 2 * r - j) * q**j, j, 0, r)
    )


def ARCH_deriv_c(r, C, W, L):
    k = var("k")
    assert L >= 1, L
    assert W >= 0, W
    assert L % 2 == 1, L
    return (-1) ** (r + W + C) * (
        W / 2 - (L - 1) / 2 * r if W % 2 == 0 else -(W + L) / 2 - (L + 1) / 2 * r
    )


def delO_for_S3_via_arch(r, l, delta, lam):
    if l % 2 == 1:
        assert lam == l, (l, lam)
    else:
        assert l < lam, (l, lam)
    if l < 0 or delta < 0:
        assert l == delta and l % 2 == 0, (l, delta)
    else:
        assert l <= 2 * delta, (l, delta)

    if l % 2 == 1:
        return ARCH_deriv_n(r, C=0, W=l + 2 * delta, H=(l - 1) / 2)
    elif l >= 0:
        n_sum = ARCH_deriv_n(r, C=0, W=lam + 2 * delta, H=l / 2)
        c_sum = ARCH_deriv_c(r, C=0, W=delta - l / 2, L=lam - l)
        return n_sum + q ** (r + l / 2) * c_sum
    else:
        assert (avd := -l / 2) > 0  # avd := abs(vd)
        if r < avd:
            return 0
        n_sum = ARCH_deriv_n(r - avd, C=-2 * avd, W=lam, H=0)
        c_sum = ARCH_deriv_c(r - avd, C=0, W=0, L=lam)
        return n_sum + q ** (r - avd) * c_sum


def delO_for_S3(r, l, delta, lam):
    if l % 2 == 1:
        assert lam == l, (l, lam)
    else:
        assert l < lam, (l, lam)
    if l < 0 or delta < 0:
        assert l == delta and l % 2 == 0, (l, delta)
    else:
        assert l <= 2 * delta, (l, delta)

    j = var("j")
    if l % 2 == 1:
        S = (-1) ** (r + 1) * sum(
            ((l + 2 * delta + 1) / 2 + 3 * r - 2 * j) * q**j, j, r + 1, r + (l - 1) / 2
        )
        S += sum(
            (-1) ** (j + 1) * ((l + 2 * delta + 1) / 2 + 2 * r - j) * q**j, j, 0, r
        )
        return S
    elif l >= 0:
        frac = (lam + 2 * delta + 1) / 2
        S = (-1) ** (r + 1) * sum((frac + 3 * r - 2 * j) * q**j, j, r + 1, r + l / 2)
        S += sum((-1) ** (j + 1) * (frac + 2 * r - j) * q**j, j, 0, r)
        S += (
            (-1) ** (r + delta - l / 2)
            * q ** (r + l / 2)
            * (
                (delta - l / 2) / 2 - (lam - l - 1) / 2 * r
                if (delta - l // 2) % 2 == 0
                else -(delta - 3 * l / 2 + lam) / 2 - (lam - l + 1) / 2 * r
            )
        )
        return S
    else:
        assert (avd := -l / 2) > 0  # avd := abs(vd)
        t = r - avd  # top-level exponent (if t < 0, then S = 0)
        S = sum((-1) ** (j + 1) * ((lam + 1) / 2 + 2 * t - j) * q**j, j, 0, t)
        S += (-1) ** (t + 1) * (lam - 1) / 2 * max(t, 0) * q**t
        return S


# Geometric side --- Gross-Keating and friends
def gross_keating_sum(n1, n2):
    j = var("j")
    assert 0 <= n1 <= n2
    if n1 % 2 == 1:
        return sum((n1 + n2 - 4 * j) * q**j, j, 0, (n1 - 1) / 2)
    else:
        return sum((n1 + n2 - 4 * j) * q**j, j, 0, n1 / 2 - 1) + (
            n2 - n1 + 1
        ) / 2 * q ^ (n1 / 2)


def GK(r, vb, vc, ve, vda):
    assert r >= 0
    assert vb + vc >= 0 and vb % 2 != vc % 2
    vxx = 2 * ve
    vxy = r + vda + ve
    vdet = 2 * r + 2 * ve + vb + vc
    vyy = min(2 * vxy, vdet) - vxx

    return gross_keating_sum(min(vxx, vxy, vyy), vdet - min(vxx, vxy, vyy))


def clean_intersection(r, vb, vc, ve, vda):
    vbeta = (vb + vc - 1) / 2
    N = min(ve, vbeta + r, vda + r)
    C = vbeta - vda
    if ve - r == vda <= vbeta:
        return (C + 1) * q**N + (C + 2) * q ** (N - 1)
    elif vbeta + r < min(ve, vda + r):
        return 2 * q**N
    else:
        return q**N + q ** (N - 1)


class RandThesisTest(unittest.TestCase):
    def get_semi_lie_params(self, r_min=0, r_max=15):
        params = {
            "r": randint(r_min, r_max),
            "vb": randint(-5, 15),
            "vda": randint(0, 15),
            "ve": randint(0, 15),
        }
        params["vc"] = randrange(1 - params["vb"], 15, 2)
        assert params["vb"] % 2 != params["vc"] % 2
        assert params["vb"] + params["vc"] >= 0
        return params

    def test_O(self):
        params = self.get_semi_lie_params()

        def O_brute_odd(r, vb, vc, ve, vda):
            assert vb % 2 != vc % 2, (vb, vc)
            assert vda >= 0, vda
            assert r >= 0, r
            assert vda * 2 > vb + vc, (vda, vb, vc)

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
            assert vb % 2 != vc % 2, (vb, vc)
            assert vda >= 0, vda
            assert r >= 0, r
            assert vda * 2 < vb + vc, (vda, vb, vc)
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
                    theta + 2 * r + 1,
                    max(n2 + vc + r, 2 * vc + 2 * r + vt) + center_plus,
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

        brute_res = (
            O_brute_odd(**params)
            if params["vb"] + params["vc"] < 2 * params["vda"]
            else O_brute_even(**params)
        )

        orb = O(**params)
        self.assertEqual(orb.subs(q_s=1), 0)
        self.assertEqual(brute_res.subs(q=17, q_s=1337), orb.subs(q=17, q_s=1337))

    def test_delO(self):
        params = self.get_semi_lie_params()
        self.assertEqual(
            delO(**params),
            derivative(O(**params), qs).subs(q_s=1),
        )

    def test_delO_combo(self):
        params = self.get_semi_lie_params(r_min=1)
        r = params.pop("r")
        self.assertEqual(
            delO(r, **params) + delO(r - 1, **params),
            delO_combo(r, **params),
        )

    def test_matrix_upper_triangular(self):
        N = 7
        params = self.get_semi_lie_params(r_min=0, r_max=N)
        del params["ve"]
        del params["r"]
        vb, vc, vda = params["vb"], params["vc"], params["vda"]
        theta = min(vb + vc, 2 * vda)

        # Here M0 = M, M1 = M', M2 = M''
        M0 = matrix(
            [
                vector(
                    (-1) ** (r + vc) * delO(r, vb, vc, i, vda) for r in range(0, N + 1)
                )
                for i in irange(0, N + theta // 2 + 1)
            ]
        )
        M1 = matrix([M0[0]] + [M0[i + 1] - M0[i] for i in range(M0.nrows() - 1)])
        M2 = matrix([M1[0], M1[1]] + [M1[i + 2] - M1[i] for i in range(M1.nrows() - 2)])

        for r in range(0, N + 1):
            t = r + theta // 2
            for i in irange(t + 2, N + theta // 2 + 1):
                self.assertEqual(M2[i][r], 0)
            C = (vb + vc - 1 - 2 * vda) / 2
            if theta % 2 == 1:
                self.assertEqual(M2[t + 1][r], q**t - (q ** (t - 1) if t > 0 else 0))
            else:
                self.assertEqual(
                    M2[t + 1][r], -C * q**t - (C + 1) * (q ** (t - 1) if t > 0 else 0)
                )

    def test_kernel_large_r(self):
        params = self.get_semi_lie_params(r_min=2)
        r = max(params.pop("r"), params["ve"] + 2)
        self.assertEqual(
            delO(r, **params) + 2 * delO(r - 1, **params) + delO(r - 2, **params), 0
        )

    def test_kernel_full(self):
        params = self.get_semi_lie_params(r_min=5)
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

    def get_S3_params(self, r_min=0, r_max=10):
        l = randint(-5, 10)
        if l < 0:
            l *= 2

        params = {
            "r": randint(r_min, r_max),
            "l": l,
            "lam": l if l % 2 == 1 else randrange(max(0, l + l % 2) + 1, 17, 2),
            "delta": l if l < 0 else randint((l + 1) // 2, 10),
        }
        return params

    def test_O_for_S3(self):
        # Brute force auxiliary functions for the inhomogeneous group AFL orbital
        def O_zero(r, l, delta):
            j = var("j")
            return sum(qs ** (2 * j), j, -r, delta + r)

        def vol_disk_single(n, vxx, rho):
            assert n >= 1 and n >= rho
            if vxx < rho:
                return 0
            elif rho <= 0:
                return q ** (-n) * (1 - q ** (-2))
            else:
                return q ** (-(n + rho)) * (1 - q ** (-1))

        def vol_disk_double(n, vxx1, vx12, rho1, rho2):
            assert rho1 >= rho2
            assert n >= 1 and n >= rho1

            if vxx1 >= rho1 and vx12 >= rho2:
                return (
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
            assert 0 <= l <= 2 * delta, (l, delta)
            assert r >= 0, r
            if lam is None:
                assert l % 2 == 1
                lam = l
            else:
                assert l % 2 == 0

            S = 0
            for n in irange(1, l + r):
                for m in irange(n - r, n + delta + r):
                    r_n = ceil((n - r) / 2)
                    r_m = m - delta - r
                    if r_n >= r_m:
                        S += vol_disk_double(
                            n, vxx1=min(l, delta), vx12=ceil(l / 2), rho1=r_n, rho2=r_m
                        ) * qs_weight(n, m)
                    else:
                        S += vol_disk_double(
                            n, vxx1=lam, vx12=ceil(l / 2), rho1=r_m, rho2=r_n
                        ) * qs_weight(n, m)
            return S

        def O_ell_odd_brute(r, l, delta):
            return O_zero(r, l, delta) + O_case_1_2_brute(r, l, delta)

        def O_case_3_4_brute(r, l, delta, lam):
            assert 0 <= l <= 2 * delta, (l, delta)
            assert r >= 0, r
            assert l % 2 == 0, l
            assert lam % 2 == 1, lam
            INFINITY = abs(r) + abs(l) + abs(delta) + abs(lam) + 1

            S = 0
            for n in irange(l + r + 1, INFINITY):
                for m in irange(n - r, n + delta + r):
                    r_n = n - l / 2 - r
                    r_m = m - delta - r

                    # Case 3+ and 4+
                    if r_n > r_m:  # Case 3+
                        S += vol_disk_double(
                            n,
                            vxx1=lam + delta - l,
                            vx12=lam - l / 2,
                            rho1=r_n,
                            rho2=r_m,
                        ) * qs_weight(n, m)
                    else:  # Case 4+
                        S += vol_disk_double(
                            n,
                            vxx1=lam,
                            vx12=lam - l / 2,
                            rho1=r_m,
                            rho2=r_n,
                        ) * qs_weight(n, m)

                    # Cases 3- and 4-
                    if r_n > r_m:  # Case 3-
                        S += vol_disk_double(
                            n,
                            vxx1=delta,
                            vx12=l / 2,
                            rho1=r_n,
                            rho2=r_m,
                        ) * qs_weight(n, m)
                    else:  # Case 4-
                        assert (
                            vol_disk_double(
                                n,
                                vxx1=lam,
                                vx12=l / 2,
                                rho1=r_m,
                                rho2=r_n,
                            )
                            == 0
                        )

            return S

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
            INFINITY = abs(r) + abs(l) + abs(delta) + abs(lam) + 1

            S = 0
            for n in irange(1, INFINITY):
                for m in irange(n - r, n + delta + r):
                    if n <= l + r:
                        S += vol_disk_single(n, vxx=lam, rho=m - delta - r) * qs_weight(
                            n, m
                        )
                    elif n > l + r:
                        rho1 = max(n - l / 2 - r, m - delta - r)
                        rho2 = min(n - l / 2 - r, m - delta - r)
                        S += vol_disk_double(
                            n, vxx1=lam, vx12=lam, rho1=rho1, rho2=rho2
                        ) * qs_weight(n, m)

            return S + O_zero(r, l, delta)

        params = self.get_S3_params()
        l = params["l"]
        orb = O_for_S3(**params)
        if l < 0:
            brute_res = O_ell_neg_brute(r=params["r"], vb=l // 2, lam=params["lam"])
        elif l % 2 == 0:
            brute_res = O_ell_even_brute(**params)
        elif l % 2 == 1:
            del params["lam"]
            brute_res = O_ell_odd_brute(**params)

        self.assertEqual(orb.subs(q_s=1), 0)
        self.assertEqual(brute_res.subs(q=17, q_s=1337), orb.subs(q=17, q_s=1337))

    def test_delO_for_S3_via_arch(self):
        params = self.get_S3_params()
        self.assertEqual(
            derivative(O_for_S3(**params), qs).subs(q_s=1),
            delO_for_S3_via_arch(**params),
        )

    def test_delO_for_S3(self):
        params = self.get_S3_params()
        self.assertEqual(delO_for_S3_via_arch(**params), delO_for_S3(**params))

    def test_ker_for_S3(self):
        params = self.get_S3_params(r_min=3)
        l = params["l"]
        r = params.pop("r")
        deriv = (
            (delO_for_S3(r, **params) - delO_for_S3(r - 1, **params))
            + 2 * q * (delO_for_S3(r - 1, **params) - delO_for_S3(r - 2, **params))
            + q**2 * (delO_for_S3(r - 2, **params) - delO_for_S3(r - 3, **params))
        )
        if l < 0:
            if r < -l // 2:
                self.assertEqual(deriv, 0)
            elif r >= -l // 2 + 3:
                self.assertEqual(deriv, -2 * q - 2)
        else:
            self.assertEqual(deriv, -2 * q - 2)

    def test_GK_to_orbital(self):
        params = self.get_semi_lie_params(r_min=1)
        omega = (-1) ** (params["r"] + params["vc"])  # transfer factor
        ve = params.pop("ve")
        self.assertEqual(
            omega * (delO(ve=ve, **params) + delO(ve=ve - 1, **params)),
            GK(ve=ve, **params),
        )

    def test_clean_intersection(self):
        params = self.get_semi_lie_params(r_min=1)
        if params["ve"] == 0:
            params["ve"] += 1
        r, vb, vc, ve, vda = (
            params["r"],
            params["vb"],
            params["vc"],
            params["ve"],
            params["vda"],
        )
        self.assertEqual(
            (GK(r, vb, vc, ve, vda) - GK(r - 1, vb, vc, ve, vda))
            - (GK(r, vb, vc, ve - 1, vda) - GK(r - 1, vb, vc, ve - 1, vda)),
            clean_intersection(**params),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "checkthesis",
        description="Checks the formulas in Evan's thesis for consistency",
    )
    parser.add_argument(
        "--trials", default=5, type=int, help="Number of trials to run."
    )
    parser.add_argument("--seed", type=int, help="Random seed passed to Sage")
    parser.add_argument(
        "--failfast", action="store_true", help="On failure, immediately stop."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--verbose", action="store_true", help="Set verbosity to 2.")
    group.add_argument("--quiet", action="store_true", help="Set verbosity to 0.")
    parser.add_argument(
        "--test",
        default="",
        type=str,
        help="The name of a specific test to run (if empty, runs all).",
    )
    args = parser.parse_args()
    if args.verbose is True:
        verbosity = 2
    elif args.quiet is True:
        verbosity = 0
    else:
        verbosity = 1

    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    if args.seed:
        set_random_seed(args.seed)
    print(f"Using random seed {initial_seed()}")
    for _ in range(args.trials):
        if args.test:
            suite.addTest(RandThesisTest(args.test))
        else:
            suite.addTest(loader.loadTestsFromTestCase(RandThesisTest))
    runner = unittest.TextTestRunner(failfast=args.failfast, verbosity=verbosity)
    runner.run(suite)
