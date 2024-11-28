# Setup

import argparse
import unittest

q = var("q")
qs = var("q_s")  # = q^s


def irange(start, stop):
    return range(start, stop + 1)


def O(r, vb, vc, ve, vda):
    """Value of the overall orbital integral"""
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


def M0(n, vb, vc, vda):
    """This is the matrix M used to prove the trivial kernel result"""
    theta = min(vb + vc, 2 * vda)
    return matrix(
        [
            vector((-1) ** (r + vc) * delO(r, vb, vc, i, vda) for r in range(0, n + 1))
            for i in irange(0, n + theta // 2 + 1)
        ]
    )


def M1(n, vb, vc, vda):
    M = M0(n, vb, vc, vda)
    return matrix([M[0]] + [M[i + 1] - M[i] for i in range(M.nrows() - 1)])


def M2(n, vb, vc, vda):
    M = M1(n, vb, vc, vda)
    return matrix([M[0], M[1]] + [M[i + 2] - M[i] for i in range(M.nrows() - 2)])


class SemiLieTest(unittest.TestCase):
    def get_params(self, r_min=0, r_max=20):
        params = {
            "r": randint(r_min, r_max),
            "vb": randint(-10, 20),
            "vda": randint(0, 20),
            "ve": randint(0, 20),
        }
        params["vc"] = randrange(1 - params["vb"], 20, 2)
        assert params["vb"] % 2 != params["vc"] % 2
        assert params["vb"] + params["vc"] >= 0
        return params

    def test_O_correct(self):
        params = self.get_params()
        brute_force_result = (
            O_brute_odd(**params)
            if params["vb"] + params["vc"] < 2 * params["vda"]
            else O_brute_even(**params)
        )

        self.assertEqual(
            brute_force_result.subs(q=17, q_s=1337), O(**params).subs(q=17, q_s=1337)
        )

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

    def test_matrix_correct(self):
        N = 7
        params = self.get_params(r_min=0, r_max=N)
        del params["ve"]
        del params["r"]
        vb, vc, vda = params["vb"], params["vc"], params["vda"]
        theta = min(vb + vc, 2 * vda)
        M = M2(n=N, **params)
        for r in range(0, N + 1):
            t = r + theta // 2
            for i in irange(t + 2, theta // 2 + N + 1):
                self.assertEqual(M[i][r], 0)
            C = (vb + vc - 1 - 2 * vda) / 2
            if theta % 2 == 1:
                self.assertEqual(M[t + 1][r], q**t - (q ** (t - 1) if t > 0 else 0))
            elif t > 0:
                self.assertEqual(
                    M[t + 1][r], -C * q**t - (C + 1) * (q ** (t - 1) if t > 0 else 0)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "thesischeck",
        description="Checks the formulas in Evan's thesis for the semi-Lie case",
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
        else:
            suite.addTest(loader.loadTestsFromTestCase(SemiLieTest))
    runner = unittest.TextTestRunner()
    runner.run(suite)
