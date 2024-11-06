# Setup

q = var("q")
qs = var("q_s")  # = q**s


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


def delO_combo(r, vb, vc, ve, vda):
    """This is 1/(log q) times the derivative of the orbital integral
    for 1_(<= r) + 1_(<= r-1)"""

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


def verify_O_correct(num_trials=1):
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
                    q ** (m - max(m - n2, ceil(m / 2), 0))
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

    for test_case_number in range(0, num_trials):
        r = randint(0, 40)
        vb = randint(-20, 40)
        vc = randrange(1 - vb, 40, 2)
        vda = randint(0, 40)
        ve = randint(0, 40)

        predicted = O(r, vb, vc, ve, vda)
        assert vb % 2 != vc % 2
        assert vb + vc >= 0
        if vb + vc < 2 * vda:
            brute_force_result = O_brute_odd(r, vb, vc, ve, vda)
        else:
            brute_force_result = O_brute_even(r, vb, vc, ve, vda)

        assert brute_force_result.subs(q=17, q_s=1337) == predicted.subs(q=17, q_s=1337)
    return True


def verify_delO_combo_correct(num_trials=1):
    for test_case_number in range(0, num_trials):
        r = randint(1, 40)
        vb = randint(-20, 40)
        vc = randrange(1 - vb, 40, 2)
        vda = randint(0, 40)
        ve = randint(0, 40)
        actual = derivative(O(r, vb, vc, ve, vda) + O(r - 1, vb, vc, ve, vda), qs).subs(
            q_s=1
        )
        assert delO_combo(r, vb, vc, ve, vda) == actual
    return True


if __name__ == "__main__":
    print(verify_O_correct(num_trials=10000))
    print(verify_delO_combo_correct(num_trials=10000))
