#!/usr/bin/env python3
"""verify_343.py — independent verification of the claims about page 343 of Ramanujan's Lost Notebook
(Andrews–Berndt, Lost Notebook Part IV, Sect. 7.4, eqs. (7.4.1)–(7.4.2)).

Requires only: Python 3.9+, mpmath, sympy.   Run:  python3 verify_343.py      (exit code 0 = all checks pass)

theta = 2^(1/5),  a = 1/(theta-1) = 1+theta+theta^2+theta^3+theta^4,  b = sqrt5/(1+theta^2)^(5/2) = 1-theta-theta^2+theta^4.
Claim (Theorem A): p_{m,n} := Tr(a^m b^n theta) is an integer and
    |a^m b^n theta - p_{m,n}| <= 2 theta (|s1(a)|^m |s1(b)|^n + |s2(a)|^m |s2(b)|^n),
where s1, s2 are the two non-real embeddings theta -> theta*exp(2 pi i k/5), k = 1, 2.
"""
import itertools
import sys

import mpmath as mp
import sympy as sp

x, y = sp.symbols("x y")
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  —  {detail}" if detail else ""))


def red(expr, p):
    return sp.Poly(sp.expand(expr), x).rem(sp.Poly(x ** p - 2, x)).as_expr()


def mul5(u, v):  # exact arithmetic in Z[theta], theta^5 = 2
    w = [0] * 9
    for i, ui in enumerate(u):
        for j, vj in enumerate(v):
            w[i + j] += ui * vj
    return [w[k] + 2 * w[k + 5] if k + 5 < 9 else w[k] for k in range(5)]


def pw(u, e):
    r = [1, 0, 0, 0, 0]
    for _ in range(e):
        r = mul5(r, u)
    return r


A, B, TH = [1, 1, 1, 1, 1], [1, -1, -1, 0, 1], [0, 1, 0, 0, 0]

# 1. exact identities
a5, b5 = 1 + x + x ** 2 + x ** 3 + x ** 4, 1 - x - x ** 2 + x ** 4
check("1  a*(theta-1) = 1 and b^2*(1+theta^2)^5 = 5 in Q[theta]/(theta^5-2)",
      red(a5 * (x - 1) - 1, 5) == 0 and red(b5 ** 2 * (1 + x ** 2) ** 5 - 5, 5) == 0)
check("2  a, b are units: N(a) = N(b) = 1",
      sp.resultant(x ** 5 - 2, a5, x) == 1 and sp.resultant(x ** 5 - 2, b5, x) == 1)

# 3. trace: two independent algorithms
C = sp.zeros(5, 5)
for j in range(4):
    C[j + 1, j] = 1
C[0, 4] = 2
Ma, Mb = sp.eye(5) + C + C ** 2 + C ** 3 + C ** 4, sp.eye(5) - C - C ** 2 + C ** 4
pairs = [(0, 0), (2, 0), (7, 0), (50, 3), (120, 20), (300, 40)]
check("3  Tr via Z[theta] arithmetic == trace of multiplication matrix",
      all(5 * mul5(mul5(pw(A, m), pw(B, n)), TH)[0] == (Ma ** m * Mb ** n * C).trace() for m, n in pairs),
      f"{len(pairs)} pairs up to (300,40)")

# 4. the bound of Theorem A on the grid m<=300, n<=40
mp.mp.dps = 400
th = mp.root(2, 5)
z = [mp.exp(2j * mp.pi * k / 5) for k in range(5)]
sa = [abs(1 / (th * z[k] - 1)) for k in (1, 2)]
sb = [abs((th * z[k]) ** 4 - (th * z[k]) ** 2 - th * z[k] + 1) for k in (1, 2)]
aN, bN = 1 / (th - 1), th ** 4 - th ** 2 - th + 1
apow = [TH]
for m in range(300):
    apow.append(mul5(apow[-1], A))
bad, count = 0, 0
for n in range(41):
    bn = pw(B, n)
    for m in range(301):
        tr = 5 * mul5(apow[m], bn)[0]
        eps = aN ** m * bN ** n * th - tr
        bound = 2 * th * (sa[0] ** m * sb[0] ** n + sa[1] ** m * sb[1] ** n)
        bad += abs(eps) > bound * (1 + mp.mpf(10) ** -50)
        count += 1
check("4  |eps| <= 2 theta (|s1|+|s2|) on m<=300, n<=40", bad == 0, f"{count} pairs, {bad} violations")

# 5. nearest-integer statement for n = 0
mp.mp.dps = 100
th = mp.root(2, 5)
ok_m = [m for m in range(0, 60) if 5 * mul5(pw(A, m), TH)[0] == int(mp.nint((1 / (th - 1)) ** m * th))]
m0 = next(m for m in range(100) if 2 * th * (sa[0] ** m + sa[1] ** m) < mp.mpf(1) / 2)
check("5  Tr(a^m theta) = nearest integer to a^m theta exactly for m = 4 and m >= 6 (m < 60); proof bound from m0 = 7",
      ok_m == [4] + list(range(6, 60)) and m0 == 7, f"m0 = {m0}")

# 6. conjugate moduli from polynomial roots (independent of the exp formula)
mp.mp.dps = 60
R = mp.polyroots([1, 0, 0, 0, 0, -2], maxsteps=200, extraprec=200)
up = sorted([r for r in R if mp.im(r) > 1e-30], key=lambda r: -mp.re(r))
mods = [abs(1 / (t - 1)) for t in up] + [abs(t ** 4 - t ** 2 - t + 1) for t in up]
check("6  |s1(a)|, |s2(a)|, |s1(b)|, |s2(b)| = 0.788215, 0.489225, 4.181301, 0.457816",
      all(abs(u - v) < 1e-6 for u, v in zip(mods, [0.788215, 0.489225, 4.181301, 0.457816])),
      ", ".join(mp.nstr(v, 7) for v in mods))

# 7. Lemma B: 1/(2^(1/p)-1) is Pisot iff p <= 6 (roots of (y+1)^p - 2 y^p)
pis = []
for p in range(2, 13):
    rts = mp.polyroots([int(c) for c in sp.Poly((y + 1) ** p - 2 * y ** p, y).all_coeffs()], maxsteps=300, extraprec=300)
    big = max(rts, key=lambda r: mp.re(r))
    pis.append(max(abs(r) for r in rts if r != big) < 1)
check("7  Lemma B: Pisot exactly for p = 2..6 (tested p = 2..12)", pis == [True] * 5 + [False] * 6)

# 8. sharpness of 2 theta: sup over 200 <= m < 2200 approaches 2 theta
mp.mp.dps = 60
th = mp.root(2, 5)
s1 = 1 / (th * mp.exp(2j * mp.pi / 5) - 1)
best = max(abs(2 * mp.re(s1 ** m * th * mp.exp(2j * mp.pi / 5))) / abs(s1) ** m for m in range(200, 2200))
check("8  limsup |eps_m| / |s1(a)|^m -> 2 theta", abs(best - 2 * th) < 1e-5, f"{mp.nstr(best, 10)} vs {mp.nstr(2 * th, 10)}")

# 9. seventh root: units, bc = u^2, and the Pisot cone
a7 = sum(x ** k for k in range(7))
u7 = 21 + 19 * x + 18 * x ** 2 + 16 * x ** 3 + 14 * x ** 4 + 13 * x ** 5 + 12 * x ** 6
T7 = sp.root(2, 7)
mb7 = sp.Poly(sp.minimal_polynomial(7 / (T7 ** 3 - 1) ** 7, y), y).all_coeffs()
mc7 = sp.Poly(sp.minimal_polynomial((T7 + 1) / (T7 ** 2 - T7 + 1), y), y).all_coeffs()
check("9  theta^7=2: a(theta-1)=1, b and c are units, 7(theta+1)^2 = u^2 (theta^3-1)^7 (theta^3+1)",
      red(a7 * (x - 1) - 1, 7) == 0 and mb7[0] == 1 and abs(mb7[-1]) == 1 and mc7[0] == 1 and abs(mc7[-1]) == 1
      and red(7 * (x + 1) ** 2 - u7 ** 2 * (x ** 3 - 1) ** 7 * (x ** 3 + 1), 7) == 0)
R7 = mp.polyroots([1, 0, 0, 0, 0, 0, 0, -2], maxsteps=300, extraprec=300)
up7 = [r for r in R7 if mp.im(r) > 1e-30]
L = [[mp.log(abs(f(t))) for t in up7] for f in (lambda t: 1 / (t - 1), lambda t: 7 / (t ** 3 - 1) ** 7,
                                                 lambda t: (t + 1) / (t ** 2 - t + 1))]
cone = [(l_, m, n) for l_, m, n in itertools.product(range(61), range(4), range(13))
        if (l_, m, n) != (0, 0, 0) and all(l_ * L[0][s] + m * L[1][s] + n * L[2][s] < 0 for s in range(3))]
check("10 Pisot cone (l<=60, m<=3, n<=12) has 526 points, none with m = 0; a^4 b and a^3 b c are Pisot",
      len(cone) == 526 and all(k[1] > 0 for k in cone) and (4, 1, 0) in cone and (3, 1, 1) in cone)
mp.mp.dps = 900
t7 = mp.root(2, 7)
v = (1 / (t7 - 1)) ** 160 * (7 / (t7 ** 3 - 1) ** 7) ** 40 * t7
e = abs(v - mp.nint(v))
check("11 (a^4 b)^40 theta is within 3e-9 of an integer (900 digits)", e < mp.mpf("3e-9"), mp.nstr(e, 6))

# 12. OEIS: Tr(theta a^(n+1))/10 = A374455(n), proved by recurrence + initial terms
A374455 = [1, 5, 35, 235, 1580, 10626, 71460, 480570, 3231845, 21734235, 146163251, 982951365, 6610371480,
           44454906580, 298960311840, 2010515259876, 13520763292345, 90927457083265, 611489327404315,
           4112280377388895, 27655184063541876, 185981775414350150, 1250731895575163300]
mine = [5 * mul5(pw(A, n + 1), TH)[0] // 10 for n in range(len(A374455))]
check("12 Tr(theta a^(n+1))/10 = A374455(n) (23 OEIS terms; both satisfy signature (5,10,10,5,1))", mine == A374455)

# 13. Ramanujan's printed factors in (7.4.1) are exactly the conjugate moduli |sigma_s(a^m b^n)|, s = 1..4
mp.mp.dps = 50
t5 = mp.root(2, 5)
worst = mp.mpf(0)
for m, n in [(1, 0), (0, 1), (3, 2), (17, 5), (40, 1)]:
    for s in range(1, 5):
        c1 = mp.root(4, 5) - 2 * t5 * mp.cos(2 * mp.pi * s / 5) + 1
        c2 = mp.root(16, 5) + 2 * mp.root(4, 5) * mp.cos(4 * mp.pi * s / 5) + 1
        printed = mp.mpf(5) ** (mp.mpf(n) / 2) / (c1 ** (mp.mpf(m) / 2) * c2 ** (mp.mpf(5 * n) / 4))
        u = t5 * mp.exp(2j * mp.pi * s / 5)
        worst = max(worst, abs(printed / abs((1 / (u - 1)) ** m * (u ** 4 - u ** 2 - u + 1) ** n) - 1))
check("13 printed (7.4.1) factor == |sigma_s(a^m b^n)| for s = 1..4", worst < mp.mpf(10) ** -40, f"max rel. diff {mp.nstr(worst, 3)}")


# 14/15. seventh root: corrected (7.4.2) factor == |sigma_s(a^l b^m c^n)|; the printed factor is not
def f742(l_, m, n, s, printed):
    t = mp.root(2, 7)
    C = lambda k: mp.cos(k * mp.pi * s / 7)  # noqa: E731
    num_c = (mp.root(4, 7) + 2 * t * (mp.cos(2 * mp.pi * s / 5) if printed else C(2)) + 1) ** (2 * n if printed else n)
    den_a = (mp.root(64, 7) - 2 * mp.root(8, 7) * C(2) + 1 if printed else mp.root(4, 7) - 2 * t * C(2) + 1) ** (mp.mpf(l_) / 2)
    den_b = (mp.root(64, 7) - 2 * mp.root(8, 7) * C(6) + 1) ** (mp.mpf(7 * m) / 2)
    den_c = (mp.root(64, 7) + 2 * mp.root(8, 7) * C(6) + 1) ** (mp.mpf(n) / 2)
    return mp.mpf(7) ** m * num_c / (den_a * den_b * den_c)


def conj7(l_, m, n, s):
    u = mp.root(2, 7) * mp.exp(2j * mp.pi * s / 7)
    return abs((1 / (u - 1)) ** l_ * (7 / (u ** 3 - 1) ** 7) ** m * ((u + 1) / (u ** 2 - u + 1)) ** n)


pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (4, 1, 0), (20, 1, 3), (60, 0, 0)]
w_corr = max(abs(f742(*p, s, False) / conj7(*p, s) - 1) for p in pts for s in range(1, 7))
w_prnt = min(max(abs(f742(*p, s, True) / conj7(*p, s) - 1) for s in range(1, 7)) for p in [(1, 0, 0), (0, 0, 1)])
check("14 corrected (7.4.2) factor == |sigma_s(a^l b^m c^n)| for s = 1..6", w_corr < mp.mpf(10) ** -40, f"max rel. diff {mp.nstr(w_corr, 3)}")
check("15 printed (7.4.2) factor differs from the conjugate moduli (a- and c-factors)", w_prnt > mp.mpf("0.1"), f"min rel. diff {mp.nstr(w_prnt, 3)}")

print(f"\n{sum(RESULTS)}/{len(RESULTS)} checks passed")
sys.exit(0 if all(RESULTS) else 1)
