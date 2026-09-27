#!/usr/bin/env python3
"""Build and verify the two OEIS drafts (new sequence + comment for A374455). Every stated formula is asserted.
Output: oeis/new-sequence.txt, oeis/A374455-comment.txt"""
import os

import mpmath as mp
import sympy as sp

os.makedirs("oeis", exist_ok=True)
x = sp.symbols("x")


def mul5(u, v):
    w = [0] * 9
    for i, ui in enumerate(u):
        for j, vj in enumerate(v):
            w[i + j] += ui * vj
    return [w[k] + 2 * w[k + 5] if k + 5 < 9 else w[k] for k in range(5)]


N = 501
A_, TH = [1, 1, 1, 1, 1], [0, 1, 0, 0, 0]
pw = [[1, 0, 0, 0, 0]]
for _ in range(N):
    pw.append(mul5(pw[-1], A_))
TrA = [5 * p[0] for p in pw]                    # Tr(alpha^n)
TrB = [5 * mul5(p, TH)[0] for p in pw]          # Tr(theta alpha^n)
assert all(t % 5 == 0 for t in TrA) and all(t % 10 == 0 for t in TrB[1:])
S = [t // 5 for t in TrA]                       # new sequence, offset 0
Aold = [TrB[n + 1] // 10 for n in range(N - 1)]  # = A374455(n)

OEIS_A374455 = [1, 5, 35, 235, 1580, 10626, 71460, 480570, 3231845, 21734235, 146163251, 982951365, 6610371480,
                44454906580, 298960311840, 2010515259876, 13520763292345, 90927457083265, 611489327404315,
                4112280377388895, 27655184063541876, 185981775414350150, 1250731895575163300]
assert Aold[:23] == OEIS_A374455, "A374455 mismatch"

# recurrence (signature 5,10,10,5,1) for both
for seq in (S, Aold):
    assert all(seq[n] == 5 * seq[n - 1] + 10 * seq[n - 2] + 10 * seq[n - 3] + 5 * seq[n - 4] + seq[n - 5] for n in range(5, 400))

# g.f. of S: (2-(1+x)^4)/(2-(1+x)^5), and of A374455: 1/(2-(1+x)^5)
ser = sp.series((2 - (1 + x) ** 4) / (2 - (1 + x) ** 5), x, 0, 40).removeO()
assert [ser.coeff(x, k) for k in range(40)] == S[:40]
ser2 = sp.series(1 / (2 - (1 + x) ** 5), x, 0, 40).removeO()
assert [ser2.coeff(x, k) for k in range(40)] == Aold[:40]

# S(n) = A(n) - 4A(n-1) - 6A(n-2) - 4A(n-3) - A(n-4)
Ag = lambda k: Aold[k] if k >= 0 else 0  # noqa: E731
assert all(S[n] == Ag(n) - 4 * Ag(n - 1) - 6 * Ag(n - 2) - 4 * Ag(n - 3) - Ag(n - 4) for n in range(300))

# nearest-integer properties (checked n < 500 with 400 digits; proof for large n: |eps| <= 2(|s1|^n+|s2|^n) or 2θ(...))
mp.mp.dps = 600  # alpha^500 ~ 10^414: 400 digits are NOT enough (artefact seen at n >= 482)
th = mp.root(2, 5)
al = 1 / (th - 1)
assert S[2] == 9 and mp.nstr(al ** 2, 6) == '45.2259'  # example line %e
nearS = [n for n in range(N - 1) if 5 * S[n] == int(mp.nint(al ** n))]
nearA = [n for n in range(N - 2) if 10 * Aold[n] == int(mp.nint(th * al ** (n + 1)))]
failS = sorted(set(range(N - 1)) - set(nearS))
failA = sorted(set(range(N - 2)) - set(nearA))
print("5*S(n) != round(alpha^n) for n =", failS)
print("10*A374455(n) != round(theta*alpha^(n+1)) for n =", failA)
s1, s2 = 0.78821487, 0.48922479
proofS = next(n for n in range(100) if 2 * (s1 ** n + s2 ** n) < 0.5)
proofA = next(m for m in range(100) if 2 * float(th) * (s1 ** m + s2 ** m) < 0.5) - 1
print("proof threshold: S from n =", proofS, "| A374455 from n =", proofA)
assert failS and max(failS) < proofS and failA and max(failA) < proofA + 1


def data_lines(seq, limit=260):
    cur = ""
    for v in seq:
        t = (cur + "," if cur else "") + str(v)
        if len(t) > limit:
            break
        cur = t
    return cur


dS = data_lines(S)
nS = dS.count(",") + 1
first_ok_S = min(n for n in nearS if all(k in nearS for k in range(n, N - 1)))
first_ok_A = min(n for n in nearA if all(k in nearA for k in range(n, N - 2)))
new_seq = f"""%I A000000
%S A000000 {dS}
%N A000000 a(n) = trace of (2^(1/5)-1)^(-n) over Q, divided by 5.
%C A000000 Let theta = 2^(1/5) and alpha = 1/(theta-1) = 1 + theta + theta^2 + theta^3 + theta^4, a unit of Z[theta] with minimal polynomial x^5 - 5x^4 - 10x^3 - 10x^2 - 5x - 1. Then 5*a(n) = alpha^n + (sum of the other four conjugates of alpha^n).
%C A000000 alpha is a Pisot number: its other conjugates have moduli 0.7882149... and 0.4892248... (two complex pairs). Hence 5*a(n) is the integer nearest to alpha^n for n = 2 and all n >= {first_ok_S}, and the distance tends to 0 like 0.788^n.
%C A000000 This is the case n=0 (without the factor theta) of the near-integer phenomenon in Ramanujan's "third manuscript" on page 343 of the Lost Notebook, which Andrews and Berndt (2013, Sect. 7.4) could not interpret.
%C A000000 Of the numbers 1/(2^(1/p)-1), p >= 2, exactly those with p <= 6 are Pisot numbers (criterion 2^(1/p) > 2*cos(2*Pi/p)).
%D A000000 G. E. Andrews and B. C. Berndt, Ramanujan's Lost Notebook, Part IV, Springer, 2013, Section 7.4, pp. 179-181.
%F A000000 a(n) = 5*a(n-1) + 10*a(n-2) + 10*a(n-3) + 5*a(n-4) + a(n-5) for n >= 5.
%F A000000 G.f.: (2 - (1+x)^4)/(2 - (1+x)^5) = (1 - 4*x - 6*x^2 - 4*x^3 - x^4)/(1 - 5*x - 10*x^2 - 10*x^3 - 5*x^4 - x^5).
%F A000000 a(n) = A374455(n) - 4*A374455(n-1) - 6*A374455(n-2) - 4*A374455(n-3) - A374455(n-4), with A374455(k) = 0 for k < 0.
%e A000000 a(2) = 9: (2^(1/5)-1)^(-2) = 45.2259..., and 5*9 = 45.
%o A000000 (Python)
%o A000000 def a_list(N):
%o A000000     s = [1, 1, 9, 61, 409]
%o A000000     while len(s) < N: s.append(5*s[-1] + 10*s[-2] + 10*s[-3] + 5*s[-4] + s[-5])
%o A000000     return s[:N]
%Y A000000 Cf. A374455 (1/(2-(1+x)^5)), A001333 (the analog for 2^(1/2): Tr((1+sqrt(2))^n)/2).
%K A000000 nonn,easy
%O A000000 0,3
%A A000000 _Otakhon U. Kenjaev_, Sep 27 2026
"""
open("oeis/new-sequence.txt", "w").write(new_seq)

comment = f"""Proposed additions to A374455 (Expansion of o.g.f. 1/(1 - 5*x - 10*x^2 - 10*x^3 - 5*x^4 - x^5)):

%C A374455 a(n) = Tr(theta*alpha^(n+1))/10, where theta = 2^(1/5) and alpha = 1/(theta-1) = 1 + theta + ... + theta^4 (a Pisot unit). Consequently 10*a(n) is the integer nearest to 2^(1/5)*(2^(1/5)-1)^(-(n+1)) for n = 3 and all n >= {first_ok_A}; this is the case n=0 of Ramanujan's claim on page 343 of the Lost Notebook (see Andrews-Berndt, Part IV, Sect. 7.4). - _Otakhon U. Kenjaev_, Sep 27 2026
%F A374455 G.f.: 1/(2 - (1+x)^5). - _Otakhon U. Kenjaev_, Sep 27 2026

Verification (make_oeis.py): our values equal the 23 terms in the DATA field of the entry; both sides satisfy the recurrence with
signature (5,10,10,5,1) and agree on the first 5 terms, hence for all n. Nearest-integer property checked for n < 499
(exceptions n = {failA}) and proved for n >= {proofA} from |eps| <= 2*2^(1/5)*(0.78822^(n+1) + 0.48923^(n+1)).
"""
open("oeis/A374455-comment.txt", "w").write(comment)
assert 2 in nearS and 3 not in nearS and 3 in nearA and 4 not in nearA  # wording of the comments
print("new sequence: terms in %S =", nS, "| first_ok_S =", first_ok_S, "| first_ok_A =", first_ok_A)
print("OK: all assertions passed")
