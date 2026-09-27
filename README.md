# Page 343 of Ramanujan's Lost Notebook — Pisot units

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23003128.svg)](https://doi.org/10.5281/zenodo.23003128)

Andrews and Berndt (*Ramanujan's Lost Notebook, Part IV*, Sect. 7.4) could not give a meaning to the "third manuscript" on
page 343 and wrote that its claims are wrong. This repository contains a short note showing that the first claim is correct:
with θ = 2^(1/5), the numbers a = 1/(θ−1) and b = √5/(1+θ²)^(5/2) are units of Z[θ], Ramanujan's error factors are the moduli of
Galois conjugates, and the integer is Tr(a^m b^n θ).

- **Note:** [`paper/Kenjaev_Ramanujan_p343.pdf`](paper/Kenjaev_Ramanujan_p343.pdf) (6 pages)
- **Calculator:** https://selloriybrendi.github.io/ramanujan-page-343/ (runs offline in the browser; built-in self-test)
- **Verification:** `python3 verify/verify_343.py` → `15/15 checks passed` (requires only SymPy and mpmath)
- **Figure:** `python3 verify/make_figure.py` (needs matplotlib)
- **OEIS drafts:** [`oeis/`](oeis/) — a new sequence Tr((2^(1/5)−1)^(−n))/5 = 1, 1, 9, 61, 409, … and a proposed comment for A374455

## Results
1. Theorem 3: |a^m b^n θ − Tr(a^m b^n θ)| ≤ 2θ(|σ₁(a^m b^n)| + |σ₂(a^m b^n)|); the bound tends to 0 exactly when m − 6.0114 n → ∞.
2. Proposition 5: the constant 2θ is sharp.
3. Lemma 6: 1/(2^(1/p) − 1) is a Pisot number if and only if p ≤ 6 — Ramanujan's two cases 5 and 7 lie on both sides.
4. Theorem 7: the printed (7.4.2) differs from the conjugate moduli in three places; the corrected form holds on an explicit cone.

## Cite
Kenjaev, O. U. (2026). *Pisot units and the third manuscript on page 343 of Ramanujan's Lost Notebook, with a verification script and calculator* (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.23003128

(10.5281/zenodo.23003128 always resolves to the latest version; v1.0.0 is 10.5281/zenodo.23003129.)

## Use of generative AI
Claude Opus 5.5 (Anthropic, `claude-opus-5-5`) via Claude Code was used for the computations (SymPy, mpmath) and to help prepare
the text, under the author's direction; every number was verified independently. The author takes full responsibility.

Author: Otakhon U. Kenjaev, ORCID [0009-0009-3566-9285](https://orcid.org/0009-0009-3566-9285)
