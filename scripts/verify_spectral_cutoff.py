#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verifier for Part XXII: absolute convergence of the global correlation
R(j) = prod_{q>3} R_q(j), the finiteness-of-resonance theorem, and the spectral
cutoff k*(epsilon).

Checks:
  (1) For a fixed lag j, the resonance primes (States A,B with O(1/q) deviation)
      are EXACTLY the primes dividing j or 3j+-1, all <= 3j+1; every prime beyond
      is State C with deviation 4/(q-2)^2 = O(1/q^2).
  (2) The full product R(1) = S_quad / S_twin^2 = 0.3969 (j=1, all State C).
  (3) The k*(epsilon) table (smallest cutoff prime within relative error epsilon).
Standard library only; also writes the data CSVs when run from the repo.
"""
import math, os, csv

def sieve(limit):
    s = bytearray([1]) * (limit + 1); s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            for k in range(i*i, limit + 1, i): s[k] = 0
    return [p for p in range(5, limit + 1) if s[p]]

def Rq(q, j):
    a = pow(6, -1, q); U = {a, (-a) % q, (a - j) % q, (-a - j) % q}
    return (q - len(U)) / q / ((q - 2) / q) ** 2

def state(q, j):
    inv3 = pow(3, -1, q)
    if j % q == 0: return "A"
    if j % q in (inv3 % q, (-inv3) % q): return "B"
    return "C"

def check_resonance(j, primes):
    res = [q for q in primes if state(q, j) != "C"]
    pred = sorted({q for q in primes if j % q == 0 or (3*j-1) % q == 0 or (3*j+1) % q == 0})
    ok = (res == pred) and all(q <= 3*j + 1 for q in res)
    print(f"(1) lag j={j}: resonance primes {res}; predicted (q|j or q|3j+-1) {pred}; "
          f"all <= 3j+1={3*j+1}: {ok}")
    return ok

def kstar_rows(limit=200000, targets=(1e-1,1e-2,1e-3,1e-4)):
    primes = sieve(limit)
    ln_full = sum(math.log(Rq(q,1)) for q in primes)
    R1 = math.exp(ln_full)
    rows, cum, kept, ti = [], 0.0, 0, 0
    for q in primes:
        cum += math.log(Rq(q,1)); kept += 1
        err = abs(math.exp(ln_full - cum) - 1)
        while ti < len(targets) and err < targets[ti]:
            rows.append((targets[ti], q, kept, err)); ti += 1
        if ti >= len(targets): break
    return R1, rows

def write_data():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(os.path.dirname(here), "data"); os.makedirs(out, exist_ok=True)
    R1, rows = kstar_rows()
    with open(os.path.join(out, "kstar_table.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["epsilon","k_star","primes_kept","error_at_kstar"])
        for eps,k,n,e in rows: w.writerow([f"{eps:.0e}", k, n, f"{e:.3e}"])
    with open(os.path.join(out, "R1_product.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["quantity","value","note"])
        w.writerow(["R(1)", f"{R1:.6f}", "global product prod_q R_q(1) = S_quad/S_twin^2 (all State C)"])
    print(f"\nwrote data/kstar_table.csv, data/R1_product.csv  (R(1)={R1:.6f})")

if __name__ == "__main__":
    P = sieve(2000)
    a = all(check_resonance(j, P) for j in (1, 7, 30, 100))
    R1, rows = kstar_rows()
    print(f"\n(2) R(1) = {R1:.6f}  (expected 0.396881)")
    print("(3) k*(epsilon):")
    for eps,k,n,e in rows: print(f"    eps={eps:.0e}: k*={k:<6} primes={n:<5} err={e:.2e}")
    write_data()
    print("\nRESONANCE-FINITENESS CHECKS PASS:", a)
