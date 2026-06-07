#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure for Part XXII: finite resonance vs O(1/q^2) tail, and the k*(eps) cutoff."""
import math, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def sieve(limit):
    s=bytearray([1])*(limit+1); s[0]=s[1]=0
    for i in range(2,int(limit**0.5)+1):
        if s[i]:
            for k in range(i*i,limit+1,i): s[k]=0
    return [p for p in range(5,limit+1) if s[p]]

def Rq(q,j):
    a=pow(6,-1,q); U={a,(-a)%q,(a-j)%q,(-a-j)%q}
    return (q-len(U))/q/((q-2)/q)**2

fig,ax=plt.subplots(1,2,figsize=(13.5,5.2))
fig.suptitle("Spectral structure of the global correlation $R(j)=\\prod_q R_q(j)$: "
             "finite resonance, $O(1/q^2)$ tail",fontsize=12.5,fontweight="bold")

# Panel A: |R_q(j)-1| vs q for j=30; resonances (State A/B) above the State-C floor
a=ax[0]; j=30; P=sieve(200)
dev=[abs(Rq(q,j)-1) for q in P]
a.loglog(P,dev,"o",ms=4,color="#1f4e79",label=f"$|R_q({j})-1|$")
qq=np.array([5]+P); a.loglog(qq,4/(qq-2)**2,"-",color="#b4341f",lw=1,
        label=r"State-C floor $4/(q-2)^2 \sim O(1/q^2)$")
# annotate resonance primes (q | j or q | 3j+-1)
res=[q for q in P if (j%q==0) or ((3*j-1)%q==0) or ((3*j+1)%q==0)]
a.loglog(res,[abs(Rq(q,j)-1) for q in res],"s",ms=9,mfc="none",mec="#2e7d32",
         label="resonance: $q\\,|\\,j$ or $q\\,|\\,3j{\\pm}1$")
a.set_xlabel("prime $q$"); a.set_ylabel(r"$|R_q(j)-1|$")
a.set_title(f"(A) $j={j}$: finite resonances above the $O(1/q^2)$ tail")
a.legend(fontsize=8); a.grid(alpha=.3,which="both")

# Panel B: truncation error vs cutoff prime k*, j=1 (all State C)
b=ax[1]; primes=sieve(20000)
ln_full=sum(math.log(Rq(q,1)) for q in primes)
ks=[]; errs=[]; cum=0.0
for q in primes:
    cum+=math.log(Rq(q,1)); ks.append(q); errs.append(abs(math.exp(ln_full-cum)-1))
b.loglog(ks,errs,"-",color="#1f4e79",lw=1.3)
table=[(1e-1,13),(1e-2,73),(1e-3,557),(1e-4,4229)]
for eps,kk in table:
    b.axhline(eps,color="0.7",ls=":",lw=.8)
    b.plot(kk,eps,"o",color="#b4341f",ms=6)
    b.annotate(f"$k^\\star={kk}$",(kk,eps),fontsize=8,xytext=(kk*1.1,eps*1.4))
b.set_xlabel("cutoff prime $k^\\star$"); b.set_ylabel("tail truncation error")
b.set_title("(B) $k^\\star(\\varepsilon)$ for $R(1)$ (all State C)")
b.grid(alpha=.3,which="both")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig("fig_spectral_cutoff.png",dpi=200); fig.savefig("fig_spectral_cutoff.pdf")
print("wrote fig_spectral_cutoff.{png,pdf}; resonance primes for j=30:",res)
