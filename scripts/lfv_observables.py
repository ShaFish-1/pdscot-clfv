# -*- coding: utf-8 -*-
"""lfv_observables.py — 带电轻子味破坏(LFV)可观测量全集合
arXiv:2605.05312 准 Dirac scotogenic, 三对模型

输入 = 拉氏量参数:  M_i, m(eta+), y_i, eps 矩阵 (Eq.30/32)
输出 = BR(mu->e gamma), BR(tau->e/mu gamma), BR(l->3l), mu-e 转换率,
       与实验上限的比值表。纯 Python (cmath), 无依赖。

用法: python3 lfv_observables.py
"""
import math

# ==================== 输入: 拉氏量参数 ====================
M   = [200.0, 500.0, 1000.0]      # 三对 Majorana 质量 [GeV]
mEp = 225.0                        # m(eta+) [GeV]
y   = [0.4, 0.4, 0.4]              # 对角 Yukawa (Eq.30)

# eps 矩阵 (Eq.32 默认; 行 = alpha: e,mu,tau;  列 = i: pair 1,2,3)
eps = [[1.24e-9, 0, (0.018 + 0.018j)*1e-9],
       [(0.024 - 0.027j)*1e-9, (0.56 - 0.002j)*1e-9, 0],
       [0, 0.089e-9, (0.20 - 0.002j)*1e-9]]
# =========================================================

alpha = 1/127.9
GF    = 1.1663787e-5               # GeV^-2
ml    = [5.11e-4, 0.10566, 1.7768] # e, mu, tau [GeV]
BRnu  = {(0,1): 1.0, (2,0): 0.1782, (2,1): 0.1739}   # BR(l_a -> l_b nu nu)

def F2(z):
    return (1/6 - z + z**2/2 + z**3/3 - z**2*math.log(z)) / (1-z)**4

# ---------- 1. 振幅 S[alpha,beta] (alpha != beta) ----------
#   S_ab = sum_i y_i*( delta_{a,i}*eps*[b,i] + delta_{b,i}*eps[a,i] ) * F2(M_i^2/mEp^2)
S = {}
for a in range(3):
    for b in range(3):
        if a == b: continue
        s = 0j
        for i in range(3):
            term = 0j
            if a == i: term += y[i] * eps[b][i].conjugate()
            if b == i: term += y[i] * eps[a][i]
            s += term * F2((M[i]/mEp)**2)
        S[(a,b)] = s

# ---------- 2. BR(l_a -> l_b gamma)  (Eq.33 推广) ----------
def BR_lg(a, b):
    return 3*alpha/(64*math.pi*GF**2*mEp**4) * abs(S[(a,b)])**2 * BRnu[(a,b)]

# ---------- 3. 三体 BR(l_a -> 3 l_b) ~ BR(l_a->l_b gamma) * (alpha/3pi)[ln(m_a^2/m_b^2)-11/4] ----------
def BR_3l(a, b):
    return BR_lg(a, b) * alpha/(3*math.pi) * (math.log(ml[a]**2/ml[b]**2) - 11/4)

# ---------- 4. mu-e 转换 (偶极主导):  CR(N) ~ BR(mu->e gamma) * D_N ----------
D_conv = {"Al": 2.6e-3, "Ti": 4.3e-3, "Au": 2.5e-3}

# ---------- 5. 输出 ----------
print("="*72)
print("LFV 可观测量 (输入: M=(200,500,1000) GeV, m_ep=225, y=0.4, eps=Eq.32)")
print("="*72)

name = {(0,1): "mu -> e", (2,0): "tau -> e", (2,1): "tau -> mu"}

bounds = {"mu -> e gamma": 3.1e-13, "tau -> e gamma": 3.3e-8, "tau -> mu gamma": 4.4e-8,
          "mu -> 3e": 1.0e-12, "tau -> 3e": 2.7e-8, "tau -> 3mu": 2.1e-8,
          "muTi -> eTi": 4.3e-12, "muAu -> eAu": 7.0e-13, "muAl -> eAl": None}

rows = []
for (a,b) in [(0,1),(2,0),(2,1)]:
    rows.append(("%s gamma" % name[(a,b)], BR_lg(a,b), bounds["%s gamma" % name[(a,b)]]))
for (a,b) in [(0,1),(2,0),(2,1)]:
    if (a,b)==(0,1):
        rows.append(("mu -> 3e", BR_3l(0,1), bounds["mu -> 3e"]))
    elif (a,b)==(2,0):
        rows.append(("tau -> 3e", BR_3l(2,0), bounds["tau -> 3e"]))
    else:
        rows.append(("tau -> 3mu", BR_3l(2,1), bounds["tau -> 3mu"]))

for nuc in ["Al","Ti","Au"]:
    rows.append(("mu%s -> e%s" % (nuc,nuc), BR_lg(0,1)*D_conv[nuc], bounds["mu%s -> e%s" % (nuc,nuc)]))

print("%-16s %-14s %-14s %s" % ("observable", "prediction", "bound", "pred/bound"))
print("-"*72)
for nm, val, bd in rows:
    ratio = "n/a" if bd is None else "%.2e" % (val/bd)
    print("%-16s %-14.2e %-14s %s" % (nm, val, ("%.1e" % bd) if bd else "future (1e-16)", ratio))

print("\n振幅 |S|: " + "  ".join("%s=%.2e" % (name[k], abs(v)) for k,v in S.items()))
print("\n注: BR(mu->e gamma) 论文值 ~1e-26; 三体/转换用偶极主导近似。")
