# -*- coding: utf-8 -*-
"""kappa_check.py — 验证论文 Eq.(10) 的 ΔM 项系数：精确单圈 vs κ=1 近似

精确:  m = Λ(M1)[(YY^T) + κ (ΔM1/M1) Y1'Y1'],  Λ = M*G(M),
       κ = d lnΛ / d lnM = 1 + M*G'(M)/G(M)
论文:  κ 被默认成 1（m_eta >> M 极限）。κ<1 时论文的 ΔM 项偏大 1/κ 倍。

用法: python3 kappa_check.py
"""
import math

def G(M, mR, mI):
    """圈函数（与 00_Conventions.md §5 的 G 相同）"""
    return (mR**2/(mR**2-M**2)*math.log(mR**2/M**2)
            - mI**2/(mI**2-M**2)*math.log(mI**2/M**2))

def Gp(M, mR, mI):
    """G'(M): f = m² ln(m²/M²)/(m²−M²)，f' = m²[−2/(M(m²−M²)) + 2M ln(m²/M²)/(m²−M²)²]"""
    def fp(m):
        u = math.log(m**2/M**2)
        d = m**2 - M**2
        return m**2*(-2/(M*d) + 2*M*u/d**2)
    return fp(mR) - fp(mI)

def kappa(M, mR, mI):
    return 1 + M*Gp(M, mR, mI)/G(M, mR, mI)

pts = {
 "论文基准  M=200,  mR=296,  mI=211":
     (200.0, 296.0, 211.0),
 "S1        M=2561, mR=2567.47, mI=2567.40":
     (2561.318289954, 2567.468979424, 2567.401443985),
 "S2        M=4974, mR=4976.39, mI=4976.36":
     (4973.658807576, 4976.391334186, 4976.359653618),
 "S3        M=3322, mR=3325.59, mI=3324.98":
     (3321.992692072, 3325.594188489, 3324.980441167),
 "S4        M=3856, mR=3860.49, mI=3860.30":
     (3855.595868081, 3860.493126902, 3860.300869081),
}

print(f"{'point':44s} {'mR/M':>8s} {'G(M)':>13s} {'kappa':>8s} {'论文项/精确项':>14s}")
print("="*95)
for k, (M, mR, mI) in pts.items():
    kap = kappa(M, mR, mI)
    g = G(M, mR, mI)
    print(f"{k:44s} {mR/M:8.4f} {g:13.5e} {kap:8.4f} {1/kap:14.3f}")

print()
print("解读:")
print(" 1/kappa = 论文的 ΔM 项相对精确结果偏大的倍数（κ<1 时 >1）")
print(" 论文 Eq.(10) 写 ΔM1/M1 等价于假设 κ=1 —— 只在 m_eta >> M 时成立")
