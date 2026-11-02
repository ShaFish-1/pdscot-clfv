# -*- coding: utf-8 -*-
"""lfv_verify_S1S4_full.py (v6) — S1-S4 全部 CLFV 过程验证（修正版）

v6 修正（按参考方检查意见逐条落实）:
  Bug1: ΔF 符号 — 用精确差分 (F2(x_n)-F2(x_p))/2（不再用带 (x_p-x_n) 的导数近似）;
  Bug2: 圈函数改 mpmath 30 位直接逐态求和（不做解析导数近似; F2 本身不再有消减误差）;
  相位: μ->eγ 改调 A_D(tag,1,0) = sum H*_e H_mu（参考约定 H*_beta H_alpha, alpha=mu,beta=e）;
  诚实标注: 独立验证 = A_D/A_ND/辐射道 BR/κ_D/LNC-LNV;
            非独立 = BR(mu->3e) 的 box 块(B 读自参考)与 CR(R_N 取自参考 CR/BR 列)。

约定 = 00_Conventions.md §5/§8; BR 用 3(4pi)^3 a/(4 GF^2)|A_D|^2。
依赖: mpmath (Ubuntu 一般自带; 若无: pip3 install mpmath)。用法: python3 lfv_verify_S1S4_full.py
"""
import math
try:
    import mpmath as mp
    mp.mp.dps = 30
    HAVE_MP = True
except ImportError:
    HAVE_MP = False

alpha0 = 1/137.035999084
GF     = 1.1663787e-5
me, mmu, mta = 5.1099895069e-4, 0.1056583755, 1.77686
BRnu = {(0,1): 1.0, (1,0): 1.0, (2,0): 0.1782, (0,2): 0.1782,
        (2,1): 0.1739, (1,2): 0.1739}

def F2(z):
    """F2 用 mpmath 30 位 (双精度兜底加近 1 级数保护)"""
    if HAVE_MP:
        z = mp.mpf(z)
        return (mp.mpf(1)/6 - z + z**2/2 + z**3/3 - z**2*mp.log(z)) / (1-z)**4
    d = 1 - z
    if abs(d) < 1e-2:
        return 1/12 + d/30 + d**2/60
    return (1/6 - z + z**2/2 + z**3/3 - z**2*math.log(z)) / (1-z)**4

def G2(z):
    if HAVE_MP:
        z = mp.mpf(z)
        return (2 - 9*z + 18*z**2 - 11*z**3 + 6*z**3*mp.log(z)) / (6*(1-z)**4)
    d = 1 - z
    if abs(d) < 1e-2:
        return 1/4 + d/20 + d**2/30
    return (2 - 9*z + 18*z**2 - 11*z**3 + 6*z**3*math.log(z)) / (6*(1-z)**4)

# ================= 四个点的输入 (参考文件 1.1 / 2.3) =================
points = {
 "S1": dict(M1=2561.318289954, M2=6403.295724886, M3=12806.591449771,
            dM1_keV=258.951, mHp=2564.409492749, y1=1.93381490957),
 "S2": dict(M1=4973.658807576, M2=12434.14701894, M3=24868.294037879,
            dM1_keV=255.878, mHp=4978.500953041, y1=2.81452893374),
 "S3": dict(M1=3321.992692072, M2=8304.98173018, M3=16609.96346036,
            dM1_keV=253.293, mHp=3324.048967899, y1=2.2740354477),
 "S4": dict(M1=3855.595868081, M2=9638.989670202, M3=19277.979340404,
            dM1_keV=263.4, mHp=3860.326346962, y1=2.43972269808),
}
y2 = y3 = 0.4

eps = {
 "S1": [[2.0927e-8+1.2079e-9j, -7.4694e-9-1.486e-8j,  -3.9696e-8-1.8102e-8j],
        [-1.7404e-9-3.4624e-9j, 1.3953e-7-1.6197e-9j,  1.2721e-7-3.0019e-10j],
        [-8.0626e-9-3.6768e-9j, 1.1089e-7-2.6168e-10j, 1.6999e-7+1.4578e-9j]],
 "S2": [[1.8784e-8+1.7681e-9j, -1.5925e-8-3.1682e-8j, -8.4672e-8-3.8613e-8j],
        [-2.5476e-9-5.0683e-9j, 2.9748e-7-3.4533e-9j,  2.7133e-7-6.4031e-10j],
        [-1.1802e-8-5.382e-9j,  2.3641e-7-5.5791e-10j, 3.6259e-7+3.1094e-9j]],
 "S3": [[1.4897e-8+1.1298e-10j, -8.22e-10-1.6353e-9j,  -4.37e-9-1.9928e-9j],
        [-1.6278e-10-3.2384e-10j, 1.5354e-8-1.7825e-10j, 1.4004e-8-3.3047e-11j],
        [-7.5411e-10-3.4389e-10j, 1.2203e-8-2.8797e-11j, 1.8712e-8+1.6048e-10j]],
 "S4": [[1.5188e-8+3.362e-10j, -2.624e-9-5.2203e-9j,    -1.3949e-8-6.3612e-9j],
        [-4.844e-10-9.6368e-10j, 4.9016e-8-5.6901e-10j,  4.4701e-8-1.0549e-10j],
        [-2.2441e-9-1.0233e-9j, 3.8955e-8-9.1928e-11j,  5.9734e-8+5.1226e-10j]],
}

# 参考值 (第 3/4 节 + 3.1 节)
ref = {
 "S1": dict(mueg=2.4776e-40, m3e=9.4695e-40, Au=8.8901e-43, Ti=1.084e-42,
            Al=6.0228e-43, taueg=9.7397e-41, taumg=1.5957e-39,
            BD=complex(1.2348e-23, -9.2167e-24), BM=complex(-7.7048e-25, -1.3012e-24),
            AD=complex(-4.7676e-26, 2.8819e-26), AND=complex(-8.4997e-26, 6.5706e-26),
            kappaD=1.3958e-7, LNC=1.03984, LNV=0.218155),
 "S2": dict(mueg=3.8627e-40, m3e=6.2641e-39, Au=1.234e-42, Ti=1.5046e-42,
            Al=8.3596e-43, taueg=8.8727e-41, taumg=2.3122e-39,
            BD=complex(2.9998e-23, -2.6543e-23), BM=complex(-6.2951e-25, -9.7339e-25),
            AD=complex(-5.1129e-26, 4.7163e-26), AND=complex(-9.8679e-26, 8.6671e-26),
            kappaD=3.0832e-7, LNC=1.01497, LNV=0.0502557),
 "S3": dict(mueg=4.0265e-44, m3e=8.7908e-44, Au=8.6758e-47, Ti=1.0579e-46,
            Al=5.8776e-47, taueg=5.1196e-44, taumg=8.2895e-44,
            BD=complex(1.6615e-25, -5.9926e-27), BM=complex(-3.6053e-26, -7.0058e-26),
            AD=complex(-6.2099e-28, -3.4458e-28), AND=complex(-7.7488e-28, 1.5042e-28),
            kappaD=2.7193e-8, LNC=0.563061, LNV=0.844819),
 "S4": dict(mueg=7.482e-43, m3e=6.9969e-42, Au=2.8431e-45, Ti=3.4667e-45,
            Al=1.9261e-45, taueg=4.841e-43, taumg=4.7258e-42,
            BD=complex(1.1277e-24, -6.8398e-25), BM=complex(-9.793e-26, -1.8222e-25),
            AD=complex(-2.9193e-27, 9.2188e-28), AND=complex(-4.8155e-27, 3.2045e-27),
            kappaD=4.9493e-8, LNC=1.01974, LNV=0.415675),
}

# 参考文件自身的 CR/BR 逐点比值 (Au 靶; 供展示 — CR 非独立验证)
R_ref_Au = {"S1": 8.8901e-43/2.4776e-40, "S2": 1.234e-42/3.8627e-40,
            "S3": 8.6758e-47/4.0265e-44, "S4": 2.8431e-45/7.482e-43}

def model(tag):
    p = points[tag]
    M = [p["M1"], p["M1"]+p["dM1_keV"]*1e-6, p["M2"], p["M2"]+1e-6,
         p["M3"], p["M3"]+1e-6]
    y = [p["y1"], y2, y3]
    ep = eps[tag]
    H = [[0j]*6 for _ in range(3)]
    for i in range(3):
        for a in range(3):
            d = 1.0 if a == i else 0.0
            H[a][2*i]   = (y[i]*d + ep[a][i]) / math.sqrt(2)
            H[a][2*i+1] = 1j*(y[i]*d - ep[a][i]) / math.sqrt(2)
    x = [(m/p["mHp"])**2 for m in M]
    return H, x, y

def A_D(tag, a, b):
    """§5: A_D^{ba} = sum_I H*_{bI} H_{aI} F2(x_I) / (2(4pi)^2 m_ep^2)  [GeV^-2]
    调用约定: (a,b)=(alpha,beta), 即 mu->e gamma 用 (1,0)。全程 mpmath 累加。"""
    H, x, _ = model(tag)
    m2 = points[tag]["mHp"]**2
    if HAVE_MP:
        S = mp.mpc(0)
        for I in range(6):
            S += (mp.mpc(H[b][I]).conjugate() * mp.mpc(H[a][I])) * F2(x[I])
        S = complex(S)
    else:
        S = 0j
        for I in range(6):
            S += H[b][I].conjugate()*H[a][I]*F2(x[I])
    return S / (2*(4*math.pi)**2*m2)

def A_ND(tag, a, b):
    H, x, _ = model(tag)
    m2 = points[tag]["mHp"]**2
    if HAVE_MP:
        S = mp.mpc(0)
        for I in range(6):
            S += (mp.mpc(H[b][I]).conjugate() * mp.mpc(H[a][I])) * G2(x[I])
        S = complex(S)
    else:
        S = 0j
        for I in range(6):
            S += H[b][I].conjugate()*H[a][I]*G2(x[I])
    return S / (6*(4*math.pi)**2*m2)

def kappa_D(tag, a, b):
    H, x, _ = model(tag)
    if HAVE_MP:
        S = mp.mpc(0)
        tot = mp.mpf(0)
        for I in range(6):
            t = (mp.mpc(H[b][I]).conjugate() * mp.mpc(H[a][I])) * F2(x[I])
            S += t
            tot += abs(t)
        return float(abs(S)/tot)
    terms = [H[b][I].conjugate()*H[a][I]*F2(x[I]) for I in range(6)]
    return abs(sum(terms))/sum(abs(t) for t in terms)

def LNC_LNV(tag, a, b):
    """3.1/§8: 每对 i: LNC=(y_b* y_a + eps_b* eps_a)Fbar, LNV=(y_b* eps_a + eps_b* y_a) dF,
    Fbar=(F2(xn)+F2(xp))/2, dF=(F2(xn)-F2(xp))/2  [Bug1 修正: 精确差分, n-p 方向]"""
    H, x, y = model(tag)
    ep = eps[tag]
    sLNC = sLNV = 0j
    for i in range(3):
        ya = y[i] if a == i else 0.0
        yb = y[i] if b == i else 0.0
        Fb = (complex(F2(x[2*i])) + complex(F2(x[2*i+1])))/2
        dF = (complex(F2(x[2*i])) - complex(F2(x[2*i+1])))/2
        sLNC += (yb.conjugate()*ya + ep[b][i].conjugate()*ep[a][i])*Fb
        sLNV += (yb.conjugate()*ep[a][i] + ep[b][i].conjugate()*ya)*dF
    S = sLNC + sLNV
    return abs(sLNC)/abs(S), abs(sLNV)/abs(S)

def BR_lg(tag, a, b):
    key = (min(a, b), max(a, b))
    return 3*(4*math.pi)**3*alpha0/(4*GF**2) * abs(A_D(tag, a, b))**2 * BRnu[key]

def photonic_dip_3body(BRlg, m1, m2):
    return BRlg * alpha0/(3*math.pi) * (math.log(m1**2/m2**2) - 11/4)

print("="*78)
print("S1-S4 全部 CLFV 过程验证 (v6 mpmath+sign-fix)")
print("mpmath 可用: %s" % HAVE_MP)
print("="*78)
for tag in ["S1","S2","S3","S4"]:
    r = ref[tag]
    mueg  = BR_lg(tag, 1, 0)          # mu->e gamma: (alpha,beta)=(mu,e) -> H*_e H_mu
    taueg = BR_lg(tag, 2, 0)          # tau->e gamma: H*_e H_tau
    taumg = BR_lg(tag, 2, 1)          # tau->mu gamma: H*_mu H_tau
    m3e_dip = photonic_dip_3body(mueg, mmu, me)
    m3e_box = math.pi**2*alpha0**2/GF**2 * abs(r["BD"]+r["BM"])**2
    m3e = m3e_dip + m3e_box
    kD = kappa_D(tag, 1, 0)
    lC, lV = LNC_LNV(tag, 1, 0)

    print("\n[%s]" % tag)
    print("  %-22s %-12s %-12s %-8s" % ("observable", "our", "ref", "ratio"))
    print("  " + "-"*58)
    for nm, our, rv in [
        ("BR(mu->e gamma)", mueg, r["mueg"]),
        ("BR(tau->e gamma)", taueg, r["taueg"]),
        ("BR(tau->mu gamma)", taumg, r["taumg"]),
        ("BR(mu->3e) [半独立]", m3e, r["m3e"])]:
        print("  %-22s %-12.3e %-12.3e %-8.3g" % (nm, our, rv, our/rv))
    print("  --- 振幅相位对照 (独立验证, §8 先比振幅) ---")
    ad, and_ = A_D(tag, 1, 0), A_ND(tag, 1, 0)
    print("  A_D : %.4e %+.4e i  (ref %.4e %+.4e i, |.| ratio %.3g)" %
          (ad.real, ad.imag, r["AD"].real, r["AD"].imag, abs(ad)/abs(r["AD"])))
    print("  A_ND: %.4e %+.4e i  (ref %.4e %+.4e i, |.| ratio %.3g)" %
          (and_.real, and_.imag, r["AND"].real, r["AND"].imag, abs(and_)/abs(r["AND"])))
    print("  --- 诊断量 (独立, mu-e 道) ---")
    print("  kappa_D        %-12.3e (ref %.3e, ratio %.3g)" % (kD, r["kappaD"], kD/r["kappaD"]))
    print("  |Sum LNC|/|A_D| %-12.4f (ref %.4f)" % (lC, r["LNC"]))
    print("  |Sum LNV|/|A_D| %-12.4f (ref %.4f)" % (lV, r["LNV"]))

print("\n" + "="*78)
print("验证口径 (诚实标注)")
print("="*78)
print(" 独立验证 (从模型参数 + §5/§8 公式从头算, 与参考文件对表):")
print("   A_D/A_ND 振幅(含相位), 三个辐射道 BR, kappa_D, LNC/LNV 比值。")
print(" 半独立: BR(mu->3e) 的偶极块自算, box 块 B_D/B_M 读自参考文件(闭式未实现);")
print("  干涉项 ∝|A_D||B| 未算(参考文件估计其对 B->-B 变体的影响 ≤0.2%)。")
print(" 非独立: CR(mu-e) — R_N 逐点取自参考文件 CR/BR 列(参考值 Au 靶):")
for tg in ["S1","S2","S3","S4"]:
    print("   %s: R_N(Au)=%.4g   (四点不恒定 → 固定 R_N 不可靠; CR 需 A_ND/Z/box 才可独立计算)" % (tg, R_ref_Au[tg]))
print(" 精度: 5 位 eps 的相对误差 1e-5; kappa_D~1e-8 的相消残差相对误差 ~1e3 (S3/S4 的已知极限)。")
