# -*- coding: utf-8 -*-
"""solve_eps.py - 由振荡数据自动求解 epsilon 矩阵 (arXiv:2605.05312 三对模型)

原理: m_nu = sum_i C_i T^(i)(eps) 对 eps 线性 -> 6 个复方程 = 线性方程组;
3 个纹理零 -> 6x6 复线性系统, numpy 直接解。

用法: python3 solve_eps.py
依赖: numpy (若无: pip3 install numpy)
"""
import numpy as np

# ==================== 输入区 (改成你的基准点/数据) ====================
M    = np.array([200., 500., 1000.])   # 三对 Majorana 质量 [GeV]
dM   = np.array([1e-6, 1e-6, 1e-6])    # 质量劈裂 [GeV] (1 keV = 1e-6 GeV)
mR, mI = 296., 211.                     # 惰性标量 [GeV]
y    = np.array([0.4, 0.4, 0.4])        # 对角 Yukawa (Eq.30)

t12, t23, t13 = 0.587, 0.819, 0.147     # 振荡角 [rad]
delta = 1.5*np.pi                       # Dirac 相位 (论文取 3pi/2)
m1, m2, m3 = 0.005, 0.01, 0.05          # 质量 [eV] (m1 为人为选的标度)

texture_zeros = [(2,0), (0,1), (1,2)]   # 纹理零 (alpha,i): 论文取 eps_tau1, eps_e2, eps_mu3
# ======================================================================

# ---------- 1. PMNS 与数据侧矩阵 ----------
c12, s12 = np.cos(t12), np.sin(t12)
c23, s23 = np.cos(t23), np.sin(t23)
c13, s13 = np.cos(t13), np.sin(t13)
U = np.array([
    [c12*c13, s12*c13, s13*np.exp(-1j*delta)],
    [-s12*c23 - c12*s23*s13*np.exp(1j*delta),
      c12*c23 - s12*s23*s13*np.exp(1j*delta), s23*c13],
    [ s12*s23 - c12*c23*s13*np.exp(1j*delta),
     -c12*s23 - s12*c23*s13*np.exp(1j*delta), c23*c13]])
UDU = U @ np.diag([m1, m2, m3]) @ U.T    # [eV]

# ---------- 2. C_i 与 epsilon_{M,i} ----------
def B(x):   # x = M^2 [GeV^2]
    return mR**2/(mR**2-x)*np.log(mR**2/x) - mI**2/(mI**2-x)*np.log(mI**2/x)
C  = M/(16*np.pi**2) * B(M**2)            # [GeV]
eM = dM/(2*M)
Cv = C*1e9                                 # [eV]

print("C_i  =", np.round(C, 4), "GeV")
print("eps_M,i =", eM)

# ---------- 3. 线性系统 A.eps = b (6 复方程, 9 复未知) ----------
# eps 展平: eps[3*i + alpha], i=pair, alpha=(e,mu,tau)
A = np.zeros((6, 9), complex); b = np.zeros(6, complex)

A[0, 3*0+0] = 2*Cv[0]*y[0];  b[0] = UDU[0,0] + Cv[0]*y[0]**2*eM[0]        # (1,1)
A[1, 3*0+1] = Cv[0]*y[0];    A[1, 3*1+0] = Cv[1]*y[1];  b[1] = UDU[0,1]   # (1,2)
A[2, 3*0+2] = Cv[0]*y[0];    A[2, 3*2+0] = Cv[2]*y[2];  b[2] = UDU[0,2]   # (1,3)
A[3, 3*1+1] = 2*Cv[1]*y[1];  b[3] = UDU[1,1] + Cv[1]*y[1]**2*eM[1]        # (2,2)
A[4, 3*1+2] = Cv[1]*y[1];    A[4, 3*2+1] = Cv[2]*y[2];  b[4] = UDU[1,2]   # (2,3)
A[5, 3*2+2] = 2*Cv[2]*y[2];  b[5] = UDU[2,2] + Cv[2]*y[2]**2*eM[2]        # (3,3)

# ---------- 4. 纹理零 -> 6x6 求解 ----------
keep = [k for k in range(9) if (k % 3, k // 3) not in texture_zeros]
x = np.linalg.solve(A[:, keep], b)

eps = np.zeros(9, complex)
for k, idx in enumerate(keep):
    eps[idx] = x[k]

print("\nepsilon 矩阵 x1e10  (行=alpha: e,mu,tau; 列=i: pair 1,2,3):")
print(np.round(eps.reshape(3,3).T * 1e10, 4))

# ---------- 5. 检验: 重建 m_nu 并与 UDU 对比 ----------
mnu = np.zeros((3,3), complex)
for i in range(3):
    Ti = np.zeros((3,3), complex)
    sig = i
    Ti[sig, sig] = y[i]*(2*eps[3*i+sig] - eM[i]*y[i])
    for beta in range(3):
        if beta != sig:
            Ti[sig, beta] = y[i]*eps[3*i+beta]
            Ti[beta, sig] = y[i]*eps[3*i+beta]
    mnu += Cv[i]*Ti

print("\n重建检验 (max|mnu-UDU|):", np.max(np.abs(mnu - UDU)), "eV")

sv = np.linalg.svd(mnu, compute_uv=False)          # 复对称矩阵用奇异值（=质量）
print("奇异值 [eV]:", np.round(np.sort(sv), 6))
print("输入质量 [eV]:", [m1, m2, m3])
print("|eps|/y 最大值:", np.max(np.abs(eps))/np.max(y), "(应 << 1, 微扰自洽)")
