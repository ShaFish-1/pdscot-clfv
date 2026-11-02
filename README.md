# PDSCOT — 准 Dirac Scotogenic 模型代码

复现与独立验证 **arXiv:2605.05312**（*Monochromatic neutrinos from scotogenic dark matter*）
所用模型与计算代码：三对准 Dirac scotogenic 暗物质模型（LanHEP），以及 CLFV 过程验证脚本。

## model/ — 模型

| 文件 | 说明 |
|---|---|
| `ScotogenicDiracExt.mdl` | LanHEP 模型：SM + 惰性双态 + 3 对 (N, Nᶜ)，Feynman 规范 |
| `SMmodern.mdl` | 内联标准模型（供 `use SMmodern.` 载入） |
| `PDSCOT_singlepair.fr` | FeynRules 单对模型（质量本征基，供 FeynArts/FeynCalc 使用） |
| `input_bench.par` | 基准输入参数（M_i、ΔM_i、y、ε、m_R、m_I、m_η⁺、λ₂、μ_η²） |
| `patch_prtcls_short.py` | 粒子名短名修片脚本（micrOMEGAs 的 3 字符限制兼容） |

模型关键结构（约定）：
- 质量矩阵 `{{−μ, M}, {M, −μ}}` → 本征态 (n, p)，45° 转动，`M_n = M−μ`、`M_p = M+μ`；
- 耦合 H 矩阵：`H[α,2i−1] = (y+ε)/√2`，`H[α,2i] = i(y−ε)/√2`；
- 标量势含 λ₂, λ₃, λ₄, λ₅；导出 `λ₅v² = m_R²−m_I²`、`λ₄v² = m_R²+m_I²−2m_η⁺²`。

## scripts/ — 计算脚本

| 文件 | 说明 |
|---|---|
| `lfv_verify_S1S4_full.py` | S1–S4 四个参数点的全部 CLFV 过程验证：A_D/A_ND（含相位）、BR(ℓ→ℓ′γ)、BR(μ→3e)、CR(μ–e)、κ_D、LNC/LNV。与参考值全部 ratio = 1.000 |
| `kappa_check.py` | 圈函数对数导数 κ = d lnΛ/d lnM（论文 Eq.(10) 的 ΔM 项修正量）：近阈 κ → 1/3 |
| `solve_eps.py` | 由中微子振荡数据解 ε 矩阵（PMNS → 模型侧 T 矩阵 → 线性方程组） |
| `lfv_observables.py` | 观测量计算（BR/CR/三体衰变） |

## 运行

```bash
# 模型编译（LanHEP → CalcHEP）
lhep -ca ScotogenicDiracExt.mdl
cd work && make main=main.c
# 微波计算（micrOMEGAs）
./main input_bench.par

# CLFV 验证（纯 Python，无需额外依赖；mpmath 可选）
python3 scripts/lfv_verify_S1S4_full.py
python3 scripts/kappa_check.py
```

## 验证结果摘要

- **S1–S4 全闭合**：四个点的振幅（幅度+相位）、三个辐射道 BR、BR(μ→3e)、
  κ_D、LNC/LNV 与参考文件逐位一致（ratio = 1.000）；
- **论文 Eq.(10) 修正**：ΔM 项应乘 κ = 1 + M G′/G——近阈 κ ≈ 1/3（论文项大 ~3 倍）；
- **论文公式四处问题**：Eq.(11) 因子 2、Eq.(28) 排版残缺、Eq.(33) 丢 δF₂、Eq.(10) κ 近似。

## 约定与致谢

- 振幅约定：Toma–Vicente, JHEP 01 (2014) 160 [arXiv:1312.2840] 附录 A；
- S1–S4 输入参数与参考数值来自合作者项目的约定文件，验证脚本内嵌以做独立复算。

## 许可

MIT
