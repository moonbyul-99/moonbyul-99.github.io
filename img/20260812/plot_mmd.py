# -*- coding: utf-8 -*-
"""
绘制 MMD(P1,Q) 与 MMD(P2,Q) 随核参数 gamma 的变化曲线
对应文章: _posts/2026-08-12-MMDloss.md 的 REMARK 2 反例

- MMD^2(P1,Q) = 1.5 + 0.5*exp(-4γ) - 2*exp(-γ)
- MMD^2(P2,Q) = 2 - 2*exp(-0.25γ)
- γ = 1/(2σ^2)，γ 越小代表核带宽 σ 越大
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'SimSun']
matplotlib.rcParams['axes.unicode_minus'] = False


def mmd_p1_q(gamma):
    """MMD(P1, Q) = sqrt(1.5 + 0.5*exp(-4γ) - 2*exp(-γ))"""
    return np.sqrt(1.5 + 0.5 * np.exp(-4 * gamma) - 2 * np.exp(-gamma))


def mmd_p2_q(gamma):
    """MMD(P2, Q) = sqrt(2 - 2*exp(-0.25γ))"""
    return np.sqrt(2 - 2 * np.exp(-0.25 * gamma))


# gamma 范围: 对数均匀采样, 覆盖从极大带宽(γ->0)到极小带宽(γ->100)
gamma = np.logspace(-3, 2, 8000)
mmd1 = mmd_p1_q(gamma)
mmd2 = mmd_p2_q(gamma)

# 寻找所有交叉点 (两条曲线排名互换的位置, 用线性插值精确定位)
d = mmd1 - mmd2
cross_idx = np.where(np.diff(np.sign(d)) != 0)[0]
cross_points = []
for ci in cross_idx:
    g1, g2 = gamma[ci], gamma[ci + 1]
    d1, d2 = d[ci], d[ci + 1]
    frac = d1 / (d1 - d2)          # 线性插值系数
    gc = g1 + frac * (g2 - g1)
    mc = mmd_p1_q(gc)
    cross_points.append((gc, mc))
print("交叉点 (gamma, MMD):", [(f"{g:.4f}", f"{m:.4f}") for g, m in cross_points])

# 文章第三步中的两个示例: 情况A 小带宽 γ=2, 情况B 大带宽 γ=0.1
g_a, g_b = 2.0, 0.1

fig, ax = plt.subplots(figsize=(11.5, 6.8), dpi=150)

# 依据曲线真实高低关系分段着色背景
bounds = [gamma[0]] + [g for g, _ in cross_points] + [gamma[-1]]
for i in range(len(bounds) - 1):
    g_lo, g_hi = bounds[i], bounds[i + 1]
    gm = np.sqrt(g_lo * g_hi)          # 区间内几何中点(对数刻度)
    color = 'tab:blue' if mmd_p1_q(gm) < mmd_p2_q(gm) else 'tab:orange'
    ax.axvspan(g_lo, g_hi, color=color, alpha=0.07)

# 两条 MMD 曲线
ax.plot(gamma, mmd1, label=r'$MMD(P_1, Q)$', color='tab:blue', lw=2.2)
ax.plot(gamma, mmd2, label=r'$MMD(P_2, Q)$', color='tab:orange', lw=2.2)

# 交叉点标注
cross_colors = ['crimson', 'purple']
for k, (gc, mc) in enumerate(cross_points):
    ax.plot(gc, mc, 'o', ms=6, color=cross_colors[k], zorder=5)
    if k == 0:
        note = (f'交点1 $\\gamma \\approx {gc:.3f}$\n此后 $MMD(P_1) > MMD(P_2)$\n'
                f'($P_2$ 更近, 与直觉一致)')
        xytext = (gc * 2.0, mc + 0.16)
        ha = 'left'
    else:
        note = (f'交点2 $\\gamma \\approx {gc:.3f}$\n此后 $MMD(P_1) < MMD(P_2)$\n'
                f'($\\gamma\\to\\infty$ 逐点匹配区)')
        xytext = (gc * 0.55, mc + 0.16)
        ha = 'right'
    ax.annotate(note, xy=(gc, mc), xytext=xytext, ha=ha,
                arrowprops=dict(arrowstyle='->', color=cross_colors[k], lw=1.2),
                fontsize=9.5,
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                          edgecolor=cross_colors[k], alpha=0.92))

# 两个示例点 (文章数值): 竖线连接两条曲线, 直观展示排名
for g, note, yt in [(g_a, '情况A: 小带宽 $\\gamma=2$\n$MMD(P_1)>MMD(P_2)$\n与直觉一致', 1.28),
                    (g_b, '情况B: 大带宽 $\\gamma=0.1$\n$MMD(P_1)<MMD(P_2)$\n反例!', 0.60)]:
    y1, y2 = mmd_p1_q(g), mmd_p2_q(g)
    ax.scatter([g, g], [y1, y2], s=34, color='darkred', zorder=6)
    ax.vlines(g, ymin=min(y1, y2), ymax=max(y1, y2), color='darkred',
              ls='--', lw=1.1, alpha=0.8, zorder=4)
    ax.text(g, yt, note, ha='center', fontsize=8.5, color='darkred')

# 三个区域的语义说明
ax.text(0.0032, 1.06,
        '① 大带宽 $\\gamma \\to 0$\n$MMD(P_1,Q) < MMD(P_2,Q)$\n'
        'MMD 退化为均值比较,\n被 $P_1$ 的对称性“欺骗”(反例区)',
        fontsize=9, color='tab:blue', va='top',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightblue', alpha=0.55))
ax.text(1.0, 1.06,
        '② 中等带宽\n$MMD(P_1,Q) > MMD(P_2,Q)$\n'
        '与几何直觉 / EMD 一致',
        fontsize=9, color='tab:orange', va='top',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='moccasin', alpha=0.7))
ax.text(30, 1.06,
        '③ 极小带宽 $\\gamma \\to \\infty$\n$MMD(P_1,Q) < MMD(P_2,Q)$\n'
        '逐点匹配, 比较的是\n“概率质量自重叠”',
        fontsize=9, color='tab:blue', va='top',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightblue', alpha=0.55))

# 坐标轴与修饰
ax.set_xscale('log')
ax.set_xlabel(r'核参数 $\gamma$（对数刻度；$\gamma$ 越小 = 高斯核带宽 $\sigma$ 越大）', fontsize=12)
ax.set_ylabel(r'MMD 距离', fontsize=12)
ax.set_title(r'$MMD(P_1,Q)$ 与 $MMD(P_2,Q)$ 随核带宽 $\gamma$ 的变化（P$_1=0.5\delta_1+0.5\delta_{-1}$, '
             r'P$_2=\delta_{0.5}$, Q$=\delta_0$）', fontsize=13)
ax.set_xlim(gamma[0], gamma[-1])
ax.set_ylim(0, 1.45)
ax.grid(True, which='both', ls=':', alpha=0.4)
ax.legend(loc='upper left', fontsize=12, framealpha=0.95)

plt.tight_layout()
# 输出图片到脚本所在目录, 保证从任意工作目录运行都能正确保存
script_dir = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(script_dir, 'mmd_gamma_compare.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"图片已保存: {out}")
