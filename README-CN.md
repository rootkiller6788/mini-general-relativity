# Mini General Relativity（迷你广义相对论）

一套**从零开始、零依赖的 C 语言实现**，涵盖爱因斯坦广义相对论的核心主题——从微分几何基础到黑洞热力学、引力波和数值相对论。每个子模块对应 MIT 8.962（广义相对论）及其他顶尖物理课程，将教科书中的方程转化为可运行、可读的 C 代码。

## 子模块

| 子模块 | 主题 | 对应课程 |
|------|------|---------|
| [mini-differential-geometry](mini-differential-geometry/) | 仿射联络、协变导数、Riemann/Ricci/Einstein/Weyl 张量、微分形式、外微积分、Hodge 对偶、Lie 导数、测地线 | MIT 8.962, Cambridge Part III, ETH 402-0891 |
| [mini-schwarzschild](mini-schwarzschild/) | Schwarzschild 度规、事件视界、测地线运动（ISCO、光子球）、引力时间膨胀、红移、光线偏折、近日点进动、Shapiro 延迟 | MIT 8.962, Wald Ch.6, Carroll Ch.5 |
| [mini-kerr-metric](mini-kerr-metric/) | Kerr 度规（Boyer-Lindquist 和 Kerr-Schild 坐标）、视界、能层、环奇点、Penrose 过程、测地线、超辐射 | MIT 8.962, Wald Ch.7, Kerr (1963) |
| [mini-einstein-equations](mini-einstein-equations/) | Einstein 场方程、Einstein-Hilbert 作用量、应力-能量张量、曲率、测地线、坐标系、FLRW 宇宙学 | MIT 8.962, Weinberg |
| [mini-gravitational-waves](mini-gravitational-waves/) | 线性化引力、TT 规范、四极矩公式、致密双星旋近/并合、匹配滤波、信噪比、脉冲星计时、引力波记忆 | MIT 8.962, Maggiore (2008) |
| [mini-black-holes](mini-black-holes/) | 黑洞度规（Schwarzschild/Kerr/RN/Kerr-Newman）、热力学、Hawking 辐射、Bekenstein-Hawking 熵、黑洞力学四定律、信息悖论、Penrose 过程、超辐射 | MIT 8.962, Hawking & Ellis (1973) |
| [mini-cosmological-model](mini-cosmological-model/) | FLRW 度规、Friedmann 方程、尺度因子演化、宇宙学距离、功率谱、宇宙学微扰 | MIT 8.962, Dodelson, Carroll Ch.8 |
| [mini-numerical-relativity](mini-numerical-relativity/) | ADM 3+1 分解、BSSN 公式化、3D 笛卡尔网格、外曲率、视界定位、初始数据 | MIT 8.962, Baumgarte & Shapiro (2010), Gourgoulhon (2012) |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅依赖 `libc` 和 `libm`
- **模块自包含** — 每个目录拥有独立的 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，记录与课程的对应关系
- **物理优先设计** — 张量运算、测地线积分和度规计算直接对应 GR 教科书（Wald、Carroll、MTW）中的数学形式

## 构建方式

每个模块独立构建。进入子模块目录运行：

```bash
cd mini-differential-geometry
make all    # 构建所有目标
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-general-relativity/
├── mini-differential-geometry/  # 仿射联络、曲率张量、微分形式、Lie 导数
├── mini-schwarzschild/          # Schwarzschild 度规、测地线、时间膨胀、光线偏折
├── mini-kerr-metric/            # Kerr 度规、视界、能层、Penrose 过程、超辐射
├── mini-einstein-equations/     # Einstein 场方程、Einstein-Hilbert 作用量、应力-能量张量
├── mini-gravitational-waves/    # 线性化引力、TT 规范、双星旋近、匹配滤波
├── mini-black-holes/            # 黑洞热力学、Hawking 辐射、熵、信息悖论
├── mini-cosmological-model/     # FLRW 度规、Friedmann 方程、尺度因子、微扰
└── mini-numerical-relativity/   # ADM 分解、BSSN 公式化、3D 网格、视界定位

```

## 许可证

MIT
