# Course Alignment — mini-differential-geometry

## Nine-School Curriculum Mapping

### MIT 8.962 (General Relativity) — Wald textbook
| Chapter | Topic | Our Implementation |
|---------|-------|-------------------|
| Sec.2.1 | Manifolds | `dg_manifold.h/c` |
| Sec.2.2 | Tangent Vectors | `TangentVector` struct |
| Sec.2.3 | Tensors | `dg_tensor.h/c` |
| Sec.3.1 | Metric, Connection | `dg_metric.h/c`, `dg_connection.h/c` |
| Sec.3.2 | Curvature | `dg_curvature.h/c` |
| Sec.3.3 | Geodesics | `dg_geodesic.h/c` |
| Sec.3.4 | Tetrad Methods | `dg_tetrad.c` |
| Sec.6.3 | Schwarzschild | `geodesic_schwarzschild_init()` |
| App.B | Differential Forms | `dg_forms.h/c` |
| App.C | Lie Derivatives | `dg_lie.h/c` |

### Stanford PHYSICS 230 (General Relativity)
| Topic | Implementation |
|-------|---------------|
| Mathematical Preliminaries | All tensor/manifold operations |
| Einstein Field Equations | `einstein_tensor_compute()` |
| Schwarzschild Solution | Schwarzschild example |
| Cosmology | FRW example |

### Cambridge Part III (General Relativity)
| Topic | Implementation |
|-------|---------------|
| Differential Geometry | Complete implementation |
| Curvature and Bianchi Identities | `bianchi_identity_*_verify()` |
| Tetrad and Spin Connection | `dg_tetrad.c` |
| Exact Solutions | Schwarzschild, FRW examples |

### Princeton PHY 535 (General Relativity)
| Topic | Implementation |
|-------|---------------|
| Riemannian Geometry | Full curvature computation |
| Geodesic Motion | RK4 integration |
| Gravitational Waves | Framework for perturbation theory |

### Caltech Ph 205 (General Relativity)
| Topic | Implementation |
|-------|---------------|
| Tensor Analysis | Complete tensor algebra |
| Differential Forms | `dg_forms.h/c` |
| Curvature Invariants | `kretschmann_scalar()`, `sectional_curvature()` |

### Oxford CMT (Condensed Matter Theory — Geometry)
| Topic | Implementation |
|-------|---------------|
| Topology and Geometry | Manifold structures |
| Differential Forms | Full forms implementation |
| Homology/Cohomology | `form_is_closed()`, `form_is_exact()` |

### ETH 402-0891 (General Relativity)
| Topic | Implementation |
|-------|---------------|
| Mathematical Foundations | All L1-L3 structures |
| Numerical Methods | RK4, matrix inversion, FD derivatives |
| Physical Applications | Schwarzschild, FRW, curvature demo |

### Berkeley/东京大学
| Topic | Implementation |
|-------|---------------|
| Theory Group Methods | Complete toolkit |
| Computational GR | Geodesic integrator |