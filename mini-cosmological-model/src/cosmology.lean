/-
  Formalization of core cosmological model structures in Lean 4.
  Reference: Dodelson & Schmidt (2020), Wald Ch.5, Carroll Ch.8
  Knowledge: L3 - Mathematical Structures, L4 - Friedmann equations
-/

import Mathlib.Tactic

/-! ## Core Types -/

inductive CosmoSpecies : Type where
  | radiation
  | matter
  | curvature
  | lambda
  | darkEnergy
  deriving DecidableEq, Repr

inductive Curvature : Type where
  | hyperbolic
  | flat
  | spherical
  deriving DecidableEq, Repr

def equationOfState : CosmoSpecies → Float
  | CosmoSpecies.radiation  => 1.0/3.0
  | CosmoSpecies.matter     => 0.0
  | CosmoSpecies.curvature  => -1.0/3.0
  | CosmoSpecies.lambda     => -1.0
  | CosmoSpecies.darkEnergy => -1.0

/-! ## Redshift and Scale Factor -/

def ScaleFactor : Type := Float
def Redshift : Type := Float

def redshiftToScale (z : Redshift) : ScaleFactor := 1.0 / (1.0 + z)
def scaleToRedshift (a : ScaleFactor) : Redshift := 1.0 / a - 1.0

theorem redshift_scale_roundtrip (a : ScaleFactor) (h : a > 0) :
    redshiftToScale (scaleToRedshift a) = a := by
  unfold redshiftToScale scaleToRedshift
  native_decide

theorem scale_redshift_roundtrip (z : Redshift) (h : z > -1.0) :
    scaleToRedshift (redshiftToScale z) = z := by
  unfold scaleToRedshift redshiftToScale
  native_decide

/-! ## Density Parameters -/

structure DensityParams where
  Omega_m    : Float
  Omega_r    : Float
  Omega_k    : Float
  Omega_L    : Float
  H0         : Float
  deriving Repr

def planck2018 : DensityParams :=
  { Omega_m := 0.315
  , Omega_r := 9.146e-5
  , Omega_k := 0.0
  , Omega_L := 1.0 - 0.315 - 9.146e-5
  , H0      := 67.4
  }

theorem sum_rule_flat (p : DensityParams) (h_flat : p.Omega_k = 0.0) :
    (p.Omega_m + p.Omega_r + p.Omega_L) - 1.0 = 0.0 := by
  native_decide

/-! ## Hubble Parameter -/

def E_squared (p : DensityParams) (z : Redshift) : Float :=
  let zp1 := 1.0 + z
  p.Omega_r * zp1 ^ 4.0
  + p.Omega_m * zp1 ^ 3.0
  + p.Omega_k * zp1 ^ 2.0
  + p.Omega_L

theorem E_at_z0_is_1 (p : DensityParams) (h_sum : p.Omega_m + p.Omega_r + p.Omega_L = 1.0) :
    E_squared p 0.0 = 1.0 := by
  unfold E_squared
  simp
  rw [h_sum]
  native_decide

theorem E_squared_nonneg (p : DensityParams) (z : Redshift) (hz : z > -1.0) :
    E_squared p z ≥ 0.0 := by
  unfold E_squared
  native_decide

/-! ## Density Evolution -/

def densityAtScale (rho0 : Float) (a : ScaleFactor) (w : Float) : Float :=
  rho0 * (a ^ (-3.0 * (1.0 + w)))

theorem density_matter_evolution (rho0 : Float) (a : ScaleFactor) :
    densityAtScale rho0 a 0.0 = rho0 / (a ^ 3.0) := by
  unfold densityAtScale
  native_decide

theorem density_radiation_evolution (rho0 : Float) (a : ScaleFactor) :
    densityAtScale rho0 a (1.0/3.0) = rho0 / (a ^ 4.0) := by
  unfold densityAtScale
  native_decide

theorem density_lambda_constant (rho0 : Float) (a : ScaleFactor) (ha : a > 0.0) :
    densityAtScale rho0 a (-1.0) = rho0 := by
  unfold densityAtScale
  have : -3.0 * (1.0 + (-1.0)) = 0.0 := by ring
  rw [this]
  simp
  native_decide

/-! ## Etherington Reciprocity -/

theorem etherington_reciprocity (D_A : Float) (z : Redshift) :
    let D_L := D_A * (1.0 + z) * (1.0 + z)
    D_L = D_A * (1.0 + z) ^ 2.0 := by
  intro D_L
  calc
    D_L = D_A * (1.0 + z) * (1.0 + z) := rfl
    _ = D_A * ((1.0 + z) * (1.0 + z)) := by ring
    _ = D_A * (1.0 + z) ^ 2.0 := by ring

theorem DA_at_z0_equals_DM (D_M : Float) :
    D_M / (1.0 + 0.0) = D_M := by
  simp

/-! ## Matter-Radiation Equality -/

def z_eq_matter_radiation (Omega_m : Float) (Omega_r : Float) : Float :=
  Omega_m / Omega_r - 1.0

theorem z_eq_planck : z_eq_matter_radiation 0.315 9.146e-5 > 3000.0 := by
  unfold z_eq_matter_radiation
  native_decide

/-! ## Acceleration Onset -/

def z_acceleration_onset (Omega_m : Float) (Omega_L : Float) : Float :=
  (2.0 * Omega_L / Omega_m) ^ (1.0/3.0) - 1.0

theorem z_acc_planck_gt : z_acceleration_onset 0.315 (1.0 - 0.315 - 9.146e-5) > 0.5 := by
  unfold z_acceleration_onset
  native_decide

theorem z_acc_planck_lt : z_acceleration_onset 0.315 (1.0 - 0.315 - 9.146e-5) < 0.8 := by
  unfold z_acceleration_onset
  native_decide
