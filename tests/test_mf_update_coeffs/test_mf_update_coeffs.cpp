// test_mf_update_coeffs.cpp — 1000-step test for MFUpdateCoeffs2D operator
// Verifies: mu(phi), chi(phi), omega_psi(phi) linear interpolation correctness
#include "freelb.h"
#include "freelb.hh"

using T = double;

// Hardcoded params matching Guo2025 ferrofluid case
constexpr int Ni = 4, Nj = 4;
constexpr T Cell_Len = 1.0;
constexpr int BlockCellLen = 4;
constexpr T InterfaceWidth = 4.0;
constexpr T Mobility = 0.01;

// Magnetic parameters (simple values for testing)
constexpr T mu_l  = 1.0;
constexpr T mu_h  = 5.0;
constexpr T chi_l = 0.0;
constexpr T chi_h = 1.0;
constexpr T H0    = 0.01;
constexpr T PsiSolver_K = 3.0;  // 1/cs², matching the legacy omega expectation

// Expected values at phi = 0.5
constexpr T phi_test  = 0.5;
constexpr T mu_exp    = mu_l + phi_test * (mu_h - mu_l);          // 3.0
constexpr T chi_exp   = chi_l + phi_test * (chi_h - chi_l);       // 0.5

// D2Q5 cs² = 1/3; nu_psi = mu (diffusion coeff = mu); tau = 0.5 + mu/cs² = 0.5 + 3*mu
// omega = 1/tau
constexpr T cs2 = 1.0/3.0;
constexpr T tau_exp  = 0.5 + mu_exp / cs2;                        // 0.5 + 3*3.0 = 9.5
constexpr T omega_exp = 1.0 / tau_exp;                             // ≈ 0.105263
constexpr T tol = 1e-10;

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);

  // ------------------ minimal geometry (4x4, single block) ------------------
  AABB<T, 2> domain(Vector<T, 2>{T(0), T(0)},
                    Vector<T, 2>{T(Ni * Cell_Len), T(Nj * Cell_Len)});
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, 1);
  GeoHelper.AdaptiveOptimization(1);
  GeoHelper.LoadBalancing(1);
  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ flag field ------------------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });

  // ------------------ PF lattice (minimal: just PHI) ------------------
  BaseConverter<T> BaseConv(D2Q9<T>::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.01), T(1.0));
  BaseConverter<T> PFBaseConv(D2Q9<T>::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), T(1.0));

  // PF: minimal — just PHI (only field needed by MFUpdateCoeffs2D)
  using PFFIELDS = TypePack<PHI<T>>;
  ValuePack PFInit(T{0.5});
  using PFCELL = Cell<T, D2Q9<T>, PFFIELDS>;
  BlockLatticeManager<T, D2Q9<T>, PFFIELDS> PFLattice(Geo, PFInit, PFBaseConv);

  // MF: magnetic fields — without POP (D2Q5 distribution not needed yet)
  using MFFIELDS = typename mfield::MFFIELDS<T>;
  using MFFIELDS_NOPOP = MFFIELDS;  // no POP needed for coeff update test
  using MFREF = TypePack<PHI<T>>;
  using MFPACK = TypePack<MFFIELDS_NOPOP, MFREF>;
  using MF_FULL = typename ExtractFieldPack<MFPACK>::mergedpack;
  using MFCELL = Cell<T, D2Q5<T>, MF_FULL>;

  BaseConverter<T> MFBaseConv(D2Q5<T>::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), T(1.0));
  ValuePack MFInit(T{}, T{}, T{}, T{}, T{}, T{}, T{},
                   mu_l, mu_h, chi_l, chi_h, H0, PsiSolver_K);
  BlockLatticeManager<T, D2Q5<T>, MFPACK> MFLattice(
    Geo, MFInit, MFBaseConv, &PFLattice.getField<PHI<T>>());

  // Broadcast block-level params
  mfield::BroadcastAllMFParams<T>(MFLattice, const_cast<T&>(mu_l), const_cast<T&>(mu_h),
                                  const_cast<T&>(chi_l), const_cast<T&>(chi_h),
                                  const_cast<T&>(H0), const_cast<T&>(PsiSolver_K));

  // ------------------ init phi = 0.5 everywhere ------------------
  PFLattice.getField<PHI<T>>().InitValue(T{0.5});

  // ------------------ coupling: PF → MF (MFUpdateCoeffs2D) ------------------
  using CoeffTask =
    tmp::Key_TypePair<BulkFlag,
      mfield::MFUpdateCoeffs2D<PFCELL, MFCELL>>;
  using CoeffSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, MFCELL, CoeffTask>;
  BlockLatManagerCoupling CoeffCoupling(PFLattice, MFLattice);

  // ------------------ run 1000 steps ------------------
  Timer timer;
  for (int step = 0; step < 1000; ++step) {
    CoeffCoupling.ApplyInnerCellDynamics<CoeffSelector>(timer(), FlagFM);
    ++timer;
  }

  // ------------------ verify at a test cell ------------------
  const auto& block = Geo.getBlock(0);
  const auto& proj = block.getProjection();
  int overlap = 0;
  std::size_t id = (1) * proj[1] + (1);  // interior cell at (1,1)

  // Access MF cell
  auto& mfBlockLat = MFLattice.getBlockLat(0);
  MFCELL mf_cell(id, mfBlockLat);

  T mu_val   = mf_cell.template get<mfield::MU_PERCELL<T>>();
  T chi_val  = mf_cell.template get<mfield::CHI_PERCELL<T>>();
  T omega_val = mf_cell.template get<mfield::OMEGA_PSI<T>>();

  bool pass = true;

  if (std::abs(mu_val - mu_exp) > tol) {
    std::cout << "FAIL mu: expected " << mu_exp << " got " << mu_val
              << " (err=" << std::abs(mu_val-mu_exp) << ")" << std::endl;
    pass = false;
  }
  if (std::abs(chi_val - chi_exp) > tol) {
    std::cout << "FAIL chi: expected " << chi_exp << " got " << chi_val
              << " (err=" << std::abs(chi_val-chi_exp) << ")" << std::endl;
    pass = false;
  }
  if (std::abs(omega_val - omega_exp) > tol) {
    std::cout << "FAIL omega: expected " << omega_exp << " got " << omega_val
              << " (err=" << std::abs(omega_val-omega_exp) << ")" << std::endl;
    pass = false;
  }

  // Also check omega is clamped (should NOT be out of range)
  if (omega_val < 0.01 || omega_val > 1.95) {
    std::cout << "FAIL omega out of range: " << omega_val << std::endl;
    pass = false;
  }

  // Bonus check: change phi to 1.0, run one more step, verify values change
  PFLattice.getField<PHI<T>>().InitValue(T{1.0});
  CoeffCoupling.ApplyInnerCellDynamics<CoeffSelector>(timer(), FlagFM);
  ++timer;

  MFCELL mf_cell2(id, mfBlockLat);
  T mu2 = mf_cell2.template get<mfield::MU_PERCELL<T>>();
  T chi2 = mf_cell2.template get<mfield::CHI_PERCELL<T>>();
  if (std::abs(mu2 - mu_h) > tol) {
    std::cout << "FAIL dynamic check: expected mu_h=" << mu_h << " got " << mu2 << std::endl;
    pass = false;
  }
  if (std::abs(chi2 - chi_h) > tol) {
    std::cout << "FAIL dynamic check: expected chi_h=" << chi_h << " got " << chi2 << std::endl;
    pass = false;
  }

  if (pass) {
    std::cout << "PASS: 1000 steps, phi=0.5 -> mu=" << mu_val << " chi=" << chi_val
              << " omega=" << omega_val << std::endl;
    std::cout << "PASS: dynamic check, phi=1.0 -> mu=" << mu2 << " chi=" << chi2 << std::endl;
  }

  return pass ? 0 : 1;
}
