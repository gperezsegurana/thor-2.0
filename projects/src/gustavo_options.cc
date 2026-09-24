// Optimise higher-order-achromat: gustavo.cc with its settings read at run time
// instead of compiled in.
//
//   gustavo_options [--config FILE] [options] LATTICE
//
// --help lists the options. sext_families and oct_families name the families to
// fit; phase = 1 or 2 selects gustavo.cc's two sets of weights. Writes
// flat_file.fit (lattice), b4.out (multipole strengths), final_config.cfg (every
// setting used) and, with flat_file_every = N, flat_file_<k>.fit every N
// iterations.

#include <algorithm>
#include <numeric>
#include <sstream>

// Map order; Lie terms whose derivatives need a higher order are left out.
#define NO 8

// CLI/config parsing and the table of Lie terms are split out to keep this file
// focused on optimization. They come before thor_lib.h, whose min/max macros
// break the standard headers Boost includes.
#include "gustavo_options_config.h"
#include "gustavo_options_terms.h"

#include "thor_lib.h"

// User defined class for parameter dependence.
#include "param_type.h"

int no_tps = NO, ndpt_tps = 5;

// One Lie term: its value and its Jacobian row (one entry per fitted family).
// Lie_gen_scl() multiplies both by cst_scl.
struct Lie_gen_class
{
  std::string label;
  std::vector<int> index_1, index_2;
  double cst_scl = 0e0, cst = 0e0;
  std::vector<double> Jacob;

  void print_label(void) const;
  void print_index_1(void) const;
  void print_index_2(void) const;
  void print_cst_scl(void) const;
  void print_cst(void) const;
  void print_Jacob(void) const;
  void print(void) const;
};

//==============================================================================

// Local re-implementation - because function in tools.cc uses global variables.

ss_vect<tps> get_map(void)
{
  //
  ss_vect<tps> M;

  M.identity();
  M.propagate(1, n_elem);
  return M;
}

void get_ab(const ss_vect<tps> &A1, double alpha[], double beta[], const long int k)
{
  // The map A1 is a globval variable for the function in tools.cc.
  ss_vect<tps> a1, A1_A1tp;

  a1 = A1;
  a1.propagate(1, k);

  A1_A1tp = a1 * tp_S(2, a1);

  alpha[X_] = -h_ijklm(A1_A1tp[x_], 0, 1, 0, 0, 0);
  alpha[Y_] = -h_ijklm(A1_A1tp[y_], 0, 0, 0, 1, 0);
  beta[X_] = h_ijklm(A1_A1tp[x_], 1, 0, 0, 0, 0);
  beta[Y_] = h_ijklm(A1_A1tp[y_], 0, 0, 1, 0, 0);
}

//==============================================================================

void Lie_gen_class::print_label(void) const
{
  printf(" %5s", label.c_str());
}

void Lie_gen_class::print_index_1(void) const
{
  printf(" [");
  for (const int i : index_1)
    printf(" %d", i);
  printf("]");
}

void Lie_gen_class::print_index_2(void) const
{
  printf(" [");
  for (const int i : index_2)
    printf(" %d", i);
  printf("]");
}

void Lie_gen_class::print_cst_scl(void) const
{
  printf(" %6.1e", cst_scl);
}

void Lie_gen_class::print_cst(void) const
{
  printf(" %10.3e", cst);
}

void Lie_gen_class::print_Jacob(void) const
{
  for (const double jac : Jacob)
    printf(" %10.3e", jac);
}

void Lie_gen_class::print(void) const
{
  print_label();
  print_index_1();
  // print_index_2();
  printf("             ");
  print_cst_scl();
  print_cst();
  print_Jacob();
  printf("\n");
}

inline double h_ijklm(const tps &h, const std::vector<int> &ind)
{
  return h_ijklm(h, ind[x_], ind[px_], ind[y_], ind[py_], ind[delta_]);
}

inline double h_ijklm_p(const tps &h, const std::vector<int> &ind)
{
  return h_ijklm_p(h, ind[x_], ind[px_], ind[y_], ind[py_], ind[delta_], 7);
}

Lie_gen_class get_Lie_gen(const std::string &label, const tps &g, const double scl, const int i, const int j, const int k, const int l, const int m)
{
  const std::vector<int> ind = {i, j, k, l, m};

  Lie_gen_class Lie_term;

  Lie_term.label = label;
  Lie_term.index_1 = ind;
  Lie_term.cst_scl = scl;
  // Don't scale yet.
  Lie_term.cst = h_ijklm(g, ind);

  return Lie_term;
}

void prt_system(const param_type &bns)
{
  printf("\n                                   scl.      cst.");
  for (const auto &name : bns.name)
    printf("      %-5s", name.c_str());
  printf("\n                                                 ");
  for (const auto &val : bns.n)
    printf("       %1d   ", val);
  printf("\n                                                 ");
  for (const auto &scl : bns.bnL_scl)
    printf("    %7.1e", scl);
}

Lie_gen_class get_Lie_K_avg_gen(const std::string &label,
                                const tps &K,
                                const double scl,
                                const int i1,
                                const int j1,
                                const int k1,
                                const int l1,
                                const int m1,
                                const int i2,
                                const int j2,
                                const int k2,
                                const int l2,
                                const int m2)
{
  const std::vector<int> ind_1 = {i1, j1, k1, l1, m1}, ind_2 = {i2, j2, k2, l2, m2};

  Lie_gen_class Lie_term;

  auto K_avg_2 = (h_ijklm(K, ind_1) + h_ijklm(K, ind_2)) / 2e0;

  Lie_term.label = label;
  Lie_term.index_1 = ind_1;
  Lie_term.index_2 = ind_2;
  Lie_term.cst_scl = scl;
  // Don't scale yet.
  Lie_term.cst = K_avg_2;

  return Lie_term;
}

void prt_K_avg_term(const Lie_gen_class &term)
{
  term.print_label();
  term.print_index_1();
  term.print_index_2();
  term.print_cst_scl();
  term.print_cst();
  term.print_Jacob();
  printf("\n");
}

// Prints each term under its section heading; defs[k] describes Lie_gen[k].
void prt_terms(const std::vector<const TermDef *> &defs, const std::vector<Lie_gen_class> &Lie_gen)
{
  const char *section = nullptr;
  for (std::size_t k = 0; k < Lie_gen.size(); k++)
  {
    if (!section || std::string(section) != defs[k]->section)
    {
      section = defs[k]->section;
      printf("\n%s\n", section);
    }
    if (defs[k]->source == TermSource::K_avg)
      prt_K_avg_term(Lie_gen[k]);
    else
      Lie_gen[k].print();
  }
}

void compute_Jacob(const param_type &bns, const int prm_no, const ss_vect<tps> &Id_scl, std::vector<Lie_gen_class> &Lie_gen)
{
  const bool debug = false;
  const std::string name = bns.name[prm_no];
  const int
    // Driving terms.
    n = bns.n[prm_no];
  const double bn_scl = (bns.L[prm_no] == 0e0) ? 1e0 : bns.bnL_scl[prm_no] / bns.L[prm_no];

  tps g_re, g_im, K_re, K_im;
  ss_vect<tps> A1, A0, M_res;

  set_bn_par(bns.Fnum[prm_no], n, 7);

  danot_(no_tps - 1);
  auto M = get_map();
  danot_(no_tps);
  auto K = MapNorm(M, g, A1, A0, M_res, 1);
  CtoR(g * Id_scl, g_re, g_im);
  CtoR(K * Id_scl, K_re, K_im);

  if (debug)
    std::cout << std::scientific << std::setprecision(3) << "\n" << std::setw(8) << name << std::setw(4) << bns.Fnum[prm_no] << ":\n" << std::scientific << std::setprecision(3) << K_re << K_im;

  clr_bn_par(bns.Fnum[prm_no], n);
  for (auto &term : Lie_gen)
  {
    const auto &label = term.label;
    const auto &index_1 = term.index_1;
    const auto &index_2 = term.index_2;

    if (label == "g")
      term.Jacob.push_back(bn_scl * h_ijklm_p(g_im, index_1));
    else if (label == "K")
      term.Jacob.push_back(bn_scl * h_ijklm_p(K_re, index_1));
    else if (label == "<K>")
      term.Jacob.push_back(bn_scl * (h_ijklm_p(K_re, index_1) + h_ijklm_p(K_re, index_2)));
    else
    {
      fprintf(stderr, "\ncompute_Jacob - undef. case: %s\n", label.c_str());
      exit(1);
    }
  }
}

inline void compute_Jacob(const ss_vect<tps> &Id_scl, const param_type &bns, std::vector<Lie_gen_class> &Lie_gen)
{
  for (int k = 0; k < bns.n_prm; k++)
    compute_Jacob(bns, k, Id_scl, Lie_gen);
}

void Lie_gen_scl(std::vector<Lie_gen_class> &Lie_gen)
{
  for (auto &term : Lie_gen)
  {
    term.cst = term.cst_scl * term.cst;
    for (auto &jac : term.Jacob)
      jac = term.cst_scl * jac;
  }
}

std::vector<Lie_gen_class> compute_Lie_gen(const ss_vect<tps> &Id_scl, const param_type &bns, const Config &cfg)
{
  double nu[3], xi[2];
  tps g, g_re, g_im, K_re, K_im;
  ss_vect<tps> A1, A0, M_res;
  std::vector<Lie_gen_class> Lie_gen;

  danot_(no_tps - 1);
  auto M = get_map();
  danot_(no_tps);
  auto K = MapNorm(M, g, A1, A0, M_res, 1);
  auto nus = dHdJ(K);
  get_nu_ksi(nus, nu, xi);

  CtoR(g * Id_scl, g_re, g_im);
  CtoR(K * Id_scl, K_re, K_im);

  if (false)
    std::cout << std::scientific << std::setprecision(3) << "\nK:\n" << K_re << "\ng:\n" << g_im;
  if (false)
    printf("\n  nu  = [%6.3f, %6.3f]\n  ksi = [%6.3f, %6.3f]\n", nu[X_], nu[Y_], xi[X_], xi[Y_]);

  std::vector<const TermDef *> defs;
  for (const auto &t : term_table())
  {
    if (!term_active(t, NO))
      continue;
    const auto &i = t.ind_1, &j = t.ind_2;
    const double scl = term_weight(t, cfg);
    defs.push_back(&t);
    switch (t.source)
    {
    case TermSource::g:
      Lie_gen.push_back(get_Lie_gen("g", g_im, scl, i[0], i[1], i[2], i[3], i[4]));
      break;
    case TermSource::K:
      Lie_gen.push_back(get_Lie_gen("K", K_re, scl, i[0], i[1], i[2], i[3], i[4]));
      break;
    case TermSource::K_avg:
      Lie_gen.push_back(get_Lie_K_avg_gen("<K>", K_re, scl, i[0], i[1], i[2], i[3], i[4], j[0], j[1], j[2], j[3], j[4]));
      break;
    }
  }

  compute_Jacob(Id_scl, bns, Lie_gen);
  // Now scale.
  Lie_gen_scl(Lie_gen);

  prt_system(bns);
  prt_terms(defs, Lie_gen);

  return Lie_gen;
}

void get_system(const int m, const int n, const std::vector<Lie_gen_class> &Lie_gen, double **A, double *b)
{
  const bool debug = false;

  auto j1 = 1;
  for (const auto &term : Lie_gen)
  {
    b[j1] = -term.cst;

    int k1 = 1;
    for (const auto &jac : term.Jacob)
    {
      A[j1][k1] = jac;
      ++k1;
    }
    ++j1;
  }

  if (debug)
    dmdump(stdout, (char *)"\nA:", A, m, n, (char *)" %10.3e");
}

std::vector<int> sort_sing_val(const int n, const double w[])
{
  const bool debug = false;

  std::vector<int> index(n);
  std::iota(index.begin(), index.end(), 1); // fill with 1..n

  std::sort(index.begin(), index.end(), [&w](int a, int b) { return w[a] > w[b]; });

  if (debug)
  {
    std::cout << "\n";
    for (const auto &idx : index)
      std::cout << " " << idx;
    std::cout << "\n";
  }

  return index;
}

void get_sing_val(const int n, double w[], const int svd_n_cut, const Config &cfg)
{
  std::vector<int> ind;

  ind = sort_sing_val(n, w);
  printf("\ninitial:\n");
  for (int k = 1; k <= n; k++)
  {
    printf("  %9.3e", w[ind[k - 1]]);
    if (k <= n - svd_n_cut)
      w[ind[k - 1]] *= cfg.scl_svd[k - 1];
    else
    {
      w[ind[k - 1]] = 0e0;
      printf(" (zeroed)");
    }
  }
  printf("\nscaled & singular values removed\n");
  for (int k = 1; k <= n; k++)
    printf("  %9.3e", w[ind[k - 1]]);
  printf("\n");
}

void set_Fam(param_type &bns, const int k, const double scl, const double *dbnL)
{
  double bnL_ext;

  bnL_ext = get_bnL(bns.Fnum[k], 1, bns.n[k]) + scl * bns.bnL_scl[k] * dbnL[k + 1];
  bns.bnL[k] = bnL_internal(bnL_ext, bns.bnL_min[k], bns.bnL_max[k]);
  set_bnL(bns.Fnum[k], bns.n[k], bnL_bounded(bns.bnL[k], bns.bnL_min[k], bns.bnL_max[k]));
}

void set_bnL(const double scl, const double *dbnL, param_type &bns)
{
  for (std::size_t k = 0; k < static_cast<std::size_t>(bns.n_prm); ++k)
    set_Fam(bns, k, scl, dbnL);
}

void prt_bend(FILE *outf, const int loc, const int n)
{
  const elem_type<double> *elemp = &elem[loc - 1];

  fprintf(outf,
          "%-8s: Multipole, L = %7.5f, Phi = %7.5f, Phi_1 = %7.5f"
          ", Phi_2 = %7.5f,\n"
          "          HOM = (%d, %12.5e, 0e0, %d, %12.5e, 0e0),\n"
          "          N = n_bend;\n",
          elemp->Name,
          elemp->L,
          elemp->L * elemp->mpole->h_bend * 180e0 / M_PI,
          elemp->mpole->edge1,
          elemp->mpole->edge2,
          Quad,
          get_bn(elem[loc - 1].Fnum, elem[loc - 1].Knum, Quad),
          n,
          get_bn(elem[loc - 1].Fnum, elem[loc - 1].Knum, n));
}

void prt_quad(FILE *outf, const int loc, const int n)
{
  fprintf(outf,
          "%-8s: Multipole, L = %7.5f,\n          HOM = (%d, %12.5e, 0e0,"
          " %d, %12.5e, 0e0),\n          N = n_quad;\n",
          elem[loc - 1].Name,
          elem[loc - 1].L,
          Quad,
          get_bn(elem[loc - 1].Fnum, elem[loc - 1].Knum, Quad),
          n,
          get_bn(elem[loc - 1].Fnum, elem[loc - 1].Knum, n));
}

void prt_single_mult(FILE *outf, const int loc, const int n)
{
  switch (n)
  {
  case Sext:
    fprintf(outf, "%-8s: Sextupole, L = %7.5f, B_3 = %12.5e, N = n_sext;\n", elem[loc - 1].Name, elem[loc - 1].L, elem[loc - 1].mpole->bn[Sext - 1]);
    break;
  case Oct:
    fprintf(outf, "%-8s: Octupole, L = %7.5f, B_4 = %12.5e, N = %d;\n", elem[loc - 1].Name, elem[loc - 1].L, elem[loc - 1].mpole->bn[Oct - 1], elem[loc - 1].mpole->n_step);
    break;
  default:
    printf("\nprt_single_mult - undefined multipole order: %d\n", n);
    break;
  }
}

int get_n_mpole(const int loc)
{
  int n_mpole = 0;

  for (int k = 0; k < elem[loc - 1].mpole->order; k++)
  {
    if ((elem[loc - 1].mpole->bn[k] != 0e0) || (elem[loc - 1].mpole->an[k] != 0e0))
      n_mpole++;
  }

  return n_mpole;
}

void prt_mult(FILE *outf, const int loc, const int n)
{
  std::string name;
  bool first = true;
  int n_step;
  double L;

  if (get_n_mpole(loc) == 1)
    prt_single_mult(outf, loc, n);
  else
  {
    name = elem[loc - 1].Name;
    L = elem[loc - 1].L;
    fprintf(outf, "%-8s: multipole, l = %7.5f, hom = (", name.c_str(), L);
    for (int k = 0; k < elem[loc - 1].mpole->order; k++)
    {
      if ((elem[loc - 1].mpole->bn[k] != 0e0) || (elem[loc - 1].mpole->an[k] != 0e0))
      {
        if (first)
        {
          fprintf(outf, "\n            %d, %12.5e, %12.5e", k + 1, elem[loc - 1].mpole->bn[k], elem[loc - 1].mpole->an[k]);
          first = false;
        }
        else
          fprintf(outf, ",\n            %d, %12.5e, %12.5e", k + 1, elem[loc - 1].mpole->bn[k], elem[loc - 1].mpole->an[k]);
      }
    }
    n_step = elem[loc - 1].mpole->n_step;
    fprintf(outf, "), n = %d;\n", n_step);
  }
}

void prt_bn(const param_type &bns)
{
  const std::string file_name = "b4.out";

  FILE *outf;

  outf = file_write(file_name.c_str());

  fprintf(outf, "\n");
  for (int k = 0; k < bns.n_prm; k++)
  {
    auto loc = (bns.Fnum[k] > 0) ? get_loc(bns.Fnum[k], 1) : bns.locs[k][0];
    if (elem[loc - 1].mpole->n_design == Dip)
      prt_bend(outf, loc, bns.n[k]);
    else if (elem[loc - 1].mpole->n_design == Quad)
      prt_quad(outf, loc, bns.n[k]);
    else
      prt_mult(outf, loc, bns.n[k]);
  }

  fclose(outf);
}

void correct(param_type &bns, const std::vector<Lie_gen_class> &Lie_gen, const int svd_n_cut, const double scl, const Config &cfg)
{
  const int m = Lie_gen.size(), n = bns.n_prm;

  double **A, **U, **V, *w, *b, *dbnL;

  printf("\nsvd:\n  m = %d n = %d\n", m, n);

  A = dmatrix(1, m, 1, n);
  U = dmatrix(1, m, 1, n);
  V = dmatrix(1, n, 1, n);
  w = dvector(1, n);
  b = dvector(1, m);
  dbnL = dvector(1, n);

  get_system(m, n, Lie_gen, A, b);

  dmcopy(A, m, n, U);
  dsvdcmp(U, m, n, w, V);
  get_sing_val(n, w, svd_n_cut, cfg);

  dsvbksb(U, w, V, m, n, b, dbnL);

  set_bnL(scl, dbnL, bns);
  bns.print();
  prt_bn(bns);

  free_dmatrix(A, 1, m, 1, n);
  free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n);
  free_dvector(w, 1, n);
  free_dvector(b, 1, m);
  free_dvector(dbnL, 1, n);
}

void no_mpoles(const int n)
{
  printf("\nzeroing multipoles: %d\n", n);
  for (auto j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}

// True if the lattice defines the family; get_Fnum() exits with status 0 instead.
bool family_defined(const std::string &name)
{
  for (int k = 0; k < max_Family; k++)
    if (name == Families[k].Name)
      return true;
  return false;
}

void get_bns(param_type &bns, const Config &cfg)
{
  for (const auto &name : cfg.sext_families)
    if (!family_defined(name))
    {
      fprintf(stderr, "\nget_bns - sextupole family not in the lattice: %s\n", name.c_str());
      exit(1);
    }
  for (const auto &name : cfg.oct_families)
    if (!family_defined(name))
    {
      fprintf(stderr, "\nget_bns - octupole family not in the lattice: %s\n", name.c_str());
      exit(1);
    }

  if (cfg.b_3_zero)
    no_mpoles(Sext);
  if (cfg.b_4_zero)
    no_mpoles(Oct);

  if (cfg.b_3_opt)
    for (const auto &name : cfg.sext_families)
      bns.add_Fam(name, Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
  if (cfg.b_4_opt)
    for (const auto &name : cfg.oct_families)
      bns.add_Fam(name, Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);

  if (bns.Fnum.empty())
  {
    fprintf(stderr, "\nget_bns - no families to fit: set sext_families/oct_families and b3_opt/b4_opt\n");
    exit(1);
  }
}

void chk_lat(void)
{
  double alpha[2], beta[2], nu[3], ksi[2];
  tps g, K;
  ss_vect<tps> nus, M, A1, A0, M_res;

  danot_(2);
  M = get_map();
  danot_(3);
  K = MapNorm(M, g, A1, A0, M_res, 1);
  nus = dHdJ(K);
  get_nu_ksi(nus, nu, ksi);
  get_ab(A1, alpha, beta, 0);
  printf("\n  alpha = [%6.3f, %6.3f]\n  beta  = [%6.3f, %6.3f]\n"
         "  nu    = [%6.3f, %6.3f]\n  ksi   = [%6.3f, %6.3f]\n",
         alpha[X_],
         alpha[Y_],
         beta[X_],
         beta[Y_],
         nu[X_],
         nu[Y_],
         ksi[X_],
         ksi[Y_]);
}

void set_state(void)
{
  rad_on = false;
  H_exact = false;
  totpath_on = false;
  cavity_on = false;
  quad_fringe_on = false;
  emittance_on = false;
  IBS_on = false;
}

int main(int argc, char *argv[])
{
  double twoJ;
  ss_vect<tps> Id_scl;
  param_type bns;
  std::vector<Lie_gen_class> Lie_gen;

  set_state();

  // Parse command-line arguments
  Config cfg = parse_arguments(argc, argv);

  rd_mfile(cfg.lattice_file.c_str(), elem);
  rd_mfile(cfg.lattice_file.c_str(), elem_tps);

  if (true)
  {
    printf("\nLattice file: %s\n", cfg.lattice_file.c_str());
    printf("Sextupole families: %s\n", join_names(cfg.sext_families).c_str());
    printf("Octupole families: %s\n", join_names(cfg.oct_families).c_str());
    printf("Phase: %d\n", cfg.phase);
    printf("SVD n_cut: %d\n", cfg.svd_n_cut);
    printf("Max iterations: %d\n", cfg.max_iter);
    printf("Step: %.6f\n", cfg.step);
    printf("Intermediate flat file frequency: %d\n", cfg.flat_file_every);
  }
  // Initialize the symplectic integrator after the energy has been defined.
  ini_si();
  // Disable log messages from TPSALib & LieLib.
  idprset(-1);

  daeps_(1e-30);

  chk_lat();

  Id_scl.identity();
  for (int k = 0; k < 2; k++)
  {
    twoJ = sqr(cfg.A_max[k]) / cfg.beta_inj[k];
    Id_scl[2 * k] *= sqrt(twoJ);
    Id_scl[2 * k + 1] *= sqrt(twoJ);
  }
  Id_scl[delta_] *= cfg.delta_max;

  get_bns(bns, cfg);
  resolve_scl_svd(cfg, std::count(bns.n.begin(), bns.n.end(), Sext), std::count(bns.n.begin(), bns.n.end(), Oct));
  write_config_to_file(cfg, "final_config.cfg");
  bns.ini_prm();
  bns.print();

  printf("\n");
  for (int k = 1; k <= cfg.max_iter; k++)
  {
    printf("\nk = %d:", k);
    Lie_gen = compute_Lie_gen(Id_scl, bns, cfg);
    correct(bns, Lie_gen, cfg.svd_n_cut, cfg.step, cfg);

    prt_mfile("flat_file.fit");
    if (cfg.flat_file_every > 0 && (k % cfg.flat_file_every == 0))
    {
      std::ostringstream file_name;
      file_name << "flat_file_" << k << ".fit";
      prt_mfile(file_name.str().c_str());
    }
  }
  Lie_gen = compute_Lie_gen(Id_scl, bns, cfg);
  return 0;
}
