// Optimise higher-order-achromat:
//
// - Runtime controls: the program accepts `--phase=N`, `--svd_n_cut=N`, and
//   `--max_iter=N` (optional) and expects a lattice file argument. `svd_n_cut`
//   controls how many singular values are removed (zeroed) during SVD-based
//   solving. `max_iter` sets the maximum number of optimization iterations
//   (default 300).
// - Outputs: writes `flat_file.fit` (lattice) and `b4.out` (multipole strenghts).

#include <algorithm>
#include <boost/any.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cassert>
#include <map>
#include <numeric>
#include <sstream>
#include <typeinfo>

#define ASSERT_MSG(cond, msg)                                                                                                                                                                          \
  do                                                                                                                                                                                                   \
  {                                                                                                                                                                                                    \
    if (!(cond))                                                                                                                                                                                       \
    {                                                                                                                                                                                                  \
      fprintf(stderr, "Assertion failed: %s\n", msg);                                                                                                                                                  \
      std::abort();                                                                                                                                                                                    \
    }                                                                                                                                                                                                  \
  } while (0)

#define NO 8

#include "thor_lib.h"

// User defined class for parameter dependence.
#include "param_type.h"

// CLI/config parsing is split out to keep this file focused on optimization.
#include "gustavo_options_config.h"

int no_tps = NO, ndpt_tps = 5;

extern double b2_max;

class Lie_gen_class
{
private:
  std::map<std::string, boost::any> g;

  template <typename T> T &get_ref(const std::string &key)
  {
    return boost::any_cast<T &>(g[key]);
  }

  template <typename T> const T &get_ref(const std::string &key) const
  {
    return boost::any_cast<const T &>(g.at(key));
  }

public:
  Lie_gen_class()
      : g({
          {"label", std::string{}},
          {"index_1", std::vector<int>{}},
          {"index_2", std::vector<int>{}},
          {"cst_scl", double{0}},
          {"cst", double{0}},
          {"Jacob", std::vector<double>{}},
        })
  {
  }

  template <typename T> T get(const std::string &key) const
  {
    return boost::any_cast<T>(g.at(key));
  }

  template <typename T> void set(const std::string &key, const T &value)
  {
    g[key] = value;
  }

  std::string get_label(void) const
  {
    return get<std::string>("label");
  }

  std::vector<int> get_index_1(void) const
  {
    return get<std::vector<int>>("index_1");
  }

  std::vector<int> get_index_2(void) const
  {
    return get<std::vector<int>>("index_2");
  }

  double get_cst_scl(void) const
  {
    return get<double>("cst_scl");
  }

  double get_cst(void) const
  {
    return get<double>("cst");
  }

  std::vector<double> get_Jacob(void) const
  {
    return get<std::vector<double>>("Jacob");
  }

  void set_Jacob(const int k, double val)
  {
    auto &jacob = get_ref<std::vector<double>>("Jacob");
    if (k >= 0 && static_cast<size_t>(k) < jacob.size())
    {
      jacob[k] = val;
    }
    else
    {
      throw std::out_of_range("set_Jacob: index " + std::to_string(k) + " out of range (size = " + std::to_string(jacob.size()) + ")");
    }
  }

  void resize_Jacob(int n)
  {
    if (n < 0)
    {
      throw std::invalid_argument("resize_Jacob: size cannot be negative (" + std::to_string(n) + ")");
    }
    get_ref<std::vector<double>>("Jacob").resize(static_cast<size_t>(n));
  }

  void append_Jacob(double value)
  {
    get_ref<std::vector<double>>("Jacob").push_back(value);
  }

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

#if 0
void print() const {
  std::cout << "Label: " << get_label() << "\n";
  std::cout << "Index_1: ";
  for (int idx : get_index_1()) std::cout << idx << " ";
  std::cout << "\nCst Scl: " << get_cst_scl()
	    << "\nCst: " << get_cst()
	    << "\nJacob: ";
  for (double val : get_Jacob()) std::cout << val << " ";
  std::cout << std::endl;
}

#else

void Lie_gen_class::print_label(void) const
{
  printf(" %5s", get_label().c_str());
}

void Lie_gen_class::print_index_1(void) const
{
  printf(" [");
  for (int k = 0; k < (int)get_index_1().size(); k++)
    printf(" %d", get_index_1()[k]);
  printf("]");
}

void Lie_gen_class::print_index_2(void) const
{
  printf(" [");
  for (int k = 0; k < (int)get_index_2().size(); k++)
    printf(" %d", get_index_2()[k]);
  printf("]");
}

void Lie_gen_class::print_cst_scl(void) const
{
  printf(" %6.1e", get_cst_scl());
}

void Lie_gen_class::print_cst(void) const
{
  printf(" %10.3e", get_cst());
}

void Lie_gen_class::print_Jacob(void) const
{
  for (int k = 0; k < (int)get_Jacob().size(); k++)
    printf(" %10.3e", get_Jacob()[k]);
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
#endif

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

  Lie_term.set("label", label);
  Lie_term.set("index_1", ind);
  Lie_term.set("cst_scl", scl);
  // Don't scale yet.
  Lie_term.set("cst", h_ijklm(g, ind));

  return Lie_term;
}

void prt_Lie_gen(const std::string &str, int &i0, const int n, const std::vector<Lie_gen_class> &Lie_gen)
{
  printf("\n%s\n", str.c_str());
  for (int k = i0; k < i0 + n; k++)
    Lie_gen[k].print();
  i0 += n;
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

std::vector<Lie_gen_class> get_rdts(const tps &g, const Config &cfg)
{
  std::vector<Lie_gen_class> rdts;

  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[0][0], 1, 0, 0, 0, 2));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[0][1], 2, 0, 0, 0, 1));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[0][2], 0, 0, 2, 0, 1));

  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[1][0], 1, 0, 1, 1, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[1][1], 2, 1, 0, 0, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[1][2], 3, 0, 0, 0, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[1][3], 1, 0, 0, 2, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[1][4], 1, 0, 2, 0, 0));

  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][0], 4, 0, 0, 0, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][1], 3, 1, 0, 0, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][2], 2, 0, 2, 0, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][3], 1, 1, 2, 0, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][4], 2, 0, 1, 1, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][5], 0, 0, 3, 1, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][6], 2, 0, 0, 2, 0));
  rdts.push_back(get_Lie_gen("g", g, cfg.scl_h[2][7], 0, 0, 4, 0, 0));

  return rdts;
}

void prt_rdts(const param_type &bns, const std::vector<Lie_gen_class> &rdts)
{
  int k = 0;
  prt_Lie_gen("3rd Order Chromatic terms:", k, 3, rdts);
  prt_Lie_gen("3rd Order Geometric terms:", k, 5, rdts);
  prt_Lie_gen("4th Order Geometric terms:", k, 8, rdts);
}

std::vector<Lie_gen_class> get_adts(const tps &K, const Config &cfg)
{
  std::vector<Lie_gen_class> adts;

  adts.push_back(get_Lie_gen("K", K, cfg.scl_a[0][0], 2, 2, 0, 0, 0));
  adts.push_back(get_Lie_gen("K", K, cfg.scl_a[0][1], 1, 1, 1, 1, 0));
  adts.push_back(get_Lie_gen("K", K, cfg.scl_a[0][2], 0, 0, 2, 2, 0));

  adts.push_back(get_Lie_gen("K", K, cfg.scl_a[1][0], 2, 2, 0, 0, 1));
  adts.push_back(get_Lie_gen("K", K, cfg.scl_a[1][1], 1, 1, 1, 1, 1));
  adts.push_back(get_Lie_gen("K", K, cfg.scl_a[1][2], 0, 0, 2, 2, 1));

  if (NO >= 6)
  {
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[2][0], 3, 3, 0, 0, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[2][1], 2, 2, 1, 1, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[2][2], 1, 1, 2, 2, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[2][3], 0, 0, 3, 3, 0));
  }

  if (NO >= 7)
  {
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[3][0], 3, 3, 0, 0, 1));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[3][1], 2, 2, 1, 1, 1));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[3][2], 1, 1, 2, 2, 1));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[3][3], 0, 0, 3, 3, 1));
  }
  // ?: CHECK IF THIS SHOULD BE NO>=8 OR NO>=9 ... Also below 10 instead of 9?
  if (NO >= 9)
  {
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[4][0], 4, 4, 0, 0, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[4][1], 3, 3, 1, 1, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[4][2], 2, 2, 2, 2, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[4][3], 1, 1, 3, 3, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[4][4], 0, 0, 4, 4, 0));
  }

  if (NO >= 11) // ? CHECK IF THIS SHOULD BE NO>=10
  {
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[5][0], 5, 5, 0, 0, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[5][1], 4, 4, 1, 1, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[5][2], 3, 3, 2, 2, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[5][3], 2, 2, 3, 3, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[5][4], 1, 1, 4, 4, 0));
    adts.push_back(get_Lie_gen("K", K, cfg.scl_a[5][5], 0, 0, 5, 5, 0));
  }

  return adts;
}

void prt_adts(const param_type &bns, const std::vector<Lie_gen_class> &adts)
{
  int k = 0;
  prt_Lie_gen("4th Order Anharmonic ADTs:", k, 3, adts);
  prt_Lie_gen("4th Order Anharmonic Cross Terms:", k, 3, adts);
  if (NO >= 6)
    prt_Lie_gen("6th Order Anharmonic ADTs:", k, 4, adts);
  if (NO >= 7)
    prt_Lie_gen("6th Order Anharmonic Cross Terms:", k, 4, adts);
  if (NO >= 9)
    prt_Lie_gen("8th Order Anharmonic ADTs:", k, 5, adts);
  if (NO >= 11)
    prt_Lie_gen("10th Order Anharmonic ADTs:", k, 6, adts);
}

std::vector<Lie_gen_class> get_xi(const tps &K, const Config &cfg)
{
  std::vector<Lie_gen_class> xi;

  xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[1][0], 1, 1, 0, 0, 1));
  xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[1][1], 0, 0, 1, 1, 1));

  xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[2][0], 1, 1, 0, 0, 2));
  xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[2][1], 0, 0, 1, 1, 2));

  if (NO >= 6)
  {
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[3][0], 1, 1, 0, 0, 3));
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[3][1], 0, 0, 1, 1, 3));
  }

  if (NO >= 7)
  {
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[4][0], 1, 1, 0, 0, 4));
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[4][1], 0, 0, 1, 1, 4));
  }

  if (NO >= 8)
  {
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[5][0], 1, 1, 0, 0, 5));
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[5][1], 0, 0, 1, 1, 5));
  }

  if (NO >= 9)
  {
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[6][0], 1, 1, 0, 0, 6));
    xi.push_back(get_Lie_gen("K", K, cfg.scl_ksi[6][1], 0, 0, 1, 1, 6));
  }

  return xi;
}

void prt_xi(const param_type &bns, const std::vector<Lie_gen_class> &xi)
{
  int k = 0;
  prt_Lie_gen("Linear Chromaticity:", k, 2, xi);
  prt_Lie_gen("2nd Order Chromaticity:", k, 2, xi);
  if (NO >= 6)
    prt_Lie_gen("3rd Order Chromaticity:", k, 2, xi);
  if (NO >= 7)
    prt_Lie_gen("4th Order Chromaticity:", k, 2, xi);
  if (NO >= 8)
    prt_Lie_gen("5th Order Chromaticity:", k, 2, xi);
  if (NO >= 9)
    prt_Lie_gen("6th Order Chromaticity:", k, 2, xi);
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

  Lie_term.set("label", label);
  Lie_term.set("index_1", ind_1);
  Lie_term.set("index_2", ind_2);
  Lie_term.set("cst_scl", scl);
  // Don't scale yet.
  Lie_term.set("cst", K_avg_2);

  return Lie_term;
}

void prt_Lie_K_avg_gen(const std::string &str, int &i0, const int n, const std::vector<Lie_gen_class> &Lie_gen)
{
  printf("\n%s\n", str.c_str());
  for (int k = i0; k < i0 + n; k++)
  {
    Lie_gen[k].print_label();
    Lie_gen[k].print_index_1();
    Lie_gen[k].print_index_2();
    Lie_gen[k].print_cst_scl();
    Lie_gen[k].print_cst();
    Lie_gen[k].print_Jacob();
    printf("\n");
  }
  i0 += n;
}

std::vector<Lie_gen_class> get_K_avg(const tps &K, const Config &cfg)
{
  tps avg;
  std::vector<Lie_gen_class> K_avg_2;

  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[0], 2, 2, 0, 0, 0, 4, 4, 0, 0, 0));

  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[1], 0, 0, 2, 2, 0, 0, 0, 4, 4, 0));

  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[2], 1, 1, 1, 1, 0, 2, 2, 1, 1, 0));
  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[3], 1, 1, 1, 1, 0, 1, 1, 2, 2, 0));

  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[4], 1, 1, 0, 0, 2, 1, 1, 0, 0, 4));
  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[5], 1, 1, 0, 0, 3, 1, 1, 0, 0, 5));

  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[6], 0, 0, 1, 1, 2, 0, 0, 1, 1, 4));
  K_avg_2.push_back(get_Lie_K_avg_gen("<K>", K, cfg.scl_K_avg[7], 0, 0, 1, 1, 3, 0, 0, 1, 1, 5));

  return K_avg_2;
}

void prt_K_avg(const param_type &bns, const std::vector<Lie_gen_class> &K_avg)
{
  int k = 0;
  prt_Lie_K_avg_gen("<K> Terms:", k, 8, K_avg);
}

void compute_Jacob(const param_type &bns, const int prm_no, const ss_vect<tps> &Id_scl, const Config &cfg, std::vector<Lie_gen_class> &Lie_gen)
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
    const auto &label = term.get_label();
    const auto &index_1 = term.get<std::vector<int>>("index_1");
    const auto &index_2 = term.get<std::vector<int>>("index_2");

    if (label == "g")
      term.append_Jacob(bn_scl * h_ijklm_p(g_im, index_1));
    else if (label == "K")
      term.append_Jacob(bn_scl * h_ijklm_p(K_re, index_1));
    else if (label == "<K>")
      term.append_Jacob(bn_scl * (h_ijklm_p(K_re, index_1) + h_ijklm_p(K_re, index_2)));
    else
    {
      fprintf(stderr, "\ncompute_Jacob - undef. case: %s\n", label.c_str());
      exit(1);
    }
  }
}

inline void compute_Jacob(const ss_vect<tps> &Id_scl, const param_type &bns, const Config &cfg, std::vector<Lie_gen_class> &Lie_gen)
{
  for (int k = 0; k < bns.n_prm; k++)
    compute_Jacob(bns, k, Id_scl, cfg, Lie_gen);
}

void Lie_gen_scl(std::vector<Lie_gen_class> &Lie_gen)
{
  for (auto j = 0; j < (int)Lie_gen.size(); j++)
  {
    Lie_gen[j].set("cst", Lie_gen[j].get_cst_scl() * Lie_gen[j].get_cst());
    for (int k = 0; k < (int)Lie_gen[j].get_Jacob().size(); k++)
      Lie_gen[j].set_Jacob(k, Lie_gen[j].get_cst_scl() * Lie_gen[j].get_Jacob()[k]);
  }
}

void get_Lie_gen(std::vector<Lie_gen_class> &Lie_gens, const std::vector<Lie_gen_class> &Lie_gen)
{
  for (const auto &lg : Lie_gen)
    Lie_gens.push_back(lg);
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

  auto rdts = get_rdts(g_im, cfg);
  compute_Jacob(Id_scl, bns, cfg, rdts);
  // Now scale.
  Lie_gen_scl(rdts);

  auto adts = get_adts(K_re, cfg);
  compute_Jacob(Id_scl, bns, cfg, adts);
  Lie_gen_scl(adts);

  auto xi_nl = get_xi(K_re, cfg);
  compute_Jacob(Id_scl, bns, cfg, xi_nl);
  Lie_gen_scl(xi_nl);

  auto K_avg = get_K_avg(K_re, cfg);
  compute_Jacob(Id_scl, bns, cfg, K_avg);
  Lie_gen_scl(K_avg);

  prt_system(bns);
  prt_rdts(bns, rdts);
  prt_adts(bns, adts);
  prt_xi(bns, xi_nl);
  prt_K_avg(bns, K_avg);

  get_Lie_gen(Lie_gen, rdts);
  get_Lie_gen(Lie_gen, adts);
  get_Lie_gen(Lie_gen, xi_nl);
  get_Lie_gen(Lie_gen, K_avg);

  return Lie_gen;
}

void get_system(const int m, const int n, const std::vector<Lie_gen_class> &Lie_gen, double **A, double *b)
{
  const bool debug = false;

  auto j1 = 1;
  for (const auto &term : Lie_gen)
  {
    b[j1] = -term.get_cst();

    int k1 = 1;
    for (const auto &jac : term.get_Jacob())
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

void set_Fam(param_type &bns, const int k, const double scl, const double *dbnL, const Config &cfg)
{
  double bnL_ext;

  bnL_ext = get_bnL(bns.Fnum[k], 1, bns.n[k]) + scl * bns.bnL_scl[k] * dbnL[k + 1];
  bns.bnL[k] = bnL_internal(bnL_ext, bns.bnL_min[k], bns.bnL_max[k]);
  set_bnL(bns.Fnum[k], bns.n[k], bnL_bounded(bns.bnL[k], bns.bnL_min[k], bns.bnL_max[k]));
}

void set_bnL(const double scl, const double *dbnL, param_type &bns, const Config &cfg)
{
  for (std::size_t k = 0; k < static_cast<std::size_t>(bns.n_prm); ++k)
    set_Fam(bns, k, scl, dbnL, cfg);
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

  double **A, **U, **V, *w, *b, *dbnL, *bnL_max, *bnL;

  printf("\nsvd:\n  m = %d n = %d\n", m, n);

  A = dmatrix(1, m, 1, n);
  U = dmatrix(1, m, 1, n);
  V = dmatrix(1, n, 1, n);
  w = dvector(1, n);
  b = dvector(1, m);
  dbnL = dvector(1, n);
  bnL_max = dvector(1, n);
  bnL = dvector(1, n);

  get_system(m, n, Lie_gen, A, b);

#if 1 // ? What is the difference between this and SVD_lim ?
  dmcopy(A, m, n, U);
  dsvdcmp(U, m, n, w, V);
  get_sing_val(n, w, svd_n_cut, cfg);

  dsvbksb(U, w, V, m, n, b, dbnL);
#else
  for (int k = 0; k < n; k++)
  {
    auto Fnum = bns.Fnum[k];
    bnL_max[k + 1] = (bns.L[k] == 0e0) ? bns.bnL_max[k] / bns.bnL_scl[k] : bns.bnL_max[k] * bns.L[k] / bns.bnL_scl[k];
    bnL[k + 1] = get_bnL(Fnum, 1, bns.n[k]) / bns.bnL_scl[k];
  }

  SVD_lim(m, n, A, b, bnL_max, 1e-11, bnL, dbnL);

#endif

  set_bnL(scl, dbnL, bns, cfg);
  bns.print();
  prt_bn(bns);
  prt_mfile("flat_file.fit");

  free_dmatrix(A, 1, m, 1, n);
  free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n);
  free_dvector(w, 1, n);
  free_dvector(b, 1, m);
  free_dvector(dbnL, 1, n);
  free_dvector(bnL_max, 1, n);
  free_dvector(bnL, 1, n);
}

void no_mpoles(const int n)
{
  printf("\nzeroing multipoles: %d\n", n);
  for (auto j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}

void get_bns(param_type &bns, const Config &cfg)
{
  const int lat = cfg.lat;

  if (cfg.b_3_zero)
    no_mpoles(Sext);
  if (cfg.b_4_zero)
    no_mpoles(Oct);

  switch (lat)
  {
  case 1:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 2:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1_f1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2_f1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3_f1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4_f1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1_f1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2_f1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3_f1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 3:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1_h1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2_h1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3_h1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4_h1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1_h1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2_h1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3_h1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 4:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1_h2", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2_h2", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3_h2", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4_h2", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      if (false)
        bns.add_Fam("s5_h2", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1_h2", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2_h2", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3_h2", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 51:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 52:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3a_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3b_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3c_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4a_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4b_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 53:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1_n1_sl", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2_n1_sl", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3_n1_sl", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 54:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("s1_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s2_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3a_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3b_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s3c_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("s4_n1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("o1_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o2_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("o3_n1", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  case 6:
    if (cfg.b_3_opt)
    {
      bns.add_Fam("sfm", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sfi", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sdqd_1", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sdqd_2", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sdqd_3", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sdqd_4", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sdqd_5", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sdendq", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
      bns.add_Fam("sfo", Sext, cfg.bnL_min[Sext], cfg.bnL_max[Sext], cfg.bnL_scl[Sext]);
    }

    if (cfg.b_4_opt)
    {
      bns.add_Fam("oxxo", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("oxyo", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
      bns.add_Fam("oyyo", Oct, cfg.bnL_min[Oct], cfg.bnL_max[Oct], cfg.bnL_scl[Oct]);
    }
    break;
  default:
    fprintf(stderr, "\nget_bns - unsupported lat value: %d\n", lat);
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

// Config parsing and serialization helpers were extracted to
// gustavo_options_config.h to keep this file focused on optimization code.

int main(int argc, char *argv[])
{
  double twoJ;
  ss_vect<tps> Id_scl;
  param_type bns;
  std::vector<Lie_gen_class> Lie_gen;

  set_state();

  // Parse command-line arguments
  Config cfg = parse_arguments(argc, argv);
  write_config_to_file(cfg, "final_config.cfg");

  rd_mfile(cfg.lattice_file.c_str(), elem);
  rd_mfile(cfg.lattice_file.c_str(), elem_tps);

  if (true)
  {
    printf("\nLattice file: %s\n", cfg.lattice_file.c_str());
    printf("Lattice family: %d\n", cfg.lat);
    printf("Phase: %d\n", cfg.phase);
    printf("Lat family: %d\n", cfg.lat);
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
  write_config_to_file(cfg, "config_output.txt");
  return 0;
}
