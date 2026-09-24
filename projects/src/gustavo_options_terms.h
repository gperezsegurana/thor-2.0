#ifndef GUSTAVO_OPTIONS_TERMS_H
#define GUSTAVO_OPTIONS_TERMS_H

// The Lie terms gustavo_options drives to zero, one table row per term.
//
// A term is included when the map carries it: its order (sum of the indices;
// for a <K> average, the higher of its two index sets) plus one for the
// derivative with respect to the magnet strength must not exceed NO. Rows
// print in table order, under their section heading.

#include <array>
#include <vector>

#include "gustavo_options_config.h"

// Where a term's value comes from.
enum class TermSource
{
  g,     // imaginary part of the generator g
  K,     // real part of the normal form K
  K_avg  // real part of K, averaged over two index sets
};

// Which Config weight row scales the term.
enum class TermWeights
{
  scl_h,
  scl_a,
  scl_ksi,
  scl_K_avg
};

struct TermDef
{
  const char *section;
  TermSource source;
  TermWeights weights;
  int row, col;  // row is ignored for scl_K_avg
  std::array<int, 5> ind_1;
  std::array<int, 5> ind_2;  // K_avg only
};

inline const std::vector<TermDef> &term_table(void)
{
  using S = TermSource;
  using W = TermWeights;

  static const std::vector<TermDef> table = {
    {"3rd Order Chromatic terms:", S::g, W::scl_h, 0, 0, {{1, 0, 0, 0, 2}}, {}},
    {"3rd Order Chromatic terms:", S::g, W::scl_h, 0, 1, {{2, 0, 0, 0, 1}}, {}},
    {"3rd Order Chromatic terms:", S::g, W::scl_h, 0, 2, {{0, 0, 2, 0, 1}}, {}},

    {"3rd Order Geometric terms:", S::g, W::scl_h, 1, 0, {{1, 0, 1, 1, 0}}, {}},
    {"3rd Order Geometric terms:", S::g, W::scl_h, 1, 1, {{2, 1, 0, 0, 0}}, {}},
    {"3rd Order Geometric terms:", S::g, W::scl_h, 1, 2, {{3, 0, 0, 0, 0}}, {}},
    {"3rd Order Geometric terms:", S::g, W::scl_h, 1, 3, {{1, 0, 0, 2, 0}}, {}},
    {"3rd Order Geometric terms:", S::g, W::scl_h, 1, 4, {{1, 0, 2, 0, 0}}, {}},

    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 0, {{4, 0, 0, 0, 0}}, {}},
    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 1, {{3, 1, 0, 0, 0}}, {}},
    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 2, {{2, 0, 2, 0, 0}}, {}},
    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 3, {{1, 1, 2, 0, 0}}, {}},
    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 4, {{2, 0, 1, 1, 0}}, {}},
    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 5, {{0, 0, 3, 1, 0}}, {}},
    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 6, {{2, 0, 0, 2, 0}}, {}},
    {"4th Order Geometric terms:", S::g, W::scl_h, 2, 7, {{0, 0, 4, 0, 0}}, {}},

    {"4th Order Anharmonic ADTs:", S::K, W::scl_a, 0, 0, {{2, 2, 0, 0, 0}}, {}},
    {"4th Order Anharmonic ADTs:", S::K, W::scl_a, 0, 1, {{1, 1, 1, 1, 0}}, {}},
    {"4th Order Anharmonic ADTs:", S::K, W::scl_a, 0, 2, {{0, 0, 2, 2, 0}}, {}},

    {"4th Order Anharmonic Cross Terms:", S::K, W::scl_a, 1, 0, {{2, 2, 0, 0, 1}}, {}},
    {"4th Order Anharmonic Cross Terms:", S::K, W::scl_a, 1, 1, {{1, 1, 1, 1, 1}}, {}},
    {"4th Order Anharmonic Cross Terms:", S::K, W::scl_a, 1, 2, {{0, 0, 2, 2, 1}}, {}},

    {"6th Order Anharmonic ADTs:", S::K, W::scl_a, 2, 0, {{3, 3, 0, 0, 0}}, {}},
    {"6th Order Anharmonic ADTs:", S::K, W::scl_a, 2, 1, {{2, 2, 1, 1, 0}}, {}},
    {"6th Order Anharmonic ADTs:", S::K, W::scl_a, 2, 2, {{1, 1, 2, 2, 0}}, {}},
    {"6th Order Anharmonic ADTs:", S::K, W::scl_a, 2, 3, {{0, 0, 3, 3, 0}}, {}},

    {"6th Order Anharmonic Cross Terms:", S::K, W::scl_a, 3, 0, {{3, 3, 0, 0, 1}}, {}},
    {"6th Order Anharmonic Cross Terms:", S::K, W::scl_a, 3, 1, {{2, 2, 1, 1, 1}}, {}},
    {"6th Order Anharmonic Cross Terms:", S::K, W::scl_a, 3, 2, {{1, 1, 2, 2, 1}}, {}},
    {"6th Order Anharmonic Cross Terms:", S::K, W::scl_a, 3, 3, {{0, 0, 3, 3, 1}}, {}},

    {"8th Order Anharmonic ADTs:", S::K, W::scl_a, 4, 0, {{4, 4, 0, 0, 0}}, {}},
    {"8th Order Anharmonic ADTs:", S::K, W::scl_a, 4, 1, {{3, 3, 1, 1, 0}}, {}},
    {"8th Order Anharmonic ADTs:", S::K, W::scl_a, 4, 2, {{2, 2, 2, 2, 0}}, {}},
    {"8th Order Anharmonic ADTs:", S::K, W::scl_a, 4, 3, {{1, 1, 3, 3, 0}}, {}},
    {"8th Order Anharmonic ADTs:", S::K, W::scl_a, 4, 4, {{0, 0, 4, 4, 0}}, {}},

    {"10th Order Anharmonic ADTs:", S::K, W::scl_a, 5, 0, {{5, 5, 0, 0, 0}}, {}},
    {"10th Order Anharmonic ADTs:", S::K, W::scl_a, 5, 1, {{4, 4, 1, 1, 0}}, {}},
    {"10th Order Anharmonic ADTs:", S::K, W::scl_a, 5, 2, {{3, 3, 2, 2, 0}}, {}},
    {"10th Order Anharmonic ADTs:", S::K, W::scl_a, 5, 3, {{2, 2, 3, 3, 0}}, {}},
    {"10th Order Anharmonic ADTs:", S::K, W::scl_a, 5, 4, {{1, 1, 4, 4, 0}}, {}},
    {"10th Order Anharmonic ADTs:", S::K, W::scl_a, 5, 5, {{0, 0, 5, 5, 0}}, {}},

    {"Linear Chromaticity:", S::K, W::scl_ksi, 0, 0, {{1, 1, 0, 0, 1}}, {}},
    {"Linear Chromaticity:", S::K, W::scl_ksi, 0, 1, {{0, 0, 1, 1, 1}}, {}},
    {"2nd Order Chromaticity:", S::K, W::scl_ksi, 1, 0, {{1, 1, 0, 0, 2}}, {}},
    {"2nd Order Chromaticity:", S::K, W::scl_ksi, 1, 1, {{0, 0, 1, 1, 2}}, {}},
    {"3rd Order Chromaticity:", S::K, W::scl_ksi, 2, 0, {{1, 1, 0, 0, 3}}, {}},
    {"3rd Order Chromaticity:", S::K, W::scl_ksi, 2, 1, {{0, 0, 1, 1, 3}}, {}},
    {"4th Order Chromaticity:", S::K, W::scl_ksi, 3, 0, {{1, 1, 0, 0, 4}}, {}},
    {"4th Order Chromaticity:", S::K, W::scl_ksi, 3, 1, {{0, 0, 1, 1, 4}}, {}},
    {"5th Order Chromaticity:", S::K, W::scl_ksi, 4, 0, {{1, 1, 0, 0, 5}}, {}},
    {"5th Order Chromaticity:", S::K, W::scl_ksi, 4, 1, {{0, 0, 1, 1, 5}}, {}},
    {"6th Order Chromaticity:", S::K, W::scl_ksi, 5, 0, {{1, 1, 0, 0, 6}}, {}},
    {"6th Order Chromaticity:", S::K, W::scl_ksi, 5, 1, {{0, 0, 1, 1, 6}}, {}},

    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 0, {{2, 2, 0, 0, 0}}, {{4, 4, 0, 0, 0}}},
    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 1, {{0, 0, 2, 2, 0}}, {{0, 0, 4, 4, 0}}},
    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 2, {{1, 1, 1, 1, 0}}, {{2, 2, 1, 1, 0}}},
    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 3, {{1, 1, 1, 1, 0}}, {{1, 1, 2, 2, 0}}},
    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 4, {{1, 1, 0, 0, 2}}, {{1, 1, 0, 0, 4}}},
    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 5, {{1, 1, 0, 0, 3}}, {{1, 1, 0, 0, 5}}},
    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 6, {{0, 0, 1, 1, 2}}, {{0, 0, 1, 1, 4}}},
    {"<K> Terms:", S::K_avg, W::scl_K_avg, 0, 7, {{0, 0, 1, 1, 3}}, {{0, 0, 1, 1, 5}}},
  };

  return table;
}

inline int term_order(const std::array<int, 5> &ind)
{
  return ind[0] + ind[1] + ind[2] + ind[3] + ind[4];
}

inline bool term_active(const TermDef &t, const int no)
{
  int order = term_order(t.ind_1);
  if (t.source == TermSource::K_avg && term_order(t.ind_2) > order)
    order = term_order(t.ind_2);
  return order + 1 <= no;
}

inline double term_weight(const TermDef &t, const Config &cfg)
{
  switch (t.weights)
  {
  case TermWeights::scl_h:
    return cfg.scl_h[t.row][t.col];
  case TermWeights::scl_a:
    return cfg.scl_a[t.row][t.col];
  case TermWeights::scl_ksi:
    return cfg.scl_ksi[t.row][t.col];
  case TermWeights::scl_K_avg:
    return cfg.scl_K_avg[t.col];
  }
  return 0e0;
}

#endif
