#ifndef GUSTAVO_OPTIONS_CONFIG_H
#define GUSTAVO_OPTIONS_CONFIG_H

#include <algorithm>
#include <array>
#include <boost/program_options.hpp>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace po = boost::program_options;

struct Config
{
  int phase = 2;        // weight preset, see apply_phase_preset()
  // Families to fit, by their names in the lattice.
  std::vector<std::string> sext_families, oct_families;
  int svd_n_cut = 0;    // smallest singular values zeroed
  int max_iter = 300;

  bool b_3_opt = true;  // fit the sextupole families
  bool b_4_opt = true;  // fit the octupole families
  bool b_3_zero = false;  // zero all sextupoles before the first iteration
  bool b_4_zero = false;  // zero all octupoles before the first iteration

  // Acceptance the Lie terms are scaled to: amplitudes [m], beta at injection [m].
  std::array<double, 2> A_max = {{6e-3, 3e-3}};
  double delta_max = 6e-2;
  std::array<double, 2> beta_inj = {{3.7, 3.9}};

  // Strength scale and limits, indexed by multipole order.
  std::array<double, 6> bnL_scl = {{0e0, 0e0, 0e0, 1e0, 5e1 / 1e2, 5 * 1e4}};
  std::array<double, 6> bnL_min = {{0e0, 0e0, 0e0, -5e2, -5.0e4, -1.5e5}};
  std::array<double, 6> bnL_max = {{0e0, 0e0, 0e0, 5e2, 5.0e4, 1.5e5}};
  // Scales the sorted singular values, largest first; one value per fitted family.
  // Empty until resolve_scl_svd() fills it from svd_scl_sext/svd_scl_oct, unless
  // given explicitly.
  std::vector<double> scl_svd;
  double svd_scl_sext = 1e0;
  double svd_scl_oct = 5e2;

  // Lie term weights, one row per term group, one value per term, in the order of
  // gustavo_options_terms.h.
  std::array<std::vector<double>, 3> scl_h = {{
    {1e-2, 1e-2, 1e-2},                                // h10002 h20001 h00201
    {1e-2, 1e-2, 1e-2, 1e-2, 1e-2},                    // h10110 h21000 h30000 h10020 h10200
    {1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2},  // h40000 h31000 h20200 h11200 h20110 h00310 h20020 h00400
  }};
  std::array<std::vector<double>, 6> scl_a = {{
    {5e0, 5e0, 5e0},                 // K22000 K11110 K00220
    {5e0, 5e0, 5e0},                 // K22001 K11111 K00221
    {5e0, 5e0, 5e0, 5e0},            // K33000 K22110 K11220 K00330
    {5e0, 5e0, 5e0, 5e0},            // K33001 K22111 K11221 K00331
    {5e0, 5e0, 5e0, 5e0, 5e0},       // K44000 K33110 K22220 K11330 K00440 (NO >= 9)
    {5e0, 5e0, 5e0, 5e0, 5e0, 5e0},  // K55000 K44110 K33220 K22330 K11440 K00550 (NO >= 11)
  }};
  // scl_ksi[r]: K1100n K0011n, the order n = r + 1 chromaticity (key scl_ksi_<n>).
  std::array<std::array<double, 2>, 6> scl_ksi = {{
    {{1e2, 1e2}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
  }};
  // <K22000,K44000> <K00220,K00440> (both NO >= 9) <K11110,K22110> <K11110,K11220>
  // <K11002,K11004> <K11003,K11005> <K00112,K00114> <K00113,K00115>
  std::array<double, 8> scl_K_avg = {{1e-3, 1e-3, 1e-3, 1e-3, 1e2, 1e2, 1e2, 1e2}};

  double step = 4 * 0.15;  // fraction of the SVD correction applied per iteration
  int flat_file_every = 0;  // write flat_file_<k>.fit every N iterations; 0: never

  std::string lattice_file;
};

// Weight presets for gustavo.cc's two phases; the defaults above are phase 2.
// Precedence: defaults < this preset < rows given in the config file or on the
// command line. Unlike gustavo.cc, phase 1 does not zero the multipoles; set
// b3_zero and b4_zero for that.
inline void apply_phase_preset(Config &cfg)
{
  if (cfg.phase != 1)
    return;
  cfg.scl_ksi = {{
    {{1e2, 1e2}},
    {{5e0, 5e0}},
    {{5e0, 5e0}},
    {{5e0, 5e0}},
    {{5e0, 5e0}},
    {{5e0, 5e0}},
  }};
  cfg.scl_a = {{
    {1e0, 1e0, 1e0},
    {1e0, 1e0, 1e0},
    {1e0, 1e0, 1e0, 1e0},
    {1e0, 1e0, 1e0, 1e0},
    {1e0, 1e0, 1e0, 1e0, 1e0},
    {1e0, 1e0, 1e0, 1e0, 1e0, 1e0},
  }};
  cfg.scl_K_avg = {{1e-3, 1e-3, 1e-3, 1e-3, 5e3, 5e3, 5e3, 5e3}};
}

// A weight row: a comma- or space-separated list of doubles in the config.
// Rows with n == 0 (scl_svd) take any number of values; the others exactly n.
struct WeightRow
{
  std::string key;
  std::string help;
  std::size_t n;
  double *values;               // fixed-length rows
  std::vector<double> *vector;  // scl_svd
};

// Every weight row, in config-file order.
inline std::vector<WeightRow> weight_rows(Config &cfg)
{
  std::vector<WeightRow> rows;
  const auto add = [&rows](const std::string &key, const std::string &help, double *values, const std::size_t n) { rows.push_back({key, help, n, values, nullptr}); };

  add("bnl_scl", "multipole strength scale per order", cfg.bnL_scl.data(), cfg.bnL_scl.size());
  rows.push_back({"scl_svd",
                  "singular value scales, largest first (at least one per fitted family); default: svd_scl_sext per sextupole family, then svd_scl_oct per octupole "
                  "family",
                  0,
                  nullptr,
                  &cfg.scl_svd});
  const char *h_rows[] = {"3rd-order chromatic", "3rd-order geometric", "4th-order geometric"};
  for (std::size_t r = 0; r < cfg.scl_h.size(); ++r)
    add("scl_h_" + std::to_string(r), std::string(h_rows[r]) + " term weights", cfg.scl_h[r].data(), cfg.scl_h[r].size());
  for (std::size_t r = 0; r < cfg.scl_a.size(); ++r)
    add("scl_a_" + std::to_string(r), "amplitude-dependent tune shift weights, row " + std::to_string(r), cfg.scl_a[r].data(), cfg.scl_a[r].size());
  for (std::size_t r = 0; r < cfg.scl_ksi.size(); ++r)
    add("scl_ksi_" + std::to_string(r + 1), "order-" + std::to_string(r + 1) + " chromaticity weights", cfg.scl_ksi[r].data(), cfg.scl_ksi[r].size());
  add("scl_K_avg", "tune footprint average weights", cfg.scl_K_avg.data(), cfg.scl_K_avg.size());

  return rows;
}

inline po::options_description get_options_description(void)
{
  const Config d;
  po::options_description desc("Allowed options (weight rows are comma- or space-separated numbers)");
  desc.add_options()
    ("help,h", "produce help message")
    ("config", po::value<std::string>(), "configuration file")
    ("lattice", po::value<std::string>()->required(), "input lattice file path (positional)")
    ("phase", po::value<int>()->default_value(d.phase), "optimization phase: 1 or 2")
    ("sext_families", po::value<std::string>(), "sextupole families to fit (comma- or space-separated names)")
    ("oct_families", po::value<std::string>(), "octupole families to fit (comma- or space-separated names)")
    ("svd_n_cut", po::value<int>()->default_value(d.svd_n_cut), "number of singular values to remove during SVD solving")
    ("max_iter", po::value<int>()->default_value(d.max_iter), "maximum optimization iterations")
    ("step", po::value<double>()->default_value(d.step), "optimization step size; must satisfy 0 < step < 1")
    ("flat_file_every", po::value<int>()->default_value(d.flat_file_every), "write flat_file_k.fit every N iterations; 0 disables")
    ("b3_opt", po::value<bool>()->default_value(d.b_3_opt), "optimize sextupoles")
    ("b4_opt", po::value<bool>()->default_value(d.b_4_opt), "optimize octupoles")
    ("b3_zero", po::value<bool>()->default_value(d.b_3_zero), "zero sextupoles initially")
    ("b4_zero", po::value<bool>()->default_value(d.b_4_zero), "zero octupoles initially")
    ("svd_scl_sext", po::value<double>()->default_value(d.svd_scl_sext), "scl_svd value for each sextupole family when scl_svd is not given")
    ("svd_scl_oct", po::value<double>()->default_value(d.svd_scl_oct), "scl_svd value for each octupole family when scl_svd is not given");

  Config rows_of = d;  // weight_rows() needs a non-const Config
  for (const auto &row : weight_rows(rows_of))
  {
    const std::string count = row.n ? std::to_string(row.n) + " values" : "list";
    desc.add_options()(row.key.c_str(), po::value<std::string>(), (row.help + " (" + count + ")").c_str());
  }

  return desc;
}

inline std::string trim_copy(const std::string &s)
{
  std::size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
    ++b;
  std::size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
    --e;
  return s.substr(b, e - b);
}

// Boost reports unknown keys too, but without the line number.
inline void validate_config_file_keys(const std::string &filename, const po::options_description &desc)
{
  std::ifstream in(filename.c_str());
  if (!in)
    throw std::runtime_error("could not open config file '" + filename + "'");

  std::string line;
  std::size_t line_no = 0;
  while (std::getline(in, line))
  {
    ++line_no;
    const std::size_t comment_pos = line.find('#');
    if (comment_pos != std::string::npos)
      line.erase(comment_pos);

    const std::string stripped = trim_copy(line);
    if (stripped.empty())
      continue;

    const std::size_t eq = stripped.find('=');
    if (eq == std::string::npos)
      throw std::runtime_error("invalid line " + std::to_string(line_no) + " in config file '" + filename + "': missing '='");

    const std::string key = trim_copy(stripped.substr(0, eq));
    if (key.empty())
      throw std::runtime_error("invalid line " + std::to_string(line_no) + " in config file '" + filename + "': empty key");

    if (key == "help" || !desc.find_nothrow(key, false))
      throw std::runtime_error("unknown key '" + key + "' at line " + std::to_string(line_no) + " in config file '" + filename + "'");
  }
}

inline std::vector<double> parse_csv(const std::string &opt_name, const std::string &s)
{
  if (s.empty())
    throw std::runtime_error("option '" + opt_name + "' cannot be empty");

  std::string tmp = s;
  for (auto &c : tmp)
    if (c == ',')
      c = ' ';

  std::istringstream iss(tmp);
  std::vector<double> vals;
  double v;
  while (iss >> v)
    vals.push_back(v);

  // Reject malformed values such as "1,2,foo".
  iss.clear();
  std::string trailing;
  if (iss >> trailing)
    throw std::runtime_error("option '" + opt_name + "' contains non-numeric token: '" + trailing + "'");

  return vals;
}

inline std::vector<std::string> parse_names(const std::string &s)
{
  std::string tmp = s;
  for (auto &c : tmp)
    if (c == ',')
      c = ' ';
  std::istringstream iss(tmp);
  std::vector<std::string> names;
  std::string name;
  while (iss >> name)
    names.push_back(name);
  return names;
}

inline std::string join_names(const std::vector<std::string> &names)
{
  std::string s;
  for (std::size_t i = 0; i < names.size(); ++i)
    s += (i ? "," : "") + names[i];
  return s;
}

inline Config parse_arguments(int argc, char *argv[])
{
  po::options_description desc = get_options_description();

  po::positional_options_description p;
  p.add("lattice", 1);

  po::variables_map vm;
  try
  {
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    if (vm.count("config"))
    {
      const std::string cfg_name = vm["config"].as<std::string>();
      validate_config_file_keys(cfg_name, desc);

      std::ifstream config_file(cfg_name);
      po::store(po::parse_config_file(config_file, desc), vm);
    }

    if (vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [options] LATTICE\n\n" << desc << "\n";
      exit(0);
    }

    po::notify(vm);
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    std::cerr << "Usage: " << argv[0] << " [options] LATTICE\n";
    std::cerr << "Try '" << argv[0] << " --help' for more information.\n";
    exit(1);
  }

  Config cfg;
  cfg.phase = vm["phase"].as<int>();
  if (vm.count("sext_families"))
    cfg.sext_families = parse_names(vm["sext_families"].as<std::string>());
  if (vm.count("oct_families"))
    cfg.oct_families = parse_names(vm["oct_families"].as<std::string>());

  if (cfg.phase != 1 && cfg.phase != 2)
  {
    std::cerr << "Error: phase must be 1 or 2 (got " << cfg.phase << ")\n";
    exit(1);
  }

  apply_phase_preset(cfg);

  cfg.svd_n_cut = vm["svd_n_cut"].as<int>();
  cfg.max_iter = vm["max_iter"].as<int>();
  cfg.step = vm["step"].as<double>();
  cfg.flat_file_every = vm["flat_file_every"].as<int>();
  cfg.b_3_opt = vm["b3_opt"].as<bool>();
  cfg.b_4_opt = vm["b4_opt"].as<bool>();
  cfg.b_3_zero = vm["b3_zero"].as<bool>();
  cfg.b_4_zero = vm["b4_zero"].as<bool>();
  cfg.lattice_file = vm["lattice"].as<std::string>();

  if (!(cfg.step > 0.0 && cfg.step < 1.0))
  {
    std::cerr << "Error: step must satisfy 0 < step < 1 (got " << cfg.step << ")\n";
    exit(1);
  }
  if (cfg.flat_file_every < 0)
  {
    std::cerr << "Error: flat_file_every must be >= 0 (got " << cfg.flat_file_every << ")\n";
    exit(1);
  }

  cfg.svd_scl_sext = vm["svd_scl_sext"].as<double>();
  cfg.svd_scl_oct = vm["svd_scl_oct"].as<double>();
  if (vm.count("scl_svd") && (!vm["svd_scl_sext"].defaulted() || !vm["svd_scl_oct"].defaulted()))
  {
    std::cerr << "Error: set either scl_svd or svd_scl_sext/svd_scl_oct, not both\n";
    exit(1);
  }

  try
  {
    for (const auto &row : weight_rows(cfg))
    {
      if (!vm.count(row.key))
        continue;
      const std::vector<double> vals = parse_csv(row.key, vm[row.key].as<std::string>());
      if (row.vector)
        *row.vector = vals;
      else if (vals.size() != row.n)
        throw std::runtime_error("option '" + row.key + "' expects exactly " + std::to_string(row.n) + " values, got " + std::to_string(vals.size()));
      else
        std::copy(vals.begin(), vals.end(), row.values);
    }
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    exit(1);
  }

  return cfg;
}

// Call once the fitted families are known.
inline void resolve_scl_svd(Config &cfg, const int n_sext, const int n_oct)
{
  const std::size_t n_prm = static_cast<std::size_t>(n_sext + n_oct);

  if (cfg.scl_svd.empty())
  {
    cfg.scl_svd.assign(n_sext, cfg.svd_scl_sext);
    cfg.scl_svd.insert(cfg.scl_svd.end(), n_oct, cfg.svd_scl_oct);
  }
  else if (cfg.scl_svd.size() < n_prm)
  {
    std::cerr << "Error: scl_svd has " << cfg.scl_svd.size() << " values but " << n_prm << " families are fitted\n";
    exit(1);
  }
}

// Shortest decimal form that reads back as the same double.
inline std::string format_double(const double x)
{
  char buf[32];
  for (int precision = 6; precision <= 17; ++precision)
  {
    std::snprintf(buf, sizeof(buf), "%.*g", precision, x);
    if (std::strtod(buf, nullptr) == x)
      break;
  }
  return buf;
}

template <typename Values> inline void write_csv(std::ostream &out, const Values &values)
{
  for (std::size_t i = 0; i < values.size(); ++i)
  {
    if (i)
      out << ",";
    out << format_double(values[i]);
  }
}

inline void write_config_to_file(const Config &cfg, const std::string &filename)
{
  std::ofstream out(filename);
  if (!out)
  {
    std::cerr << "Error: Could not write to config file: " << filename << "\n";
    return;
  }

  out << "phase=" << cfg.phase << "\n";
  out << "sext_families=" << join_names(cfg.sext_families) << "\n";
  out << "oct_families=" << join_names(cfg.oct_families) << "\n";
  out << "svd_n_cut=" << cfg.svd_n_cut << "\n";
  out << "max_iter=" << cfg.max_iter << "\n";
  out << "step=" << format_double(cfg.step) << "\n";
  out << "flat_file_every=" << cfg.flat_file_every << "\n";
  out << "b3_opt=" << (cfg.b_3_opt ? "true" : "false") << "\n";
  out << "b4_opt=" << (cfg.b_4_opt ? "true" : "false") << "\n";
  out << "b3_zero=" << (cfg.b_3_zero ? "true" : "false") << "\n";
  out << "b4_zero=" << (cfg.b_4_zero ? "true" : "false") << "\n";
  out << "lattice=" << cfg.lattice_file << "\n";

  Config copy = cfg; // weight_rows() hands out writable pointers
  for (const auto &row : weight_rows(copy))
  {
    out << row.key << "=";
    if (row.vector)
      write_csv(out, *row.vector);
    else
      write_csv(out, std::vector<double>(row.values, row.values + row.n));
    out << "\n";
  }
}

#endif
