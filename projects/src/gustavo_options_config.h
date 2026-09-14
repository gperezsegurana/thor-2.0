#ifndef GUSTAVO_OPTIONS_CONFIG_H
#define GUSTAVO_OPTIONS_CONFIG_H

#include <array>
#include <boost/program_options.hpp>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace po = boost::program_options;

struct Config
{
  int phase = 2;
  int lat = 51;
  int svd_n_cut = 0;
  int max_iter = 300;

  bool b_3_opt = true;
  bool b_4_opt = true;
  bool b_3_zero = false;
  bool b_4_zero = false;

  std::array<double, 2> A_max = {{6e-3, 3e-3}};
  double delta_max = 6e-2;
  std::array<double, 2> beta_inj = {{3.7, 3.9}};

  std::array<double, 6> bnL_scl = {{0e0, 0e0, 0e0, 1e0, 5e1 / 1e2, 5 * 1e4}};
  std::array<double, 6> bnL_min = {{0e0, 0e0, 0e0, -5e2, -5.0e4, -1.5e5}};
  std::array<double, 6> bnL_max = {{0e0, 0e0, 0e0, 5e2, 5.0e4, 1.5e5}};
  std::array<double, 10> scl_svd = {{1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 5e2, 5e2, 5e2}};

  // Per-term Lie weights; rows map directly to get_rdts/get_adts/get_xi blocks.
  std::array<std::array<double, 8>, 3> scl_h = {{
    {{1e-2, 1e-2, 1e-2, 0e0, 0e0, 0e0, 0e0, 0e0}},
    {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0e0, 0e0, 0e0}},
    {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2}},
  }};
  std::array<std::array<double, 6>, 6> scl_a = {{
    {{5e0, 5e0, 5e0, 0e0, 0e0, 0e0}},
    {{5e0, 5e0, 5e0, 0e0, 0e0, 0e0}},
    {{5e0, 5e0, 5e0, 5e0, 0e0, 0e0}},
    {{5e0, 5e0, 5e0, 5e0, 0e0, 0e0}},
    {{5e0, 5e0, 5e0, 5e0, 5e0, 0e0}},
    {{5e0, 5e0, 5e0, 5e0, 5e0, 5e0}},
  }};
  std::array<std::array<double, 2>, 7> scl_ksi = {{
    {{0e0, 0e0}},
    {{1e2, 1e2}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
    {{1e0, 1e0}},
  }};
  std::array<double, 8> scl_K_avg = {{1e-3, 1e-3, 1e-3, 1e-3, 1e2, 1e2, 1e2, 1e2}};
  std::array<double, 2> scl_k_sum = {{0e2, 0e2}};

  double step = 4 * 0.15;
  int flat_file_every = 0;

  std::string lattice_file;
};

inline po::options_description get_options_description(void)
{
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")("phase", po::value<int>()->default_value(2), "optimization phase: 1 or 2")("lat",
                                                                                                                                  po::value<int>()->default_value(51),
                                                                                                                                  "lattice family selector used by get_bns (e.g. 1,2,3,4,51,52,53,6)")(
    "svd_n_cut",
    po::value<int>()->default_value(0),
    "number of singular values to remove during SVD solving")("max_iter", po::value<int>()->default_value(300), "maximum optimization iterations")(
    "step",
    po::value<double>()->default_value(4 * 0.15),
    "optimization step size; must satisfy 0 < step < 1")("flat_file_every", po::value<int>()->default_value(0), "write flat_file_k.fit every N iterations; 0 disables")(
    "b3_opt",
    po::value<bool>()->default_value(true),
    "optimize sextupoles")("b4_opt", po::value<bool>()->default_value(true), "optimize octupoles")("b3_zero", po::value<bool>()->default_value(false), "zero sextupoles initially")(
    "b4_zero",
    po::value<bool>()->default_value(false),
    "zero octupoles initially")("config", po::value<std::string>(), "configuration file")("lattice", po::value<std::string>()->required(), "input lattice file path (positional)")

    ("bnl_scl", po::value<std::string>(), "comma- or space-separated bnl_scl values (6 values)")("scl_svd", po::value<std::string>(), "comma- or space-separated scl_svd values (7 values)")(
      "scl_h_0",
      po::value<std::string>(),
      "comma- or space-separated scl_h row 0 values (8 values)")("scl_h_1", po::value<std::string>(), "comma- or space-separated scl_h row 1 values (8 values)")(
      "scl_h_2",
      po::value<std::string>(),
      "comma- or space-separated scl_h row 2 values (8 values)")("scl_a_0", po::value<std::string>(), "comma- or space-separated scl_a row 0 values (6 values)")(
      "scl_a_1",
      po::value<std::string>(),
      "comma- or space-separated scl_a row 1 values (6 values)")("scl_a_2", po::value<std::string>(), "comma- or space-separated scl_a row 2 values (6 values)")(
      "scl_a_3",
      po::value<std::string>(),
      "comma- or space-separated scl_a row 3 values (6 values)")("scl_a_4", po::value<std::string>(), "comma- or space-separated scl_a row 4 values (6 values)")(
      "scl_a_5",
      po::value<std::string>(),
      "comma- or space-separated scl_a row 5 values (6 values)")("scl_ksi_0", po::value<std::string>(), "comma- or space-separated scl_ksi row 0 values (2 values)")(
      "scl_ksi_1",
      po::value<std::string>(),
      "comma- or space-separated scl_ksi row 1 values (2 values)")("scl_ksi_2", po::value<std::string>(), "comma- or space-separated scl_ksi row 2 values (2 values)")(
      "scl_ksi_3",
      po::value<std::string>(),
      "comma- or space-separated scl_ksi row 3 values (2 values)")("scl_ksi_4", po::value<std::string>(), "comma- or space-separated scl_ksi row 4 values (2 values)")(
      "scl_ksi_5",
      po::value<std::string>(),
      "comma- or space-separated scl_ksi row 5 values (2 values)")("scl_ksi_6", po::value<std::string>(), "comma- or space-separated scl_ksi row 6 values (2 values)")(
      "scl_K_avg",
      po::value<std::string>(),
      "comma- or space-separated scl_K_avg values (8 values)")("scl_k_sum", po::value<std::string>(), "comma- or space-separated scl_k_sum values (2 values)");

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

inline void validate_config_file_keys(const std::string &filename)
{
  static const char *allowed[] = {
    "phase",     "lat",       "svd_n_cut", "max_iter",  "step",      "flat_file_every", "b3_opt",    "b4_opt",    "b3_zero",   "b4_zero", "lattice",
    "bnl_scl",   "scl_svd",   "scl_h_0",   "scl_h_1",   "scl_h_2",   "scl_a_0",         "scl_a_1",   "scl_a_2",   "scl_a_3",   "scl_a_4", "scl_a_5",
    "scl_ksi_0", "scl_ksi_1", "scl_ksi_2", "scl_ksi_3", "scl_ksi_4", "scl_ksi_5",       "scl_ksi_6", "scl_K_avg", "scl_k_sum", "config",
  };

  std::ifstream in(filename.c_str());
  if (!in)
    return;

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

    bool is_allowed = false;
    for (const char *allowed_key : allowed)
      if (key == allowed_key)
      {
        is_allowed = true;
        break;
      }

    if (!is_allowed)
      throw std::runtime_error("unknown key '" + key + "' at line " + std::to_string(line_no) + " in config file '" + filename + "'");
  }
}

template <std::size_t N> inline void parse_csv_to_array_exact(const std::string &opt_name, const std::string &s, std::array<double, N> &arr)
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

  if (vals.size() != static_cast<std::size_t>(N))
    throw std::runtime_error("option '" + opt_name + "' expects exactly " + std::to_string(static_cast<unsigned long long>(N)) + " values, got " +
                             std::to_string(static_cast<unsigned long long>(vals.size())));

  for (std::size_t i = 0; i < static_cast<std::size_t>(N); ++i)
    arr[i] = vals[i];
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
      validate_config_file_keys(cfg_name);

      std::ifstream config_file(cfg_name);
      if (config_file)
        po::store(po::parse_config_file(config_file, desc), vm);
      else
        std::cerr << "Warning: Could not open config file: " << cfg_name << "\n";
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
  cfg.lat = vm["lat"].as<int>();

  if (cfg.phase != 1 && cfg.phase != 2)
  {
    std::cerr << "Error: phase must be 1 or 2 (got " << cfg.phase << ")\n";
    exit(1);
  }

  if (cfg.phase == 1)
  {
    cfg.scl_ksi = {{
      {{0e0, 0e0}},
      {{1e2, 1e2}},
      {{5e0, 5e0}},
      {{5e0, 5e0}},
      {{5e0, 5e0}},
      {{5e0, 5e0}},
      {{5e0, 5e0}},
    }};
    cfg.scl_a = {{
      {{1e0, 1e0, 1e0, 0e0, 0e0, 0e0}},
      {{1e0, 1e0, 1e0, 0e0, 0e0, 0e0}},
      {{1e0, 1e0, 1e0, 1e0, 0e0, 0e0}},
      {{1e0, 1e0, 1e0, 1e0, 0e0, 0e0}},
      {{1e0, 1e0, 1e0, 1e0, 1e0, 0e0}},
      {{1e0, 1e0, 1e0, 1e0, 1e0, 1e0}},
    }};
    cfg.scl_K_avg = {{1e-3, 1e-3, 1e-3, 1e-3, 5e3, 5e3, 5e3, 5e3}};
  }

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

  if (vm.count("bnl_scl"))
    parse_csv_to_array_exact<6>("bnl_scl", vm["bnl_scl"].as<std::string>(), cfg.bnL_scl);
  if (vm.count("scl_svd"))
    parse_csv_to_array_exact<10>("scl_svd", vm["scl_svd"].as<std::string>(), cfg.scl_svd);

  if (vm.count("scl_h_0"))
    parse_csv_to_array_exact<8>("scl_h_0", vm["scl_h_0"].as<std::string>(), cfg.scl_h[0]);
  if (vm.count("scl_h_1"))
    parse_csv_to_array_exact<8>("scl_h_1", vm["scl_h_1"].as<std::string>(), cfg.scl_h[1]);
  if (vm.count("scl_h_2"))
    parse_csv_to_array_exact<8>("scl_h_2", vm["scl_h_2"].as<std::string>(), cfg.scl_h[2]);

  if (vm.count("scl_a_0"))
    parse_csv_to_array_exact<6>("scl_a_0", vm["scl_a_0"].as<std::string>(), cfg.scl_a[0]);
  if (vm.count("scl_a_1"))
    parse_csv_to_array_exact<6>("scl_a_1", vm["scl_a_1"].as<std::string>(), cfg.scl_a[1]);
  if (vm.count("scl_a_2"))
    parse_csv_to_array_exact<6>("scl_a_2", vm["scl_a_2"].as<std::string>(), cfg.scl_a[2]);
  if (vm.count("scl_a_3"))
    parse_csv_to_array_exact<6>("scl_a_3", vm["scl_a_3"].as<std::string>(), cfg.scl_a[3]);
  if (vm.count("scl_a_4"))
    parse_csv_to_array_exact<6>("scl_a_4", vm["scl_a_4"].as<std::string>(), cfg.scl_a[4]);
  if (vm.count("scl_a_5"))
    parse_csv_to_array_exact<6>("scl_a_5", vm["scl_a_5"].as<std::string>(), cfg.scl_a[5]);

  if (vm.count("scl_ksi_0"))
    parse_csv_to_array_exact<2>("scl_ksi_0", vm["scl_ksi_0"].as<std::string>(), cfg.scl_ksi[0]);
  if (vm.count("scl_ksi_1"))
    parse_csv_to_array_exact<2>("scl_ksi_1", vm["scl_ksi_1"].as<std::string>(), cfg.scl_ksi[1]);
  if (vm.count("scl_ksi_2"))
    parse_csv_to_array_exact<2>("scl_ksi_2", vm["scl_ksi_2"].as<std::string>(), cfg.scl_ksi[2]);
  if (vm.count("scl_ksi_3"))
    parse_csv_to_array_exact<2>("scl_ksi_3", vm["scl_ksi_3"].as<std::string>(), cfg.scl_ksi[3]);
  if (vm.count("scl_ksi_4"))
    parse_csv_to_array_exact<2>("scl_ksi_4", vm["scl_ksi_4"].as<std::string>(), cfg.scl_ksi[4]);
  if (vm.count("scl_ksi_5"))
    parse_csv_to_array_exact<2>("scl_ksi_5", vm["scl_ksi_5"].as<std::string>(), cfg.scl_ksi[5]);
  if (vm.count("scl_ksi_6"))
    parse_csv_to_array_exact<2>("scl_ksi_6", vm["scl_ksi_6"].as<std::string>(), cfg.scl_ksi[6]);

  if (vm.count("scl_K_avg"))
    parse_csv_to_array_exact<8>("scl_K_avg", vm["scl_K_avg"].as<std::string>(), cfg.scl_K_avg);
  if (vm.count("scl_k_sum"))
    parse_csv_to_array_exact<2>("scl_k_sum", vm["scl_k_sum"].as<std::string>(), cfg.scl_k_sum);

  return cfg;
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
  out << "lat=" << cfg.lat << "\n";
  out << "svd_n_cut=" << cfg.svd_n_cut << "\n";
  out << "max_iter=" << cfg.max_iter << "\n";
  out << "step=" << cfg.step << "\n";
  out << "flat_file_every=" << cfg.flat_file_every << "\n";
  out << "b3_opt=" << (cfg.b_3_opt ? "true" : "false") << "\n";
  out << "b4_opt=" << (cfg.b_4_opt ? "true" : "false") << "\n";
  out << "b3_zero=" << (cfg.b_3_zero ? "true" : "false") << "\n";
  out << "b4_zero=" << (cfg.b_4_zero ? "true" : "false") << "\n";
  out << "lattice=" << cfg.lattice_file << "\n";

  out << "bnl_scl=";
  for (size_t i = 0; i < cfg.bnL_scl.size(); ++i)
  {
    if (i)
      out << ",";
    out << cfg.bnL_scl[i];
  }
  out << "\n";

  out << "scl_svd=";
  for (size_t i = 0; i < cfg.scl_svd.size(); ++i)
  {
    if (i)
      out << ",";
    out << cfg.scl_svd[i];
  }
  out << "\n";

  for (size_t r = 0; r < cfg.scl_h.size(); ++r)
  {
    out << "scl_h_" << r << "=";
    for (size_t c = 0; c < cfg.scl_h[r].size(); ++c)
    {
      if (c)
        out << ",";
      out << cfg.scl_h[r][c];
    }
    out << "\n";
  }

  for (size_t r = 0; r < cfg.scl_a.size(); ++r)
  {
    out << "scl_a_" << r << "=";
    for (size_t c = 0; c < cfg.scl_a[r].size(); ++c)
    {
      if (c)
        out << ",";
      out << cfg.scl_a[r][c];
    }
    out << "\n";
  }

  for (size_t r = 0; r < cfg.scl_ksi.size(); ++r)
  {
    out << "scl_ksi_" << r << "=";
    for (size_t c = 0; c < cfg.scl_ksi[r].size(); ++c)
    {
      if (c)
        out << ",";
      out << cfg.scl_ksi[r][c];
    }
    out << "\n";
  }

  out << "scl_K_avg=";
  for (size_t i = 0; i < cfg.scl_K_avg.size(); ++i)
  {
    if (i)
      out << ",";
    out << cfg.scl_K_avg[i];
  }
  out << "\n";

  out << "scl_k_sum=";
  for (size_t i = 0; i < cfg.scl_k_sum.size(); ++i)
  {
    if (i)
      out << ",";
    out << cfg.scl_k_sum[i];
  }
  out << "\n";
}

#endif
