#ifndef GUSTAVO_OPTIONS_TERMS_H
#define GUSTAVO_OPTIONS_TERMS_H

// FUTURE: Per-term weight override infrastructure.
//
// Currently, all Lie term weights are driven by Config's scl_h, scl_a, scl_ksi, scl_K_avg arrays.
// This file is a placeholder for future fine-grained per-term control:
//
//   - Build a stable term catalog (key = "group:label:i,j,k,l,m")
//   - Load per-term overrides from JSON (e.g., term_weights.json per run)
//   - Apply overrides on top of Config defaults
//
// This keeps the current behavior unchanged while preparing the foundation for individual
// term weight tuning via Python sweep script or manual config.

#endif
