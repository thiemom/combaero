#ifndef COMBAERO_REGISTRY_H
#define COMBAERO_REGISTRY_H

#include "solver_interface.h"
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace combaero {

class Registry {
public:
  using FrictionFunc =
      std::function<solver::CorrelationResult<std::tuple<double, double>>(
          double Re, double e_D)>;

  // Get the global singleton instance
  static Registry &get() {
    static Registry instance;
    return instance;
  }

  // Retrieve a friction function by name
  FrictionFunc get_friction_model(const std::string &name) const;

  // Check if name exists
  bool has_friction_model(const std::string &name) const;

  // Get a list of all available friction models
  std::vector<std::string> available_friction_models() const;

private:
  Registry(); // Private constructor builds the mappings

  std::map<std::string, FrictionFunc> friction_models_;
};

} // namespace combaero

#endif // COMBAERO_REGISTRY_H
