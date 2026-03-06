#include "registry.h"
#include <stdexcept>

namespace combaero {

Registry::Registry() {
  friction_models_["haaland"] = solver::friction_and_jacobian_haaland;
  friction_models_["serghides"] = solver::friction_and_jacobian_serghides;
  friction_models_["colebrook"] = solver::friction_and_jacobian_colebrook;
  friction_models_["petukhov"] = [](double Re, double) {
    return solver::friction_and_jacobian_petukhov(Re);
  };
}

Registry::FrictionFunc
Registry::get_friction_model(const std::string &name) const {
  auto it = friction_models_.find(name);
  if (it != friction_models_.end()) {
    return it->second;
  }
  throw std::invalid_argument("Unknown friction model registry_id: " + name);
}

bool Registry::has_friction_model(const std::string &name) const {
  return friction_models_.find(name) != friction_models_.end();
}

std::vector<std::string> Registry::available_friction_models() const {
  std::vector<std::string> models;
  models.reserve(friction_models_.size());
  for (const auto &pair : friction_models_) {
    models.push_back(pair.first);
  }
  return models;
}

} // namespace combaero
