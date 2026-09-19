#include "../include/cooling_correlations.h"
#include "../include/heat_transfer.h"
#include "../include/humidair.h"
#include "../include/transport.h"

#include <iomanip>
#include <iostream>
#include <vector>

using namespace combaero;

// Cooling example: simple convective + effusion-cooled liner estimate.

int main() {
    using combaero::cooling::cooled_wall_heat_flux;

    std::cout << std::fixed << std::setprecision(3);

    const std::vector<double> x_air = dry_air();

    // Representative combustor conditions
    const double p = 15.0e5;        // Pa
    const double t_hot = 1850.0;    // K
    const double t_cool = 760.0;    // K
    const double v_hot = 28.0;      // m/s
    const double v_cool = 42.0;     // m/s
    const double d_hot = 0.085;     // m (hot-side hydraulic diameter)
    const double d_cool = 0.012;    // m (cooling channel hydraulic diameter)
    const double t_wall = 0.0025;   // m
    const double k_wall = 19.0;     // W/(m*K)

    const double re_hot = reynolds(t_hot, p, x_air, v_hot, d_hot);
    const double pr_hot = prandtl(t_hot, p, x_air);
    const double k_hot = thermal_conductivity(t_hot, p, x_air);
    const double nu_hot = nusselt_dittus_boelter(re_hot, pr_hot, false);
    const double h_hot = htc_from_nusselt(nu_hot, k_hot, d_hot);

    const double re_cool = reynolds(t_cool, p, x_air, v_cool, d_cool);
    const double pr_cool = prandtl(t_cool, p, x_air);
    const double k_cool = thermal_conductivity(t_cool, p, x_air);
    const double nu_cool = nusselt_dittus_boelter(re_cool, pr_cool, true);
    const double h_cool = htc_from_nusselt(nu_cool, k_cool, d_cool);

    std::cout << "==============================================================\n";
    std::cout << "Simple Cooling Example (C++)\n";
    std::cout << "==============================================================\n";
    std::cout << "Hot side:  Re=" << re_hot << ", h=" << h_hot << " W/(m^2*K)\n";
    std::cout << "Cool side: Re=" << re_cool << ", h=" << h_cool << " W/(m^2*K)\n\n";

    std::cout << "eta_eff   q_cooled [kW/m^2]   q_uncooled [kW/m^2]\n";
    std::cout << "---------------------------------------------------------\n";

    // The effusion effectiveness correlation this example used was removed in
    // 0.7.0: its momentum flux ratio was defined as M^2 * DR where the
    // definition is M^2 / DR, and its form could not be traced to the cited
    // source. See issue #339. Until a provenanced replacement lands, the
    // example sweeps effectiveness directly rather than predicting it.
    const double q_uncooled =
        cooled_wall_heat_flux(t_hot, t_cool, h_hot, h_cool, 0.0, t_wall, k_wall) / 1000.0;

    for (double eta = 0.0; eta <= 0.6; eta += 0.1) {
        const double q_cooled =
            cooled_wall_heat_flux(t_hot, t_cool, h_hot, h_cool, eta, t_wall, k_wall) / 1000.0;

        std::cout << std::setw(7) << eta << "   "
                  << std::setw(15) << q_cooled << "   "
                  << std::setw(16) << q_uncooled << "\n";
    }

    return 0;
}
