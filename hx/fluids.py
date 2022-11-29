from dataclasses import dataclass

from pint import Measurement

from hx import pintil


@dataclass
class Coolant:
    density: str
    dynamic_viscosity: str
    name: str
    prandtl_number: str
    specific_heat: str
    thermal_conductivity: str

    def temperature_change(
        self,
        energy_transfer: str,
        mass_flow_rate: str,
    ) -> str:
        M_ = Measurement

        delta_temperature = (
            M_(*pintil.mtargs(energy_transfer)).to_root_units()
            / M_(*pintil.mtargs(mass_flow_rate)).to_root_units()
            / M_(*pintil.mtargs(self.specific_heat)).to_root_units()
        ).to("kelvin")
        return f"{delta_temperature:C}"

    def reynolds_number(
        self,
        fluid_velocity: str,
        hydraulic_diameter: str,
    ) -> str:
        M_ = Measurement

        u = M_(*pintil.mtargs(fluid_velocity)).to_root_units()
        dh = M_(*pintil.mtargs(hydraulic_diameter)).to_root_units()
        rho = M_(*pintil.mtargs(self.density)).to_root_units()
        mu = M_(*pintil.mtargs(self.dynamic_viscosity)).to_root_units()
        re = rho * u * dh / mu

        return f"{re:C}"

    def nusselt_number(
        self,
        reynolds_number: str,
        nusselt_coefficient: float,
        reynolds_exponent: float,
        prandtl_exponent: float,
    ) -> str:
        M_ = Measurement

        c, m, n = (nusselt_coefficient, reynolds_exponent, prandtl_exponent)
        re = M_(*pintil.mtargs(reynolds_number)).to_root_units()
        pr = M_(*pintil.mtargs(self.prandtl_number)).to_root_units()

        nu = c * re**m * pr**n

        return f"{nu:C}"

    def convective_coefficient(
        self,
        nusselt_number: str,
        hydraulic_diameter: str,
    ) -> str:
        M_ = Measurement

        nu = M_(*pintil.mtargs(nusselt_number)).to_root_units()
        dh = M_(*pintil.mtargs(hydraulic_diameter)).to_root_units()
        k = M_(*pintil.mtargs(self.thermal_conductivity)).to_root_units()

        h = (nu * k / dh).to("W/m**2/degC")

        return f"{h:C}"
