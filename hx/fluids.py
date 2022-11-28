from dataclasses import dataclass

from pint import Measurement

from hx import pintutil


@dataclass
class Coolant:
    density: str
    dynamic_viscosity: str
    prandtl_number: float
    specific_heat: str
    thermal_conductivity: str

    def temperature_change(self, energy_transfer: str, mass_flow_rate: str) -> str:
        M_ = Measurement

        temperature = (
            M_(*pintutil.mtargs(energy_transfer)).to_root_units()
            / M_(*pintutil.mtargs(mass_flow_rate)).to_root_units()
            / M_(*pintutil.mtargs(self.specific_heat)).to_root_units()
        ).to("kelvin")
        return f"{temperature:C}"

    def convective_coefficient(
        self,
        nu_params: tuple[float],
        fluid_velocity: str,
        characteristic_length: str,
    ) -> str:
        M_ = Measurement

        c, m, n = nu_params
        rho = M_(*pintutil.mtargs(self.density)).to_root_units()
        mu = M_(*pintutil.mtargs(self.dynamic_viscosity)).to_root_units()
        lc = M_(*pintutil.mtargs(characteristic_length)).to_root_units()
        k = M_(*pintutil.mtargs(self.thermal_conductivity)).to_root_units()
        u = M_(*pintutil.mtargs(fluid_velocity)).to_root_units()

        reynolds_number = rho * u * lc / mu
        nusselt_number = c * reynolds_number**m * self.prandtl_number**n
        htcoeff = (nusselt_number * k / lc).to("W/m**2/degC")

        return f"{htcoeff:C}"
