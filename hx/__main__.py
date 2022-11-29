import math
from typing import Union

import toml
from pint import UnitRegistry
from uncertainties import ufloat, umath

from hx import colors, pintil
from hx.fluids import Coolant


def lmtd_analysis(
    heat_exchanger: dict[str, Union[str, int, float]],
    hot_coolant: Coolant,
    cold_coolant: Coolant,
) -> None:

    # Define pint utilities
    ureg = UnitRegistry()
    M_ = ureg.Measurement
    Q_ = ureg.Quantity

    # Calculate hot coolant temperatures
    hot_coolant_delta_temperature = M_(
        *pintil.mtargs(
            hot_coolant.temperature_change(
                "-" + heat_exchanger["energy"],
                heat_exchanger["hot_mass_flow_rate"],
            )
        )
    )
    hot_coolant_inlet_temperature = M_(
        *pintil.mtargs(heat_exchanger["hot_inlet_temperature"])
    )
    hot_coolant_outlet_temperature = (
        hot_coolant_inlet_temperature.to("kelvin") + hot_coolant_delta_temperature
    ).to("degC")
    hot_coolant_delta_temperature = hot_coolant_delta_temperature.to("delta_degC")

    # Calculate cold coolant temperatures
    cold_coolant_delta_temperature = M_(
        *pintil.mtargs(
            cold_coolant.temperature_change(
                heat_exchanger["energy"],
                heat_exchanger["cold_mass_flow_rate"],
            )
        )
    )
    cold_coolant_inlet_temperature = M_(
        *pintil.mtargs(heat_exchanger["cold_inlet_temperature"])
    )
    cold_coolant_outlet_temperature = (
        cold_coolant_inlet_temperature.to("kelvin") + cold_coolant_delta_temperature
    ).to("degC")
    cold_coolant_delta_temperature = cold_coolant_delta_temperature.to("delta_degC")

    # Calculate log mean temperature difference
    delta_t1 = hot_coolant_inlet_temperature - cold_coolant_inlet_temperature
    delta_t2 = hot_coolant_outlet_temperature - cold_coolant_outlet_temperature
    lmtd_cf = (delta_t1 - delta_t2) / umath.log(
        ufloat(
            (delta_t1 / delta_t2).value.magnitude,
            (delta_t1 / delta_t2).error.magnitude,
        )
    )

    # Calculate heat exchanger dimensions
    channel_volume = M_(
        *pintil.mtargs(heat_exchanger["channel_volume"])
    ).to_root_units()
    plate_width = M_(*pintil.mtargs(heat_exchanger["plate_width"])).to_root_units()
    plate_height = M_(*pintil.mtargs(heat_exchanger["plate_height"])).to_root_units()

    plate_surface_area = plate_width * plate_height
    channel_width = channel_volume / plate_surface_area
    channel_area = plate_width * channel_width
    hydraulic_diameter = 4 * channel_area / (2 * plate_width + 2 * channel_width)

    # Calculate fluid velocity
    max_hot_fluid_velocity = (
        M_(*pintil.mtargs(heat_exchanger["hot_mass_flow_rate"])).to_root_units()
        / M_(*pintil.mtargs(hot_coolant.density)).to_root_units()
        / channel_area
    )
    max_cold_fluid_velocity = (
        M_(*pintil.mtargs(heat_exchanger["cold_mass_flow_rate"])).to_root_units()
        / M_(*pintil.mtargs(cold_coolant.density)).to_root_units()
        / channel_area
    )

    # Create csv file for data logging
    with open("data/results.csv", "w") as csvfile:
        csvfile.write(
            "P_used, P_req, V_hot, V_cold, Re_hot, Re_cold, Nu_hot, Nu_cold, h_hot [W/(m*°C)], h_cold [W/(m*°C)], U [W/(m*°C)]\n",
        )

    # Preload data from parameters file
    nusselt_coefficient: float = heat_exchanger["nusselt_coefficient"]
    nusselt_coefficient: float = heat_exchanger["nusselt_coefficient"]
    energy = M_(*pintil.mtargs(heat_exchanger["energy"])).to_root_units()
    nusselt_coefficient: float = heat_exchanger["nusselt_coefficient"]
    reynolds_exponent: float = heat_exchanger["reynolds_exponent"]
    prandtl_exponent: float = heat_exchanger["prandtl_exponent"]

    # Determine plate count based on defined parameters
    solution_identified = False
    for plates_used in range(1, heat_exchanger["plate_max_count"]):
        hot_fluid_velocity = max_hot_fluid_velocity / plates_used
        cold_fluid_velocity = max_cold_fluid_velocity / plates_used

        hot_reynolds_number = M_(
            *pintil.mtargs(
                hot_coolant.reynolds_number(
                    f"{hot_fluid_velocity:C}",
                    f"{hydraulic_diameter:C}",
                )
            )
        )
        cold_reynolds_number = M_(
            *pintil.mtargs(
                cold_coolant.reynolds_number(
                    f"{cold_fluid_velocity:C}",
                    f"{hydraulic_diameter:C}",
                )
            )
        )

        hot_nusselt_number = M_(
            *pintil.mtargs(
                hot_coolant.nusselt_number(
                    f"{hot_reynolds_number:C}",
                    nusselt_coefficient,
                    reynolds_exponent,
                    prandtl_exponent,
                )
            )
        )
        cold_nusselt_number = M_(
            *pintil.mtargs(
                cold_coolant.nusselt_number(
                    f"{cold_reynolds_number:C}",
                    nusselt_coefficient,
                    reynolds_exponent,
                    prandtl_exponent,
                )
            )
        )

        hot_convective_heat_transfer_coefficient = M_(
            *pintil.mtargs(
                hot_coolant.convective_coefficient(
                    f"{hot_nusselt_number:C}",
                    f"{hydraulic_diameter:C}",
                )
            )
        )
        cold_convective_heat_transfer_coefficient = M_(
            *pintil.mtargs(
                cold_coolant.convective_coefficient(
                    f"{cold_nusselt_number:C}",
                    f"{hydraulic_diameter:C}",
                )
            )
        )

        u_value = (
            1 / hot_convective_heat_transfer_coefficient
            + 1 / cold_convective_heat_transfer_coefficient
        ).to("m**2*degC/W")
        u_value = 1 / u_value
        required_surface_area = (energy / plates_used) / (
            u_value.to_root_units() * lmtd_cf.to_root_units()
        )
        plates_required = required_surface_area / plate_surface_area
        max_plates_required = math.ceil(
            plates_required.value.magnitude + plates_required.error.magnitude
        )

        with open("data/results.csv", "a", newline="\n") as csvfile:
            csvfile.write(
                "{p_used}, {p_req}, {v_hot}, {v_cold}, {re_hot:.0f}, {re_cold:.0f}, {nu_hot:.0f}, {nu_cold:.0f}, {h_hot:.0f}, {h_cold:.0f}, {u}\n".format(
                    p_used=plates_used,
                    p_req=plates_required.magnitude,
                    v_hot=hot_fluid_velocity.magnitude,
                    v_cold=cold_fluid_velocity.magnitude,
                    re_hot=hot_reynolds_number.magnitude,
                    re_cold=cold_reynolds_number.magnitude,
                    nu_hot=hot_nusselt_number.magnitude,
                    nu_cold=cold_nusselt_number.magnitude,
                    h_hot=hot_convective_heat_transfer_coefficient.magnitude,
                    h_cold=cold_convective_heat_transfer_coefficient.magnitude,
                    u=u_value.magnitude,
                )
            )

        if plates_used >= max_plates_required and not solution_identified:
            solution_identified = True
            print(
                f"\n{colors.underline}Hot Coolant :: {hot_coolant.name}{colors.reset}",
            )
            pintil.mprint(
                f"{colors.fg.blue}Inlet Temperature {colors.fg.lightgreen}[Ti] {colors.fg.darkgrey}::{colors.reset} ",
                hot_coolant_inlet_temperature,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Outlet Temperature {colors.fg.lightgreen}[To] {colors.fg.darkgrey}::{colors.reset} ",
                hot_coolant_outlet_temperature,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Delta Temperature {colors.fg.lightgreen}[ΔT] {colors.fg.darkgrey}::{colors.reset} ",
                hot_coolant_delta_temperature,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Max Velocity {colors.fg.lightgreen}[u_max] {colors.fg.darkgrey}::{colors.reset} ",
                max_hot_fluid_velocity,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Channel Velocity {colors.fg.lightgreen}[u] {colors.fg.darkgrey}::{colors.reset} ",
                hot_fluid_velocity,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Reynolds Number {colors.fg.lightgreen}[Re] {colors.fg.darkgrey}::{colors.reset}",
                hot_reynolds_number,
                ".0fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Nusselt Number {colors.fg.lightgreen}[Nu] {colors.fg.darkgrey}::{colors.reset}",
                hot_nusselt_number,
                ".0fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Convective Coefficient {colors.fg.lightgreen}[h] {colors.fg.darkgrey}::{colors.reset} ",
                hot_convective_heat_transfer_coefficient,
                ".0fP~",
            )
            print(
                f"\n{colors.underline}Cold Coolant :: {cold_coolant.name}{colors.reset}",
            )
            pintil.mprint(
                f"{colors.fg.blue}Inlet Temperature {colors.fg.lightgreen}[Ti] {colors.fg.darkgrey}::{colors.reset} ",
                cold_coolant_inlet_temperature,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Outlet Temperature {colors.fg.lightgreen}[To] {colors.fg.darkgrey}::{colors.reset} ",
                cold_coolant_outlet_temperature,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Delta Temperature {colors.fg.lightgreen}[ΔT] {colors.fg.darkgrey}::{colors.reset} ",
                cold_coolant_delta_temperature,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Max Velocity {colors.fg.lightgreen}[u_max] {colors.fg.darkgrey}::{colors.reset} ",
                max_cold_fluid_velocity,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Channel Velocity {colors.fg.lightgreen}[u] {colors.fg.darkgrey}::{colors.reset} ",
                cold_fluid_velocity,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Reynolds Number {colors.fg.lightgreen}[Re] {colors.fg.darkgrey}::{colors.reset}",
                cold_reynolds_number,
                ".0fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Nusselt Number {colors.fg.lightgreen}[Nu] {colors.fg.darkgrey}::{colors.reset}",
                cold_nusselt_number,
                ".0fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Convective Coefficient {colors.fg.lightgreen}[h] {colors.fg.darkgrey}::{colors.reset} ",
                cold_convective_heat_transfer_coefficient,
                ".0fP~",
            )
            print(
                f"\n{colors.underline}Heat Exchanger :: {heat_exchanger.get('name')}{colors.reset}",
            )
            pintil.mprint(
                f"{colors.fg.blue}LMTD {colors.fg.darkgrey}::{colors.reset} ",
                lmtd_cf,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Channel Width {colors.fg.lightgreen}[V] {colors.fg.darkgrey}::{colors.reset} ",
                channel_width,
                ".6fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Channel Volume {colors.fg.lightgreen}[V] {colors.fg.darkgrey}::{colors.reset} ",
                channel_volume,
                ".6fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Channel Area {colors.fg.lightgreen}[Ac] {colors.fg.darkgrey}::{colors.reset} ",
                channel_area,
                ".6fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Hydraulic Diameter {colors.fg.lightgreen}[Dh] {colors.fg.darkgrey}::{colors.reset} ",
                hydraulic_diameter,
                ".6fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}U-value {colors.fg.lightgreen}[U] {colors.fg.darkgrey}::{colors.reset} ",
                u_value,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Required Surface Area {colors.fg.lightgreen}[A] {colors.fg.darkgrey}::{colors.reset}",
                required_surface_area,
                ".3fP~",
            )
            pintil.mprint(
                f"{colors.fg.blue}Plates Required {colors.fg.darkgrey}::{colors.reset}",
                plates_required,
                ".3fP~",
            )
            print(
                f"{colors.fg.blue}Plates Used {colors.fg.darkgrey}::{colors.reset}",
                plates_used,
            )

    if not solution_identified:
        print(
            f"\n{colors.underline}Hot Coolant :: {hot_coolant.name}{colors.reset}",
        )
        pintil.mprint(
            f"{colors.fg.blue}Inlet Temperature {colors.fg.lightgreen}[Ti] {colors.fg.darkgrey}::{colors.reset} ",
            hot_coolant_inlet_temperature,
            ".3fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Outlet Temperature {colors.fg.lightgreen}[To] {colors.fg.darkgrey}::{colors.reset} ",
            hot_coolant_outlet_temperature,
            ".3fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Delta Temperature {colors.fg.lightgreen}[ΔT] {colors.fg.darkgrey}::{colors.reset} ",
            hot_coolant_delta_temperature,
            ".3fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Max Velocity {colors.fg.lightgreen}[u_max] {colors.fg.darkgrey}::{colors.reset} ",
            max_hot_fluid_velocity,
            ".3fP~",
        )
        print(
            f"\n{colors.underline}Cold Coolant :: {cold_coolant.name}{colors.reset}",
        )
        pintil.mprint(
            f"{colors.fg.blue}Inlet Temperature {colors.fg.lightgreen}[Ti] {colors.fg.darkgrey}::{colors.reset} ",
            cold_coolant_inlet_temperature,
            ".3fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Outlet Temperature {colors.fg.lightgreen}[To] {colors.fg.darkgrey}::{colors.reset} ",
            cold_coolant_outlet_temperature,
            ".3fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Delta Temperature {colors.fg.lightgreen}[ΔT] {colors.fg.darkgrey}::{colors.reset} ",
            cold_coolant_delta_temperature,
            ".3fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Max Velocity {colors.fg.lightgreen}[u_max] {colors.fg.darkgrey}::{colors.reset} ",
            max_cold_fluid_velocity,
            ".3fP~",
        )
        print(
            f"{colors.underline}Heat Exchanger :: {heat_exchanger.get('name')}{colors.reset}",
        )
        pintil.mprint(
            f"{colors.fg.blue}LMTD {colors.fg.darkgrey}::{colors.reset} ",
            lmtd_cf,
            ".3fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Channel Width {colors.fg.lightgreen}[V] {colors.fg.darkgrey}::{colors.reset} ",
            channel_width,
            ".6fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Channel Volume {colors.fg.lightgreen}[V] {colors.fg.darkgrey}::{colors.reset} ",
            channel_volume,
            ".6fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Channel Area {colors.fg.lightgreen}[Ac] {colors.fg.darkgrey}::{colors.reset} ",
            channel_area,
            ".6fP~",
        )
        pintil.mprint(
            f"{colors.fg.blue}Channel Hydraulic Diameter {colors.fg.lightgreen}[Dh] {colors.fg.darkgrey}::{colors.reset} ",
            hydraulic_diameter,
            ".4P~",
        )
        print(
            f"{colors.fg.blue}Plates Required {colors.fg.darkgrey}::{colors.reset}",
            f"{colors.fg.red}No Solution{colors.reset}",
        )


if __name__ == "__main__":
    # Load heat exchanger parameters
    parameters = toml.load("data/parameters.toml")
    coolant_params = (
        "density",
        "dynamic_viscosity",
        "name",
        "prandtl_number",
        "specific_heat",
        "thermal_conductivity",
    )
    hot_coolant_params = {
        k: v for k, v in parameters["hot-coolant"].items() if k in coolant_params
    }
    cold_coolant_params = {
        k: v for k, v in parameters["cold-coolant"].items() if k in coolant_params
    }

    # Run heat transfer analysis with the specified coolants
    dielectric_fluid = Coolant(**hot_coolant_params)
    water = Coolant(**cold_coolant_params)
    lmtd_analysis(parameters["plate-heat-exchanger"], dielectric_fluid, water)
    print()
