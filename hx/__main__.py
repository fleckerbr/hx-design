import csv
import math
from typing import Union

import toml
from pint import UnitRegistry
from uncertainties import ufloat, umath

from hx import colors, pintutil
from hx.fluids import Coolant


def lmtd_analysis(
    heat_exchanger: dict[str, Union[str, int]],
    hot_coolant: Coolant,
    cold_coolant: Coolant,
) -> None:

    # Define pint utilities
    ureg = UnitRegistry()
    M_ = ureg.Measurement
    Q_ = ureg.Quantity

    # Calculate hot coolant temperatures
    hot_coolant_delta_temperature = M_(
        *pintutil.mtargs(
            hot_coolant.temperature_change(
                "-" + heat_exchanger["energy"],
                heat_exchanger["hot_mass_flow_rate"],
            )
        )
    )
    hot_coolant_inlet_temperature = M_(
        *pintutil.mtargs(heat_exchanger["hot_inlet_temperature"])
    )
    hot_coolant_outlet_temperature = (
        hot_coolant_inlet_temperature.to("kelvin") + hot_coolant_delta_temperature
    ).to("degC")
    hot_coolant_delta_temperature = hot_coolant_delta_temperature.to("delta_degC")

    pintutil.mprint(
        f"{colors.fg.blue}Hot Coolant Ti {colors.fg.darkgrey}::{colors.reset} ",
        hot_coolant_inlet_temperature,
        ".3fP~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Hot Coolant To {colors.fg.darkgrey}::{colors.reset} ",
        hot_coolant_outlet_temperature,
        ".3fP~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Hot Coolant ΔT {colors.fg.darkgrey}::{colors.reset} ",
        hot_coolant_delta_temperature,
        ".3fP~",
    )

    # Calculate cold coolant temperatures
    cold_coolant_delta_temperature = M_(
        *pintutil.mtargs(
            cold_coolant.temperature_change(
                heat_exchanger["energy"],
                heat_exchanger["cold_mass_flow_rate"],
            )
        )
    )
    cold_coolant_inlet_temperature = M_(
        *pintutil.mtargs(heat_exchanger["cold_inlet_temperature"])
    )
    cold_coolant_outlet_temperature = (
        cold_coolant_inlet_temperature.to("kelvin") + cold_coolant_delta_temperature
    ).to("degC")
    cold_coolant_delta_temperature = cold_coolant_delta_temperature.to("delta_degC")

    pintutil.mprint(
        f"{colors.fg.blue}Cold Coolant Ti {colors.fg.darkgrey}::{colors.reset} ",
        cold_coolant_inlet_temperature,
        ".3fP~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Cold Coolant To {colors.fg.darkgrey}::{colors.reset} ",
        cold_coolant_outlet_temperature,
        ".3fP~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Cold Coolant ΔT {colors.fg.darkgrey}::{colors.reset} ",
        cold_coolant_delta_temperature,
        ".3fP~",
    )

    # Calculate log mean temperature difference
    delta_t1 = hot_coolant_inlet_temperature - cold_coolant_inlet_temperature
    delta_t2 = hot_coolant_outlet_temperature - cold_coolant_outlet_temperature
    lmtd_cf = (delta_t1 - delta_t2) / umath.log(
        ufloat(
            (delta_t1 / delta_t2).value.magnitude,
            (delta_t1 / delta_t2).error.magnitude,
        )
    )

    pintutil.mprint(
        f"{colors.fg.blue}LMTD {colors.fg.darkgrey}::{colors.reset} ",
        lmtd_cf,
        ".3fP~",
    )

    # Calculate heat exchanger dimensions
    channel_volume = M_(
        *pintutil.mtargs(heat_exchanger["channel_volume"])
    ).to_root_units()
    plate_width = M_(*pintutil.mtargs(heat_exchanger["plate_width"])).to_root_units()
    plate_height = M_(*pintutil.mtargs(heat_exchanger["plate_height"])).to_root_units()
    plate_spacing = M_(
        *pintutil.mtargs(heat_exchanger["plate_spacing"])
    ).to_root_units()
    plate_thickness = M_(
        *pintutil.mtargs(heat_exchanger["plate_thickness"])
    ).to_root_units()

    channel_area = plate_width * plate_spacing
    plate_surface_area = plate_width * plate_height
    plate_volume = plate_surface_area * plate_thickness
    plate_characteristic_length = plate_volume / plate_surface_area

    pintutil.mprint(
        f"{colors.fg.blue}Channel Area {colors.fg.darkgrey}::{colors.reset} ",
        channel_area,
        ".6fP~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Channel Volume {colors.fg.darkgrey}::{colors.reset} ",
        channel_volume,
        ".6fP~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Plate Volume {colors.fg.darkgrey}::{colors.reset} ",
        plate_volume,
        ".4P~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Plate Characteristic Length {colors.fg.darkgrey}::{colors.reset} ",
        plate_characteristic_length,
        ".4P~",
    )

    # Calculate fluid velocity
    hot_fluid_velocity = (
        M_(*pintutil.mtargs(heat_exchanger["hot_mass_flow_rate"])).to_root_units()
        / M_(*pintutil.mtargs(hot_coolant.density)).to_root_units()
        / channel_area
    )
    cold_fluid_velocity = (
        M_(*pintutil.mtargs(heat_exchanger["cold_mass_flow_rate"])).to_root_units()
        / M_(*pintutil.mtargs(cold_coolant.density)).to_root_units()
        / channel_area
    )

    pintutil.mprint(
        f"{colors.fg.blue}Hot Fluid Velocity {colors.fg.darkgrey}::{colors.reset} ",
        hot_fluid_velocity,
        ".6fP~",
    )
    pintutil.mprint(
        f"{colors.fg.blue}Cold Fluid Velocity {colors.fg.darkgrey}::{colors.reset} ",
        cold_fluid_velocity,
        ".6fP~",
    )

    # Create csv file for data logging
    with open("data/results.csv", "w") as csvfile:
        csvfile.write("P_used, P_req, V_h, V_c, Re_h, Re_c, U [W/(m*degC)]\n")

    # Determine plate count based on defined parameters
    energy = M_(*pintutil.mtargs(heat_exchanger["energy"])).to_root_units()
    for plates in range(1, heat_exchanger["plate_max_count"]):
        split_hot_fluid_velocity = hot_fluid_velocity / plates
        split_cold_fluid_velocity = cold_fluid_velocity / plates

        hot_convective_heat_transfer_coefficient = M_(
            *pintutil.mtargs(
                hot_coolant.convective_coefficient(
                    (0.1381, 0.75, 0.333),
                    f"{split_hot_fluid_velocity:C}",
                    f"{plate_characteristic_length:C}",
                )
            )
        )
        cold_convective_heat_transfer_coefficient = M_(
            *pintutil.mtargs(
                cold_coolant.convective_coefficient(
                    (0.1381, 0.75, 0.333),
                    f"{split_cold_fluid_velocity:C}",
                    f"{plate_characteristic_length:C}",
                )
            )
        )

        rho = M_(*pintutil.mtargs(hot_coolant.density)).to_root_units()
        mu = M_(*pintutil.mtargs(hot_coolant.dynamic_viscosity)).to_root_units()
        hot_reynolds_number = (
            rho * hot_fluid_velocity * plate_characteristic_length / mu
        ).to_root_units()

        rho = M_(*pintutil.mtargs(cold_coolant.density)).to_root_units()
        mu = M_(*pintutil.mtargs(cold_coolant.dynamic_viscosity)).to_root_units()
        cold_reynolds_number = (
            rho * cold_fluid_velocity * plate_characteristic_length / mu
        ).to_root_units()

        ufactor = (
            1 / hot_convective_heat_transfer_coefficient
            + 1 / cold_convective_heat_transfer_coefficient
        ).to("m**2*degC/W")
        ufactor = 1 / ufactor
        surface_area = (energy / plates) / (ufactor * lmtd_cf)
        plates_required = round((surface_area / plate_surface_area).value.magnitude)

        with open("data/results.csv", "a", newline="\n") as csvfile:
            csvfile.write(
                "{p_used}, {p_req}, {v_hot}, {v_cold}, {re_hot:.0f}, {re_cold:.0f}, {u}\n".format(
                    p_used=plates,
                    p_req=plates_required,
                    v_hot=split_hot_fluid_velocity.value.magnitude,
                    v_cold=split_cold_fluid_velocity.value.magnitude,
                    re_hot=hot_reynolds_number.value.magnitude,
                    re_cold=cold_reynolds_number.value.magnitude,
                    u=ufactor.value.magnitude,
                )
            )

        if plates >= plates_required:
            pintutil.mprint(
                f"{colors.fg.blue}Overall Heat Transfer Coefficient {colors.fg.darkgrey}::{colors.reset} ",
                ufactor,
                ".6fP~",
            )
            print(
                f"{colors.fg.blue}Hot Coolant Re {colors.fg.darkgrey}::{colors.reset}",
                hot_reynolds_number.value.magnitude,
            )
            print(
                f"{colors.fg.blue}Cold Coolant Re {colors.fg.darkgrey}::{colors.reset}",
                cold_reynolds_number.value.magnitude,
            )
            print(
                f"{colors.fg.blue}Plates Required {colors.fg.darkgrey}::{colors.reset}",
                plates_required,
            )
            break
    else:
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
