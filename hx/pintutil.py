from pint import Measurement


def mtargs(__arg: str, /) -> tuple:
    """Convert string into tuple of arguments for pint Measurement.

    Args:
        __arg: Formatted string for pint.

    Returns:
        tuple: (value, error, units)
    """

    value, units = __arg.rsplit(" ", 1)
    try:
        return float(value), 0.0, units
    except ValueError:
        if "e+" in value:
            value, power = value.split("e+")
            power = "e+" + power
        elif "e-" in value:
            value, power = value.split("e-")
            power = "e-" + power
        else:
            power = "e+0"
        try:
            value, uncertainty = value.split("+/-", 1)
        except ValueError:
            value, uncertainty = value.split("±", 1)
        return (
            float(value.strip("( ") + f"{power}"),
            float(uncertainty.strip(" )") + f"{power}"),
            units,
        )


def mprint(text: str, measurement: Measurement, fmt: str = "P~") -> None:
    """Print pint Measurement with Quantity formatting codes.

    Args:
        text: Text to prefix measurement.
        measurement: Measurement object from pint.
        fmt: Quantity formatting code. Defaults to "P~".
    """

    _, units = f"{measurement.value:{fmt}}".rsplit(" ", 1)
    value, _ = f"{measurement:{fmt.replace('~','')}}".rsplit(" ", 1)
    print(f"{text}{value} {units}")
