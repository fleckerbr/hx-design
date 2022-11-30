# Heat Exchanger Design Problem
This repository contains calculations for a [WCR A216B](http://wcrhx.com/wp-content/uploads/2022/01/WCR-A216B-Spec-Sheet-FINAL.pdf) plate heat exchanger from Sondex.

Calculations are based on the Log Mean Temperature Difference method for calculating heat transfer through a heat exchanger.
The parameters used in the calculations are stored in the [parameters.toml](data/parameters.toml) file.

The friction factor for the plates in the Heat Exchanger were calculated using the [Correlation of Mulley](https://powderprocess.net/Tools_html/Thermodynamics/Plate_Heat_Exchanger_Pressure_Drop_Calculation.html) with an assumed angle of corregation.
