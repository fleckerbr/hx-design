# Heat Exchanger Design Problem
This repository contains calculations for a [WCR A216B](http://wcrhx.com/wp-content/uploads/2022/01/WCR-A216B-Spec-Sheet-FINAL.pdf) plate heat exchanger from Sondex.
Since the given spec sheet did not contain information on the spacing between plates the [Sondex Free Flow Heat Exchanger Manual](https://sme-llc.com/dev/wp-content/uploads/2015/12/FREE-FLOW-2007.pdf) was referenced to find a reasonable value.

Calculations are based on the Log Mean Temperature Difference method for calculating heat transfer through a heat exchanger.
The parameters used in the calculations are stored in the [parameters.toml](data/parameters.toml) file.
