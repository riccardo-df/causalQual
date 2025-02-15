# Causal Inference for Qualitative Outcomes  <a href="https://riccardo-df.github.io/causalQual/"><img src="man/figures/causalQual.svg" align="right" height="130" /></a>

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT) [![CRAN](https://www.r-pkg.org/badges/version/causalQual)](https://CRAN.R-project.org/package=causalQual) [![Downloads](https://cranlogs.r-pkg.org/badges/causalQual)](https://CRAN.R-project.org/package=causalQual)

## Overview

The `causalQual` package provides a suite of tools for estimating causal effects when the outcome of interest is qualitative - i.e., multinomial or ordered. Standard causal inference methods such as instrumental variables (IV), regression discontinuity (RD), and difference-in-differences (DiD) are typically designed for numeric outcomes. Their direct application to qualitative outcomes leads to ill-defined estimands, rendering results arbitrary and uninterpretable.

This package implements the framework introduced in Di Francesco and Mellace (2025), shifting the focus to well-defined and interpretable estimands that quantify how treatment affects the probability distribution over outcome categories. The methods remain compatible with conventional research designs, ensuring ease of implementation for applied researchers.

------------------------------------------------------------------------

## Why Use `causalQual`?

| Feature                             | Benefit                                                                                                                                        |
|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| **Avoids misleading conclusions**   | Conventional estimands are often undefined or depend on arbitrary outcome coding. `causalQual` targets interpretable and meaningful estimands. |
| **Provides well-defined estimands** | Instead of relying on average effects, `causalQual` models how treatment shifts probabilities over outcome categories.                         |
| **Extensible and open-source**      | Actively developed with planned support for staggered adoption, fuzzy regression discontinuity, and more.                                      |

------------------------------------------------------------------------

## Key Features

✔ **Implements causal inference methods for qualitative outcomes**\
✔ **Supports selection-on-observables, IV, RD, DiD**\
✔ **Open-source and extensible**

------------------------------------------------------------------------

## Installation

To install the latest stable release from CRAN, run:

```         
install.packages("causalQual")
```

To install the latest development version from GitHub, run:

```         
devtools::install_github("riccardo-df/causalQual")
```

------------------------------------------------------------------------

## Functionality

| Identification strategy       | Function           |
|-------------------------------|--------------------|
| **Selection-on-Observables**  | `causalQual_soo()` |
| **Instrumental Variables**    | `causalQual_iv()`  |
| **Regression Discontinuity**  | `causalQual_rd()`  |
| **Difference-in-Differences** | `causalQual_did()` |

------------------------------------------------------------------------

## Contributing

We welcome contributions! If you encounter issues, have feature requests, or want to contribute to the package, please follow the guidelines below:

📌 **Report an Issue:** If you encounter a bug or have a suggestion, please open an issue on GitHub:\
[Submit an Issue](https://github.com/riccardo-df/causalQual/issues)

📌 **Contribute Code:** We encourage contributions via pull requests. Before submitting, please:\
1. Fork the repository and create a new branch.\
2. Ensure that your code follows the existing style and documentation conventions.\
3. Run tests and check for package integrity.\
4. Submit a pull request with a clear description of your changes.

📌 **Feature Requests:** If you have ideas for new features or extensions, feel free to discuss them by opening an issue.

------------------------------------------------------------------------

## Citation

If you use `causalQual` in your research, please cite the corresponding paper:

> **Author(s).** *Title of Paper.* arXiv, 2025
