# Causal inference for qualitative outcomes <a href="https://riccardo-df.github.io/causalQual/"><img src="man/figures/logo.svg" align="right" height="130"/></a>

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT) [![CRAN](https://www.r-pkg.org/badges/version/causalQual)](https://CRAN.R-project.org/package=causalQual) [![Downloads](https://cranlogs.r-pkg.org/badges/causalQual)](https://CRAN.R-project.org/package=causalQual)

The `causalQual` package provides a suite of tools for estimating causal effects when the outcome targeted by the treatment is qualitative - i.e., multinomial or ordered.

Standard causal inference methods such as instrumental variables (IV), regression discontinuity (RD), and difference-in-differences (DiD) are typically designed for numeric outcomes. Their direct application to qualitative outcomes leads to ill-defined estimands, rendering results arbitrary and uninterpretable.

This package implements the framework introduced in Di Francesco and Mellace (2025), shifting the focus to well-defined and interpretable estimands that quantify how treatment affects the probability distribution over outcome categories. The methods remain compatible with conventional research designs, ensuring ease of implementation for applied researchers.

------------------------------------------------------------------------

## Why use `causalQual`?

| Feature | Benefit |
|--------------------|----------------------------------------------------|
| **Avoids misleading conclusions** | Conventional estimands are often undefined or depend on arbitrary outcome coding. `causalQual` targets interpretable and meaningful estimands. |
| **Provides well-defined estimands** | Instead of relying on average effects, `causalQual` models how treatment shifts probabilities over outcome categories. |
| **Wide applicability** | Supports selection-on-observables, IV, RD, and DiD. |
| **Extensible and open-source** | Actively developed with planned support for staggered adoption, fuzzy regression discontinuity, and more. |

------------------------------------------------------------------------

## 🚀 Installation

To install the latest stable release from CRAN, run:

```         
install.packages("causalQual")
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```         
devtools::install_github("riccardo-df/causalQual")
```

------------------------------------------------------------------------

## Contributing

We welcome contributions! If you encounter issues, have feature requests, or want to contribute to the package, please follow the guidelines below.

📌 **Report an issue:** If you encounter a bug or have a suggestion, please open an issue on GitHub: [Submit an issue](https://github.com/riccardo-df/causalQual/issues)

📌 **Contribute code:** We encourage contributions via pull requests. Before submitting, please: 1. Fork the repository and create a new branch. 2. Ensure that your code follows the existing style and documentation conventions. 3. Run tests and check for package integrity. 4. Submit a pull request with a clear description of your changes.

📌 **Feature requests:** If you have ideas for new features or extensions, feel free to discuss them by opening an issue.

------------------------------------------------------------------------

## Citation

If you use `causalQual` in your research, please cite the corresponding paper:

> **Author(s).** *Title of Paper.* arXiv, 2025
