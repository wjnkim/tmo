# TMO: Thresholding Multiple Outcomes

Version 1.0.0. TMO estimates standard errors using auxiliary outcomes to measure correlation across locations.

## Install

Run in Stata:

    ssc install gtools
    net install tmo, from("https://raw.githubusercontent.com/wjnkim/tmo/main/src") replace
    help tmo

If Stata cannot access GitHub directly, download this repository and run
`net install tmo, from("path/to/extracted/repository/src") replace` using the local
directory containing `tmo.pkg`.

Stata installation includes the compiled Java kernel. The repository also contains its source and build scripts. Java is optional at runtime; TMO uses its Mata implementation when the kernel is unavailable. Stata 13 or later is required. The helpfile lists the supported regression commands and options.

## Optional dependencies

The scpc and scpc_cmd() options require the TMO-compatible SCPC companion. Install the tested version with:

    net install scpc, from("https://raw.githubusercontent.com/wjnkim/tmo/main/scpc_tmo")

If Stata cannot access GitHub directly, download https://github.com/wjnkim/tmo and run
`net install scpc, from("path/to/extracted/repository/scpc_tmo")`.

This companion is based on Ulrich K. Müller's MIT-licensed SCPC package (https://github.com/ukmueller/SCPC). TMO checks the installed SCPC source before using it. If another SCPC version is already installed, decide whether to replace it before installing the companion.

The geodist package is needed only for legacy distance dyads under a custom threshold override:

    ssc install geodist

The reghdfe, ivreghdfe, and ivreg2 packages are needed only for specifications that use those estimators. See help tmo for the versions tested with TMO.

## Paper specification

Automatic selection uses the grid 0.1, 0.15, ..., 1 and, as in the paper, keeps at most 5% of eligible off-diagonal pairs by default; `maxeligiblepct(100)` removes the restriction. The helpfile defines protected and eligible pairs.

## Files

- src/tmo.ado: Stata command and Mata implementation.
- src/tmo.sthlp: syntax, options, examples, and returned results.
- src/tmo.jar: optional compiled Java kernel.
- src/tmo.pkg and src/stata.toc: Stata installation files.
- src/tmo.java and src/build.ps1: kernel source and Windows PowerShell build recipe.
- LICENSE: MIT License.
- scpc_tmo/: the SCPC companion, with its own MIT License in scpc_tmo/LICENSE.

The source build is optional for users. See src/README.md for build and test instructions.

## License

TMO is released under the MIT License; see LICENSE.
