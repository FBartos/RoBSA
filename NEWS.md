## development version
### fixes
- `parallel = TRUE` on machines with fewer cores than chains
- spurious BayesTools version warning with development versions of BayesTools

## version 1.1.0
### changes
- update to BayesTools 0.4.0 (seeded results differ from earlier versions and models fitted with earlier versions need to be refitted)

### fixes
- installation on macOS when the pkg-config of the R toolchain does not search `/usr/local/lib/pkgconfig`

## version 1.0.4
### changes
- update to BayesTools 0.3.0
- add `restarts`, `max_extend`, and `check_indicators` to `set_autofit_control()`
- add `check_indicators` to `set_convergence_checks()`

### fixes
- fix `plot_models()` compatibility with RoBSA fits saved under older BayesTools versions

## version 1.0.3
### fixes
- fixing CRAN compilation issues

## version 1.0.2
### features
- adding meandif contrasts to the `prior_factor()` function 

### changes
- update to BayesTools 0.2.14

## version 1.0.1
### changes
- update to R 4.2 and JAGS 4.3.1
- suppressing start-up message 
- cleaning up imports

## version 1.0.0
This is a release version of the package. Please note that it is incompatible to the previous osf version and it is not possible to exactly reproduce previous output.
