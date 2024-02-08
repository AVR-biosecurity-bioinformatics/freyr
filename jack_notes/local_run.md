# Working notes on how to run the R pipeline on local computer

Install `R/4.2.0`.

Install `renv/1.0.3`.

Create the `piperline_test` directory on laptop using RStudio to clone Github project, and run next steps in terminal in RStudio. 

Load `renv` R package:

    > library(renv)

Restore packages from the `./renv.lock` file using `renv`:

    > renv::restore(rebuild = T)
    # type 'y' to proceed; packages will install from download or cache (takes a while)

> **NOTE:** If you don't use the `rebuild = T` command (at least on initial run?), only one package will be installed for some reason, and you need to delete the project directory and re-clone a fresh project from Github.

> **NOTE:** Needed to install [Rtools42](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html), as `make` was not found and is needed to install packages from source.


