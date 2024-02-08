# Working notes on how to run the R pipeline on local computer

Install `R/4.2.0`.

Install `renv/1.0.3`.

Create the `piperline_test` directory on laptop using RStudio to clone Github project, and run next steps in terminal in RStudio. 

Load `renv` R package:

    > library(renv)

Restore packages from the `./renv.lock` file using `renv`:

    > renv::restore(rebuild = T)
    # type 'y' to proceed; packages will install from download or cache (takes a while)
    
    Successfully downloaded 257 packages in 2700 seconds.

> **NOTE:** If you don't use the `rebuild = T` command (at least on initial run?), only one package will be installed for some reason, and you need to delete the project directory and re-clone a fresh project from Github.

> **NOTE:** Needed to install [Rtools42](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html), as `make` was not found and is needed to install packages from source.

Threw this error:

    Error: Error installing package 'GenomicRanges':
    =========================================

    '\\internal.vic.gov.au\DEPI\HomeDirs1\js7t\Documents\code\piperline_test'
    CMD.EXE was started with the above path as the current directory.
    UNC paths are not supported.  Defaulting to Windows directory.
    Error in if (any(diff)) { : missing value where TRUE/FALSE needed
    install of package 'GenomicRanges' failed [error code 1]
    Traceback (most recent calls last):
    13: renv::restore(rebuild = T)
    12: renv_restore_run_actions(project, diff, current, lockfile, rebuild)
    11: renv_install_impl(records)
    10: renv_install_staged(records)
    9: renv_install_default(records)
    8: handler(package, renv_install_package(record))
    7: renv_install_package(record)
    6: withCallingHandlers(renv_install_package_impl(record), error = function(e) writef("FAILED"))
    5: renv_install_package_impl(record)
    4: r_cmd_install(package, path)
    3: r_exec_error(package, output, "install", status)
    2: abort(all)
    1: stop(fallback)

After this, running just `> renv:restore()` seems to be working?

    Successfully downloaded 174 packages in 2300 seconds.
    # This doesn't seem like enough

