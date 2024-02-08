# Working notes on how to run the R pipeline in BASC

This will run in the `dros_test` directory in `/group/pathogens/IAWS/Personal/JackS/piperline_tests`. 

    # Change into the main directory you wish to make the project in
    cd /group/pathogens/IAWS/Personal/JackS/piperline_tests

    # Clone the repository
    git clone https://github.com/alexpiper/piperline.git dros_test
    cd dros_test

Load interactive job, modules and R session:

    # Create new interactive SLURM session
    sinteractive --ntasks=1 --cpus-per-task=10 --mem-per-cpu=10GB --time=72:00:00

    module load R/4.2.0-foss-2021b
    module load pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2
    module load GDAL/3.3.0-foss-2021a
    module load BLAST+/2.11.0-gompi-2020a
    module load Pandoc/2.5
    module load ZeroMQ/4.3.2-GCCcore-9.3.0
    module load CMake/3.26.3-GCCcore-12.3.0 # trying this to get 'nanonext' package to install correctly through renv
    # 'CMake/3.10.3-intel-2019a' didn't work
    # ' CMake 3.13 or higher is required.  You are running version 3.10.3 '

    # Load R
    R

Install `renv` package (version `1.0.3`) in my R/4.2.0 library in home directory:

    > install.packages("renv", version = "1.0.3", repos = "http://cran.rstudio.com/", lib = "~/R_libs/4.2.0")

Load `renv` R package:

    > library(renv, lib.loc = "~/R_libs/4.2.0")

Restore packages from `renv.lock` file using `renv`:

    > renv::restore()
    # type 'y' to proceed; packages will install from download or cache (takes a while)

Currently running into this error: 

    - Installing nanonext ...                       FAILED
    Error: Error installing package 'nanonext':
    ====================================

    * installing *source* package ‘nanonext’ ...
    ** package ‘nanonext’ successfully unpacked and MD5 sums checked
    ** using staged installation
    No existing 'libmbedtls' >= 2.5 found
    Detecting 'cmake'...
    /usr/local/eb/software/skylake/CMake/3.26.3-GCCcore-12.3.0/bin/cmake
    Detecting 'xz'...
    /usr/local/eb/software/skylake/XZ/5.4.2-GCCcore-12.3.0/bin/xz
    Compiling 'libmbedtls' from source ...
    cmake: /usr/local/eb/software/skylake/GCCcore/11.2.0/lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found (required by cmake)
    cmake: /usr/local/eb/software/skylake/GCCcore/11.2.0/lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found (required by cmake)
    No existing 'libnng' >= 1.5 found
    Detecting 'cmake'...
    /usr/local/eb/software/skylake/CMake/3.26.3-GCCcore-12.3.0/bin/cmake
    Detecting 'xz'...
    /usr/local/eb/software/skylake/XZ/5.4.2-GCCcore-12.3.0/bin/xz
    Compiling 'libnng' from source ...
    cmake: /usr/local/eb/software/skylake/GCCcore/11.2.0/lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found (required by cmake)
    cmake: /usr/local/eb/software/skylake/GCCcore/11.2.0/lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found (required by cmake)
    ** libs
    gcc -I"/usr/local/eb/software/skylake/R/4.2.0-foss-2021b/lib64/R/include" -DNDEBUG   -I/usr/local/eb/software/skylake/PostgreSQL/13.4-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/OpenSSL/1.1/include -I/usr/local/eb/software/skylake/libgit2/1.1.1-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/MPFR/4.1.0-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/GDAL/3.3.2-foss-2021b/include -I/usr/local/eb/software/skylake/nodejs/14.17.6-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/GLPK/5.0-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/ImageMagick/7.1.0-4-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/GSL/2.7-GCC-11.2.0/include -I/usr/local/eb/software/skylake/UDUNITS/2.2.28-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/HDF5/1.12.1-gompi-2021b/include -I/usr/local/eb/software/skylake/ICU/69.1-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/libsndfile/1.0.31-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/FFTW/3.3.10-gompi-2021b/include -I/usr/local/eb/software/skylake/NLopt/2.7.0-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/GMP/6.2.1-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/libxml2/2.9.10-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/cURL/7.78.0-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/Tk/8.6.11-GCCcore-11.2.0/include -I/usr/local/eb/software/generic/Java/11.0.2/include -I/usr/local/eb/software/skylake/LibTIFF/4.3.0-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/libjpeg-turbo/2.0.6-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/libpng/1.6.37-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/PCRE2/10.37-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/SQLite/3.36-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/zlib/1.2.11-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/XZ/5.2.5-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/bzip2/1.0.8-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/ncurses/6.2-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/libreadline/8.1-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/cairo/1.16.0-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/libGLU/9.0.2-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/Mesa/21.1.7-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/X11/20210802-GCCcore-11.2.0/include -I/usr/local/eb/software/skylake/FlexiBLAS/3.0.4-GCC-11.2.0/include -I/usr/local/eb/software/skylake/FlexiBLAS/3.0.4-GCC-11.2.0/include/flexiblas  -fvisibility=hidden -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno  -c aio.c -o aio.o
    In file included from aio.c:20:
    nanonext.h:22:10: fatal error: nng/nng.h: No such file or directory
    22 | #include <nng/nng.h>
        |          ^~~~~~~~~~~
    compilation terminated.
    make: *** [aio.o] Error 1
    ERROR: compilation failed for package ‘nanonext’
    * removing ‘/group/pathogens/IAWS/Personal/JackS/piperline_tests/dros_test/renv/staging/2/nanonext’
    install of package 'nanonext' failed [error code 1]
    Traceback (most recent calls last):
    13: renv::restore()
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