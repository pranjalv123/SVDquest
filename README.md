# SVDquest

Building accurate SVDquartets-based species trees

To build SVDquest you will first need [phylokit](https://github.com/pranjalv123/phylokit) and [phylonaut](https://github.com/pranjalv123/phylonaut),

    mkdir build
    cd build
    cmake ../src
    make


Ideally, you will have wine installed, and the Windows version of PAUP*. The Windows version has better numerical routines and returns substantially more accurate trees than the Linux version.

To check if this is installed correctly, if you run `wine paup4c` from your command line, it should start PAUP*.

Alternatively, you can use the Linux version of PAUP*. In this case you should run SVDquest with the `--nowine` option.

If you are running on Windows or using the Linux version of PAUP* on Linux, you can specify a custom PAUP* executable using the --paup-exe option followed by the path to the executable.

You also need ASTRAL in the same folder as SVDquest. So if you are in the folder with the SVDquest executable, Astral/astral.4.7.8.jar (or any other version of ASTRAL) should exist. 


To run SVDquest*, run

    SVDquest -i <input gene trees> -a <input alignment> -o <output file>


SVDquest* uses [SIESTA](https://link.springer.com/chapter/10.1007/978-3-319-67979-2_13) to generate strict, greedy, and majority consensus trees of all optimal trees, as well as a single arbitrarily chosen optimal tree.
If for some reason this is a problem, you can disable these with --nogreedy, --nostrict, and --nomajority options.



You can also use the `--score` option to score a tree like this:

     SVDquest --score -a <input alignment> -o <output file>
