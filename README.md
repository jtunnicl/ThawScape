# ThawScape

## Building and Running ThawScape

ThawScape uses the CMake build system and requires a C++11 compiler.

```
mkdir build && cd build     # create a build directory and switch to it
cmake ..                    # run CMake to configure the build
make                        # build the code
ctest                       # run the tests
```

After the build has completed there should be a `ThawScape` executable in the
build directory.

The `ThawScape.ini` input file must be present in the directory you run the
executable from.
