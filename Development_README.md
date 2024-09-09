### Development(Changing source code and rebuilding packages)
If you make change to the source code and want to rebuild the package. In the folder containing the LevSeq repo, build the package again.

```
python setup.py sdist bdist_wheel
pip install dist/levseq-0.1.0.tar.gz
```

### Steps to rebuild the C++ executables

This is to run the code locally, rather than via the docker instance, note this will be dependent on your computer and may not 
function as above, we highly recomend using the docker version for this reason.

### Mac intel chip
To rebuild on mac move into the `source/source` folder and execute the following command:

```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/local/bin/gcc-13 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-13 ../source
```

Note it expects c to be installed in `/usr/local/bin` if this is not the case on your machine you will need to update 
accordingly. 

After building you need to make the make file with the following command:

```
make -j
```

The demultiplex file should now function!
