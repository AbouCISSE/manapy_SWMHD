[![Build Status](https://travis-ci.org/imadki/manapy.svg?branch=manapy_hpc)](https://travis-ci.org/imadki/manapy)


# manapy

# install


run 

```python
pip3 install -r requirements.txt
```

```python
python3 -m pip install --user -e .
```

pyccelize functions

```python
./run_pyccel.sh

```

##**********************************************************************##

### Using manapy 
A suite of 1D and 2D examples for the SWMHD equations is available in the `test` directory, more precisely (Test62, Test63, Test64, Test65, Test66,Test68) : 

And the tests (Test1 and Test67) are in models



### Most important :

Viewing is done in paraview

# Test 6.1 : (Experimental order of convergence)
```sh
$ cd /manapy/models/Test61
$ python Test61.py
```

# Test 6.2 : (One dimensional dam break problem in SWMHD)

```sh
$ cd /manapy/tests/Test62/
$ python Test62.py
```

# Test 6.3 : (Two dimensional explosion problem)

```sh
$ cd /manapy/tests/Test63/
$ python Test63.py
```


# Test 6.5 : (C-property ).

```sh
$ cd /manapy/tests/Test65/
$ python Test65.py
```

# Test 6.6 : (C-property test with magnetic field ).

```sh
$ cd /manapy/tests/Test66/
$ python Test66.py
```

# Test 6.7 : (Perturbation of stationary state).

```sh
$ cd /manapy/models/Test67
$ python Test67.py

Result : 
Case 1 : Folder : kx(2,0)
Case 2 : Folder : ky(0,1)
```


# Test 6.8 : (Orszag–Tang-like turbulence problem).

```sh
$ cd /manapy/tests/Test66/
$ python Test68.py
```

