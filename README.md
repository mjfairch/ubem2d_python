# Unsteady Boundary-Element Method in 2-D (UBEM2D)
This is a Python-3 implementation of the unsteady boundary-element method of Basu & Hancock for the flow of an inviscid incompressible fluid past a moving airfoil.  The method is described in the following paper:

B.C. Basu and G.J. Hancock, _The unsteady motion of a two-dimensional aerofoil in incompressible inviscid flow_, Journal of Fluid Mechanics, **87** (1), pp. 159-178, 1978.

## Authors of the Code
Michael J. Fairchild and Clancy W. Rowley, Princeton University.

## License
This code may be freely used and modified, and it is distributed under no particular license.  Instead, users of this code are merely asked to acknowledge the authors in any scientific talk or publication whose results make use of this code.

## Dependencies
This code has a few standard dependencies, namely: 

* `numpy`
* `matplotlib`
* `scipy`
* `ffmpeg` (only needed when saving movies to filesystem)

Each of `numpy`, `matplotlib`, and `scipy` are Python packages, which may be installed in the usual way (e.g. `conda install` or `pip install`).  The `ffmpeg` dependency is a library the user installs on their operating system (e.g. with `brew install ffmpeg`, or by visiting the [ffmpeg homepage](http://ffmpeg.org/)).

## Documentation
LaTeX documentation of the underlying fluid dynamics, as well as details of the Python implementation, may be found in the `docs` subdirectory.

## Installation and Verification
After installing the aforementioned dependencies into a Python-3 environment, the user may wish to append the installation directory to one's `PYTHON_PATH` environment variable.  The following examples and test cases may be used to verify the installation and the validity of the code.

### Pre-defined examples
Several examples of increasing complexity are provided in order to demonstrate the basic usage of the library, including:

* Test steady flow past a cylinder: `python examples/cylinder.py`
* Test steady flow past a NACA airfoil: `python examples/naca.py`
* Test steady flow past several NACA airfoils: `python examples/naca_multiple.py`
* Test unsteady flow past a moving NACA airfoil: `python examples/unsteady.py`

### Validation
Several validation cases are provided:

* Validate the steady flow past a cylinder with arbitrary circulation vs. the analytical solution: `python validation/cylinder/cylinder.py`
* Validate the steady flow past a NACA 2412 cambered airfoil with 120 panels at 8-degrees angle of attack against the linear-strength vortex-panel method digitized from Figure 5.25 on p.164 of the book _Foundations of Aerodynamics_, 5th Ed., by Kuethe and Chow: `python validation/naca_kc/naca_kc.py`
* Validate the lift curve (i.e. lift coefficient vs. angle of attack) for steady flow past symmetric and cambered NACA airfoils against experimental data and thin-airfoil theory: `python validation/naca_avd/naca_avd.py`
* Validate rapid pitching motion of an 8.4%-thick von Mises airfoil vs. the U2DIIF Fortran code of the Naval Postgraduate School: `python validation/u2diif_pitch/pitch.py`
* Validate rapid plunging motion of a NACA 0015 airfoil vs. the U2DIIF Fortran code of the Naval Postgraduate School: `python validation/u2diif_plunge/plunge.py`
* Validate a pitch-up maneuver for an 8.4%-thick von Mises airfoil vs. the U2DIIF Fortran code of the Naval Postgraduate School: `python validation/u2diif_ramp/ramp.py`

The experimental data in the lift-curve validation case were digitized from images taken of pages in the book _Theory of Wing Sections_, by Ira H. Abbott and Albert E. von Doenhoff, Dover Publications, 1959.

The `U2DIIF` Fortran code referenced in several of the above cases is described in the following thesis, which may be downloaded from the Naval Postgraduate School's website:

Ngai-Huat Teng, _THE DEVELOPMENT OF A COMPUTER CODE (U2DIIF) FOR THE NUMERICAL SOLUTION OF UNSTEADY, INVISCID AND INCOMPRESSIBLE FLOW OVER AN AIRFOIL_, Naval Postgraduate School, 1987.

### Unit tests
Several unit tests are available in the `ubem2d/tests` directory, which the user may invoke with Python's `unittest` module in the usual manner.  For example:

* `python -m unittest ubem2d/test/test_cylinder.py`
* `python -m unittest ubem2d/test/test_hess_smith.py`
* `python -m unittest ubem2d/test/test_steady_lift_model.py`

In addition, many modules have a `__main__` method which may be used to demonstrate their basic usage.
