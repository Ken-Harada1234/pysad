"""
pysad - Python translation of the SAD (Strategic Accelerator Design) code

SAD was originally written by K. Oide (KEK) in SADScript (a Mathematica-like language)
with a C/Fortran computational core. This package provides equivalent functionality
in pure Python with NumPy.

Excluded (no GUI, no device control):
  - Tcl/Tk display
  - X-window graphics
  - EPICS accelerator control

Modules
-------
constants   Physical and accelerator constants (from SAD's init.n)
elements    Accelerator element definitions and linear transfer matrices
optics      Beam optics: symplectic algebra, Twiss, beam matrix, coupling
lattice     Lattice: Twiss propagation, radiation integrals, emittance
tracking    Exact symplectic particle tracking (matching SAD's Fortran core)
parser      SADScript (.n file) interpreter

Quick start
-----------
>>> from pysad import fodo_cell, track_particles
>>> lat = fodo_cell(L=2.0, theta=0.1, k1=1.2)
>>> tw = lat.twiss()
>>> print(lat.summary())
"""

from .constants import (
    SpeedOfLight, ElectronMass, ElectronRadius,
    ProtonMass, Cgamma, Cq,
)

from .elements import (
    Element, Drift, Quadrupole, Bend, Sextupole,
    Cavity, Monitor, Marker, Multipole,
    make_element,
)

from .optics import (
    symplectic_j, symplectic_inverse,
    normal_coordinates, transfer_matrix,
    beam_matrix, periodic_r, periodic_optics,
    dispersion_from_6d, geo_base,
)

from .lattice import (
    Lattice,
    fodo_cell,
)

from .parser import (
    SADScript,
    load_sad_file,
    load_sad_string,
    load_sad_files,
    SADEnvironment,
)

from .tracking import (
    pxy2dpz,
    track_drift_exact,
    track_thin_quad,
    track_thin_sext,
    track_multipole_kick,
    track_thick_quad,
    track_thick_bend,
    track_cavity,
    radiation_loss_per_turn,
    damping_times,
    track_element,
    track_particles,
    TrackResult,
    dynamic_aperture,
    tune_by_fft,
)

from .matching import (
    Vary,
    Constraint,
    Match,
    MatchResult,
    DynamicApertureSurvey,
)

from .ffs import (
    FFSState,
    SADList,
)

__version__ = "0.1.0"
__author__ = "SAD Python translation (KEK-Acc6th)"
__all__ = [
    # constants
    "SpeedOfLight", "ElectronMass", "ElectronRadius",
    "ProtonMass", "Cgamma", "Cq",
    # elements
    "Element", "Drift", "Quadrupole", "Bend", "Sextupole",
    "Cavity", "Monitor", "Marker", "Multipole", "make_element",
    # optics
    "symplectic_j", "symplectic_inverse",
    "normal_coordinates", "transfer_matrix",
    "beam_matrix", "periodic_r", "periodic_optics",
    "dispersion_from_6d", "geo_base",
    # lattice
    "Lattice", "fodo_cell",
    # tracking
    "pxy2dpz", "track_drift_exact",
    "track_thin_quad", "track_thin_sext", "track_multipole_kick",
    "track_thick_quad", "track_thick_bend",
    "track_element", "track_particles", "TrackResult",
    "dynamic_aperture", "tune_by_fft",
    # parser
    "SADScript", "load_sad_file", "load_sad_string", "load_sad_files", "SADEnvironment",
    # matching
    "Vary", "Constraint", "Match", "MatchResult",
    "DynamicApertureSurvey",
    # ffs
    "FFSState", "SADList",
]
