using MagneticFieldSources
using Base.Test

length = .1
radius = .01
p = Point3D(100,0,0)
s = Solenoid(length,radius)
d = Dipole()

eps = abs(0.01*hfield(d,p)[3])       # match within 1%

@test_approx_eq_eps hfield(d,p)[3] hfield(s,p)[3] eps
