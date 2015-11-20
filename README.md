# MagneticFieldSources
A Julia module for modeling basic magnetic field sources

Models are provided for 
* ideal point dipoles
* [finite continuous solenoids](https://en.wikipedia.org/wiki/Solenoid#Finite_continuous_solenoid)

The
[FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl)
library is used for source and observer positions, and the
[Quaternions](https://github.com/forio/Quaternions.jl) library is used
for rotation.
