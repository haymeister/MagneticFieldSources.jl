module MagneticFieldSources

using StaticArrays
const Vector3D = MVector{3}
const Point3D = Vector3D

Vector3D(x::T,y::T,z::T) where T = Vector3D{T}(x,y,z)
Vector3D(x) = Vector3D(x,x,x)
Vector3D(x,y,z) = Vector3D(promote(x,y,z)...)

using Elliptic: E,K,Pi
using Quaternions: rotationmatrix, Quaternion

export Dipole,Solenoid,hfield,Vector3D,Point3D

abstract type MagneticFieldSource end

"""
    Dipole{T<:AbstractFloat} <: MagneticFieldSource

Magnetic dipole with 3D position and moment (orientation)
"""
struct Dipole{T<:AbstractFloat} <: MagneticFieldSource
    position::MVector{3,T}
    moment::MVector{3,T}
end
Dipole() = Dipole(Vector3D(0.),Vector3D(0.,0.,1.))
Dipole(position::Point3D) = Dipole(position,Vector3D(0.,0.,1.))

"""
    hfield(d::Dipole,p::Vector3D)

The magnetic field from an ideal magnetic dipole
"""
function hfield(d::Dipole,p::Vector3D)
    r = p-d.position
    return (3r*dot(d.moment,r)/norm(r)^2 - d.moment)/(4π*norm(r)^3)
end

"""
    Solenoid{T<:AbstractFloat} <: MagneticFieldSource

A finite-length solenoid modeled as a cylindrical current sheet
"""
type Solenoid{T<:AbstractFloat} <: MagneticFieldSource
    position::MVector{3,T}
    moment::MVector{3,T}
    length::T
    radius::T
end

Solenoid(length::Number,radius::Number) = Solenoid(Point3D(0.,0.,0.),Vector3D(0.,0.,1.),length,radius)
Solenoid(position::Point3D,length::Number,radius::Number) = Solenoid(position,Vector3D(0.,0.,1.),length,radius)

"""
    hfield(s::Solenoid,p::Point3D)

The magnetic field from a current sheet. See [here](https://en.wikipedia.org/wiki/Solenoid#Finite_continuous_solenoid).
"""
function hfield(s::Solenoid,p::Point3D)
    NI = norm(s.moment)/(π*s.radius^2) # turn-current product
    r = p-s.position
    # get rotation angles
    φ = atan2(s.moment[1],s.moment[2]) # angle about z axis
    θ = atan2(sqrt(s.moment[1]^2+s.moment[2]^2),s.moment[3]) # angle about x axis
    r = rotationmatrix(quatfromeuler(0,θ,φ))*r
    ρ = max(sqrt(r[1]^2+r[2]^2),eps(Float64))
    z = r[3]

    ζp = z+s.length/2; ζm = z-s.length/2
    Hz = (intz(ρ,ζp,s.radius) - intz(ρ,ζm,s.radius))/(4π*s.length*sqrt(s.radius*ρ))

    Hρ = (sqrt(s.radius/ρ)/(2π*s.length))*(intρ(ρ,ζp,s.radius) - intρ(ρ,ζm,s.radius))

    α=atan2(r[2],r[1])          # angle about z-axis
    H = Vector3D(Hρ*cos(α),Hρ*sin(α),Hz)*NI

    # now rotate field based on the orientation of the moment
    H = rotationmatrix(quatfromeuler(-φ,-θ,0))*H
    return H
end

function intz(ρ,ζ,a)
    h²=4a*ρ/(a+ρ)^2
    k² = 4a*ρ/((a+ρ)^2 + ζ^2)
    return ζ*√k²*(K(k²) + ((a-ρ)/(a+ρ))Pi(h²,π/2,k²))
end

function intρ(ρ,ζ,a)
    k² = (4a*ρ)/((a+ρ)^2 + ζ^2)
    return ((k²-2)K(k²) + 2E(k²))/√k²
end

"""
    quatfromeuler(φ,θ,ψ)
Returns a quaternion based on the z-x-z Euler angles. See
[here](https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.E2.86.94_Quaternion)
"""
function quatfromeuler(φ,θ,ψ)
    Quaternion(cos((φ+ψ)/2)*cos(θ/2),cos((φ-ψ)/2)*sin(θ/2),sin((φ-ψ)/2)*sin(θ/2),sin((φ+ψ)/2)*cos(θ/2))
end

end
