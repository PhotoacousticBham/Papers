#Used to calculate figure 2
using Jolab, ThreadTools, Interpolations

# Focal length of the collimator
f_col = 10.9E-3

# spatial frequencies
nsr = range(0, .3, length = 1500)

# Design of the system
fibre = SingleModeFibre(10.4e-6, 1, 1, ReferenceFrame(0,0,0))
col = Lens(f_col, 1, ReferenceFrame(0,0,f_col))

# field initialization
field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(fibre, nsr, 1550E-9)
(tmp, field) = lightinteraction(col, field)
r_R = field.r_R
nsr2 = range(0, 0.10, length = 1000)

# refractive index
zns = Jolab.refractiveindex_zincsulfide()
naalf = Jolab.refractiveindex_cryolite()
sio = Jolab.refractiveindex_fusedsilica()
air = Jolab.refractiveindex_air()

# fp design
λ_des = 1402E-9
hcav = 102E-6
n = cat(air, repeat([zns, naalf], 6), sio, repeat([naalf, zns], 6), air, dims = 1)
h = real.(λ_des ./ 4 ./ map(i -> i(λ_des), n[2:end-1]))
h[13] = hcav
fp = MultilayerStructure(n, h, ReferenceFrame(0,0, 2f_col))

# get the reflectivity using bessel beam illumination
function int_bessel(λ, z)
    fp.ref.z = 2f_col + z
    field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(fibre, nsr, λ)
    (tmp, field) = lightinteraction(col, field)
    (tmp, field) = lightinteraction(axicon, field)

    (rfield, tfield) = lightinteraction(fp, field)
    (rfield, tmp) = lightinteraction(axicon, rfield)
    (rfield, tmp) = lightinteraction(col, rfield)
    return (signal(fibre, rfield), intensity(tfield))
end

# used to find the focal plane of the axicon
function findzmax(f, λ, z)
    i_z = map(i -> f(λ, i), z)
    return z[findmax(i_z)[2]]
end

# spatial frequencies
nsr = range(0, .3, length = 1500)

# Design of the optical systems
f_obj = 35E-3
obj = Lens(f_obj, 1, ReferenceFrame(0,0,2f_col + f_obj))
axicon = AxiconFourier(r_R, r_R, nsr2, nsr2, 2.4π/180, ReferenceFrame(0,0,2f_col))
fp_g = MultilayerStructure(n, h, ReferenceFrame(0,0, 2f_col + 2f_obj))

# Used to find the focal plane of the axicon by maximising reflectivity
z_m = range(0, 30E-2, length = 100)
z_max = findzmax(int_bessel, 1550E-9, z_m)

# radial range
r = range(0, 200E-6, length = 2000)

# fourier transform for angular spectrum to field in space
fourier = FourierTransform(r, r, nsr, nsr)

# calculating the angular spectrum of the fields
λ = 1550E-9
field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(fibre, nsr, λ)
(tmp, field) = lightinteraction(col, field)
(tmp, field) = lightinteraction(axicon, field)
changereferenceframe!(field, ReferenceFrame(0,0,2z_max))
(tmp, space_bessel) = lightinteraction(fourier, field)
# electric field: abs2.(space_bessel.e_SXY)

field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(fibre, nsr, λ)
(tmp, field) = lightinteraction(col, field)
(tmp, field) = lightinteraction(obj, field)
(tmp, space_gauss) = lightinteraction(fourier, field)
# plotted field: abs2.(space_gauss.e_SXY