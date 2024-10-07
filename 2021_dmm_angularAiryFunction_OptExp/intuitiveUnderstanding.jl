using Jolab, Plots

# Initialization of both mirrors that form the FP cavity
mirror1 = Mirror(.97,1,1,ReferenceFrame(0,0,0.))
mirror2 = Mirror(.97,1,1,ReferenceFrame(0,0,10E-6))

# Initialization of the FP etalon based on the mirrors
fp = [mirror1, mirror2]

## Create data from fig. 3b
# Sampling of the plane wave directions of propagation
# sx = sin(θ)cos(ϕ)
# sy = sin(θ)sin(ϕ)
sx = range(-.3, .3, length = 512)
sy = range(-.3, .3, length = 512)

# wavelength
λ = 1596E-9

# Initialization of the incident field up the FP etalon assumed to be Gaussian
field = FieldAngularSpectrumScalar_gaussian(sx,sy,10E-6,λ,1,1,ReferenceFrame(0,0,0.))


# Calculation of the reflected and transmitted angular spectra
(fieldr,fieldt) = lightinteraction(fp, field)

# Data shown in fig. 3b
# The first dimension is only valid for vectorial simulations in Jolab
# Absolute value of the incident, reflected and transmitted plane waves components.
idata = abs.(reshape(field.e_SXY, 512, 512))
tdata = abs.(reshape(fieldt.e_SXY, 512, 512) ./ maximum(idata)
rdata = abs.(reshape(fieldr.e_SXY, 512, 512) ./ maximum(idata)
idata = idata ./ maximum(idata)

## Creating data for fig. 3a

# Sampling of the plane wave directions of propagation
# sx = sin(θ)cos(ϕ)
# sy = sin(θ)sin(ϕ)
sx = range(-.15, .15, length = 200)
sy = range(-.15, .15, length = 200)

# Wavelength range considered
λ = range(1570E-9, 1610E-9, length = 512)

# Initialization of the itf array
itf = zeros(length(λ))

# Calculating the transmitted intensity by creating a for cycle changing the light wavelength
for i in eachindex(λ)
    # Initialization of the incident field
    field = FieldAngularSpectrumScalar_gaussian(sx,sx,10E-6,λ[i],1,1,ReferenceFrame(0,0,0.))

    # Calculation of the reflected and transmitted fields
    (fieldr, fieldt) = lightinteraction(fp, field)

    # Calculation of the transmitted intensity
    itf[i] = intensity(fieldt)
end

# Plotting the ITF of the inset of fig.3a
plot!(λ, itf)
cosθ = λ .* (13 .- 1/2) ./ 2 ./ 1 ./ 10E-6

# Increase the sampling for better plots sampling
sx = range(-.3, .3, length = 512)

# Choose the wavlength to plot the angular spectrum
i = 321

# Initialization of the initial field
field = FieldAngularSpectrumScalar_gaussian(sx,[0.],10E-6,λ[i],1,1,ReferenceFrame(0,0,0.))

# Calculation of the transmitted and reflected field
(fieldr, fieldt) = lightinteraction(fp, field)

# Plotting the angular spectrum of fig 3.a
plot(field.nsx_X, abs.(field.e_SXY[1,:,1]))

# To create the data from fig.4, these script only needs to be edited to change the mirror reflectivity
