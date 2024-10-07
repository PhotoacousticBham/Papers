using Jolab, Plots

# This script uses Jolab to compare two Fabry-Pérot (FP) etalon models.
# The two models are a multilayer model and the angular Airy function.
# The multilayer model assumes the FP etalon as a stack of dieletric materials
# and requires the thickness and refractive index of each layer to be know.
# The angular Airy function describes the FP etalon based only on its mirror
# reflectivities and optical cavity thickness.

function data(λ, ω)
	# This function computes the ITF for both models.
	# λ - wavelength range
	# ω - Gausian spot size.

	## Initialization of the multilayer structure
	# Initialization of refractive indexes to model the FP etalon
	n_fusedsilica = Jolab.refractiveindex_fusedsilica(printBool = false)
	n_cryo = Jolab.refractiveindex_cryolite(printBool = false)
	n_zns = Jolab.refractiveindex_zincsulfide(printBool = false)
	n_air = Jolab.refractiveindex_air(printBool = false)

	# Sampling of the plane wave directions of propagation
	# sx = sin(θ)cos(ϕ)
	# sy = sin(θ)sin(ϕ)
	nsx = range(-0.15, 0.15, length = 150)
	nsy = range(-0.15, 0.15, length = 150)

	# Direction of the field in Jolab
	dir = 1

	# Cavity thickness
	h = 102E-6

	# Design wavelength of the dielectric mirror
	λ_des = 1402E-9

	# Creation of the FP etalon as a multilayer structure
	ref_1 = ReferenceFrame(0.,0.,0.,0.,0.)
	n_M = cat(n_air, repeat([n_cryo, n_zns], 6), n_fusedsilica, repeat([n_zns, n_cryo], 6), dims = 1)
	h_M = λ_des ./ 4 ./ [n_M[i](λ_des) for i in 2:length(n_M)-1]
	h_M[13] = h
	fp_mls = MultilayerStructure(n_M, h_M, ref_1)

	## Initialization of the FP etalon as the angular Airy function

	# We use the multilayer model to compute the FP mirror reflectivities
	# Assuming a nearly collimate incident beam
	nsx_col = range(-0.05, 0.05, length = 100)
	nsy_col = range(-0.05, 0.05, length = 100)
	field_aux = FieldAngularSpectrumScalar_gaussian(nsx_col, nsy_col, 200E-6, λ[1], n_air(λ[1]), dir, ref_1)

	# Creates a multilayer structure with the layers of the FP etalon that form the mirror
	mirror1 = MultilayerStructure(n_M[1:14], h_M[1:12], ref_1)
	(fieldR,tmp) = lightinteraction(mirror1, field_aux)

	# The mirror reflectivity is the ratio between the intensity of the reflected
	# and incident field
	R = intensity(fieldR) / intensity(field_aux)


	# Creating the FP etalon based on mirrors with reflectivity R
	ref_2 = ReferenceFrame(0.,0.,h)
	mirror1_airy = Mirror(R, n_air, n_fusedsilica, ref_1)
	mirror2_airy = Mirror(R, n_fusedsilica, n_air, ref_2)
	fp_airy = [mirror1_airy, mirror2_airy]

	# Initialization of the arrays to store the results
	itf_t_mls = zeros(length(λ))
	itf_r_mls = zeros(length(λ))
	itf_r_airy = zeros(length(λ))
	itf_t_airy = zeros(length(λ))

	Threads.@threads for i in eachindex(λ)
		## Propagation of the fields by the optical system in Jolab
		field = FieldAngularSpectrumScalar_gaussian(nsx, nsy, ω, λ[i], n_air(λ[i]), dir, ref_1)
    		(rfield, tfield) = lightinteraction(fp_mls, field)
    		itf_t_mls[i] = intensity(tfield)
    		itf_r_mls[i] = intensity(rfield)

    		(rfield, tfield) = lightinteraction(fp_airy, field)
    		itf_r_airy[i] = intensity(rfield)
		itf_t_airy[i] = intensity(tfield)
	end
	return (itf_r_airy, itf_t_airy, itf_r_mls, itf_t_mls)
end

# Number of wavelengths per ITF
sizeλ = 1000

# Initialization of the arrays to store results
itf_t_mls = zeros(4, sizeλ)
itf_r_mls = zeros(4, sizeλ)
itf_r_airy = zeros(4, sizeλ)
itf_t_airy = zeros(4, sizeλ)
λ = zeros(4, sizeλ)

# Compute the ITF on a specific wavelength range and spot size
λ[1,:] = range(1400E-9, 1407E-9, length = sizeλ)
(itf_r_airy[1,:], itf_t_airy[1,:], itf_r_mls[1,:], itf_t_mls[1,:]) = data(λ[1,:], 30E-6)

# Compute the ITF on a specific wavelength range and spot size
λ[2,:] = range(1500E-9, 1507E-9, length = sizeλ)
(itf_r_airy[2,:], itf_t_airy[2,:], itf_r_mls[2,:], itf_t_mls[2,:]) = data(λ[2,:], 50E-6)

# Compute the ITF on a specific wavelength range and spot size
λ[3,:] = range(1549.5E-9, 1556.5E-9, length = sizeλ)
(itf_r_airy[3,:], itf_t_airy[3,:], itf_r_mls[3,:], itf_t_mls[3,:]) = data(λ[3,:], 8.5E-6)

# Compute the ITF on a specific wavelength range and spot size
λ[4,:] = range(1600E-9, 1607E-9, length = sizeλ)
(itf_r_airy[4,:], itf_t_airy[4,:], itf_r_mls[4,:], itf_t_mls[4,:]) = data(λ[4,:], 250E-6)

# Ploting the data for visualization. The offset in the wavelength axis is
# because both models do not agree in the fring position but agree very well in
# the fringe shape
plot(λ[2,:] .- 4.6E-9, itf_r_airy[2,:])
plot!(λ[2,:] .- 4.6E-9, itf_t_airy[2,:])
plot!(λ[2,:], itf_r_mls[2,:])
plot!(λ[2,:], itf_t_mls[2,:])
