using Jolab, Plots, Interpolations, StatsBase
plotly()

# This script uses Jolab to compute the data of fig.1 of [1].

function data(λ, mfd, f_col, f_obj, R1, R2)
	# Function computes the ITF using the angular Airy function
	# λ - vector with the wavelength range
	# mfd - mode field diameter of the single mode fibre
	# f_col - focal length of the collimator
	# f_obj - focal length of the objective
	# R1 - mirror reflectivity of the first FP etalon mirror
	# R2 - mirror reflectivity of the second FP etalon mirror

	# Initialization of refractive indexes to model the FP etalon
	n_fusedsilica = Jolab.refractiveindex_fusedsilica(printBool = false)
	n_air = Jolab.refractiveindex_air(printBool = false)

	## Defining the optical system in Jolab. The optical system and FP etalon
	# is based on [2]
	fibre = SingleModeFibre(mfd, n_air, 1, ReferenceFrame(0,0,0.))
	collimator = Lens(f_col, 1, ReferenceFrame(0,0,f_col))
	objective = Lens(f_obj, 1, ReferenceFrame(0,0,2f_col + f_obj))

	# Sampling of the plane wave directions of propagation
	# sx = sin(θ)cos(ϕ)
	# sy = sin(θ)sin(ϕ)
	nsx = range(-0.2, 0.2, length = 200)
	nsy = range(-0.2, 0.2, length = 200)

	# Direction of the field in Jolab
	dir = 1

	# Cavity thickness
	h = 102E-6

	# Creation of the FP etalon as a multilayer structure
	ref_1 = ReferenceFrame(0.,0.,2f_obj + 2f_col)

	# Creating the FP etalon based on mirrors with reflectivity R
	ref_2 = ReferenceFrame(0.,0.,h)
	mirror1_airy = Mirror(R1, n_air, n_fusedsilica, ref_1)
	mirror2_airy = Mirror(R2, n_fusedsilica, n_air, ref_1 + ref_2)
	fp_airy = [mirror1_airy, mirror2_airy]

	# Initialization of the arrays to store the results
	itf_r_airy = zeros(length(λ))
	itf_t_airy = zeros(length(λ))

	Threads.@threads for i in eachindex(λ)
		## Propagation of the fields by the optical system in Jolab
		field = Jolab.FieldAngularSpectrumScalar_fromfibre(fibre, nsx, nsy, λ[i])
    		(~, field) = lightinteraction(collimator, field)
    		(~, field) = lightinteraction(objective, field)

    		(rfield, tfield) = lightinteraction(fp_airy, field)
    		itf_t_airy[i] = intensity(tfield)

    		(rfield, ~) = lightinteraction(objective, rfield)
    		(rfield, ~) = lightinteraction(collimator, rfield)
    		itf_r_airy[i] = signal(fibre, rfield)
	end

	return (itf_r_airy, itf_t_airy)
end

# Number of wavelengths per ITF
sizeλ = 100

# Initialization of the arrays to store results
itf_r_airy = zeros(4, sizeλ)
itf_t_airy = zeros(4, sizeλ)
λ = zeros(4, sizeλ)

## First ITF

# Defines the wavelength range
λ[1,:] = range(1520.25E-9, 1521.75E-9, length = sizeλ)

# mirror reflectivity values given in fig 6c of [2]
R1_val = 0.9941
R2_val = 0.9933

# The focal length are given by table 1 of [2]
# The values of the mode field diameter were calculated such that the
# field in the back focal plane matches the results measured in [2]
(itf_r_airy[1,:], itf_t_airy[1,:]) = data(λ[1,:], 10.4E-6, 12.56E-3, 300E-3, R1_val, R2_val)

plot(λ[1,:], itf_r_airy[1,:])
plot!(λ[1,:], itf_t_airy[1,:])

## Second ITF

λ[2,:] = range(1543.5E-9, 1545E-9, length = sizeλ)
R1_val = 0.9927
R2_val = 0.9916
(itf_r_airy[2,:], itf_t_airy[2,:]) = data(λ[2,:], 10.4E-6, 18.75E-3, 125E-3, R1_val, R2_val)

plot(λ[2,:], itf_r_airy[2,:])
plot!(λ[2,:], itf_t_airy[2,:])


## Third ITF

λ[3,:] = range(1568.25E-9, 1569.75E-9, length = sizeλ)
R1_val = 0.9909
R2_val = 0.9891
(itf_r_airy[3,:], itf_t_airy[3,:]) = data(λ[3,:], 11.5E-6, 18.75E-3, 75E-3, R1_val, R2_val)

plot(λ[3,:], itf_r_airy[3,:])
plot!(λ[3,:], itf_t_airy[3,:])

## Forth ITF

λ[4,:] = range(1593.0E-9, 1594.5E-9, length = sizeλ)
R1_val = 0.9879
R2_val = 0.9856
(itf_r_airy[4,:], itf_t_airy[4,:]) = data(λ[4,:], 12E-6, 18.75E-3, 50E-3, R1_val, R2_val)

plot(λ[4,:], itf_r_airy[4,:])
plot!(λ[4,:], itf_t_airy[4,:])

# [1] - Marques, Dylan M., et al. "Angular Airy function: a model of Fabry-Perot etalons illuminated by arbitrary beams"
# [2] - Marques, Dylan M., et al. "Modelling fabry-pérot etalons illuminated by focussed beams."
# Optics express 28.5 (2020): 7691-7706. https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-28-5-7691&id=427957
