using JLD2
using Dates
using Optim
using NLsolve
using Jolab

# We do not recommend to rerun this function without prior study. This script took 7 days to run.

function itfpoint(R, hcav, n, ω, θ, λ)
#This functions calculates an ITF point assuming the mirror reflectivities R, cavity thickness hcav, refractive index n, Gaussian spot size 2ω0 ω, FP mirrors tilted by θ and the ITF point is for the wavelength λ

	# Create the Fabry Perot etalon by setting the position of the mirrors
    ref1 = ReferenceFrame(0., 0, 0, 0, 0)
    ref2 = ReferenceFrame(0, 0, hcav, θ, 0)
    mirror1 = Mirror(R, n, n, ref1)
    mirror2 = Mirror(R, n, n, ref2)
    fp = [mirror1, mirror2]

	#Defines the maximum range of the spatial directions to consider
    srmax = √(4 * 16 / (ω^2 * (2π / λ)^2))
    thresold = 1E-4

	# If very high mirror reflectivity is used the thresold for the ray trace program must be higher
    (R > .99) && (thresold = 1E-5)

	# if the a small spot size, the sampling points of sx and sy must be higher
    length = 100
    (ω < 50E-6) && (length = 150)
    sx = range(-srmax, srmax, length = length)
    sy = range(-srmax, srmax, length = length)

	# Calculates the incident field
    field = FieldAngularSpectrum_gaussian(sx, sy, ω, λ, n, 1, ref1)
    if θ != 0
		# If the mirrors are not parallel the field propagation must be evaluated iteratively
        (fiedlr, fieldt) = lightinteraction_recursive(fp, field, thresold = thresold)
    else
		# If the FP mirrors are parallel the ITF can be calculated with interatively evaluating the algorithm
        (fiedlr, fieldt) = lightinteraction(fp, field)
    end
	#return the intensity of the transmitted field
    return intensity(fieldt)
end

function good_itf(R, hcav, n, ω, θ, λi)
	# Return true or false if the ITF is good, i.e, visibility higher then 0,05
    N = round(2 * n * hcav / λi)
    λaux = (4π * n * hcav) / (N * 2π - π)
    λaux = λaux .+ range(-1E-9, 2E-9, length = 200)
    itf = map((λi) -> itfpoint(R, hcav, n, ω, θ, λi), λaux)
    (maximum(itf) > 0.05) ? true : false
end

function maximize_sens(R, hcav, n, ω, θ, λi)
	# return the wavelength that maximize sensitivity and respective sensitivity
	eps = 1E-12
    i_point(λ) = itfpoint(R, hcav, n, ω, θ, λ)
	# derivative of the itf, i.e the sensitivity
    sens(λ) = (i_point(λ[1] + eps) - i_point(λ[1] - eps)) / 2eps
	# derivative of the sensitivity
    diff_sens!(G, λ) = G[1] = (sens(λ[1] + eps) - sens(λ[1] - eps)) / 2eps
	# second derivative of the sensitivity
    diff2_sens!(G, λ) = G[1] = (sens(λ[1] + eps) - 2sens(λ[1]) + sens(λ[1] - eps)) / eps^2
	# use of Optim.jl to find the wavelength that maximizes the sensitivity
    aux = optimize(sens, diff_sens!, diff2_sens!, [λi], x_tol = 1E-12)
    return (aux.minimizer[1], aux.minimum)
end

function f_θmax(R, hcav, n, ω, λi, θi)
# Function that calculates the α_max

	# Check if the ITF is good, i.e. visibility > 0.05
    bool_good = good_itf(R, hcav, n, ω, 0., λi)
    bool_good || return (-1., false)

	# Calculates the resonance wavelength
    N = round(2 * n * hcav / λi)
    λaux = (4π * n * hcav) / (N * 2π - π)
    ref = minimize_sens(R, hcav, n, ω, 0., λaux)

	# itf point
    i_point(θ,λ) = itfpoint(R, hcav, n, ω,(θ[1]), λ)
    eps = 1E-12
	# sensitivity
    sens(θ) = (i_point(θ[1], ref[1] + eps) - i_point(θ[1], ref[1] - eps)) / 2eps
	# value of the solution
    f!(F, θ) =
             F[1] = sens(θ[1]) - ref[2] * .95
    end
    epsθ = 1E-6
	# derivative of the sensitivity
    g!(G, θ) = G[1] = (sens(θ[1] + epsθ) - sens(θ[1] - epsθ)) / 2epsθ
    initial_θ = [θi];
	# Utilization of NLSolve.jl to find α_max
    aux = nlsolve(f!, g!, initial_θ, xtol = 1E-6, ftol = abs(ref[2] / 1E2), iterations = 10)
    return (aux.zero[1], aux.f_converged)
end

# Range of parameters to simulate

#range of mirror reflectivities
R = range(.9, .999, length = 50)

# range of cavity thicnesses
hcav = range(10E-6, 300E-6, length = 50)

# range of gaussian spot size (ω is 2ω0 in the paper)
ω = range(30E-6, 200E-6, length = 50)

θmax = zeros(length(R), length(ω), length(hcav))
bool = zeros(Bool, length(R), length(ω), length(hcav))

# Refractive index of the cavity
n = 1.
# Wavelength were the interferometric fringe is evaluated
λi = 1550E-9

hcav = repeat(hcav', size(θmax, 2), 1)
ω = repeat(ω, 1, size(θmax, 3))

θs = zeros(length(hcav))
θs .= 0.00005

for iR in eachindex(R)
	# calculates the maximum angle as a function of R, h and 2ω0
    aux = map((hi, ωi, θi) -> f_θmax(R[iR], hi, n, ωi, λi, θi), vec(hcav), vec(ω), vec(θs))

	# Sets the initial condition to the previous interation to be faster to find the solution
    θmax[iR,:,:] = [abs(aux[i][1]) for i in eachindex(aux)]
    bool[iR,:,:] = [abs(aux[i][2]) for i in eachindex(aux)]
    θs .= vec(θmax[iR,:,:])
    @save "alphamax_function.jld2" iR, θmax, bool, R, hcav, ω
end
