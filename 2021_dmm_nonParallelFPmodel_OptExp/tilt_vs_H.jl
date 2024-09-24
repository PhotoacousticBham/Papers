cd("/home/dylan/Documents/TiltPaper/")
using Pkg
Pkg.activate("/home/dylan/Documents/TiltPaper/TiltedJolab/")
using Jolab
using Polynomials
using Plots; plotly()
using JLD2


function itfpoint(R, hcav, n, ω, θ, λ)
#This functions calculates an ITF point

	# Create the Fabry Perot etalon by setting the position of the mirrors
        ref1 = ReferenceFrame(0., 0, 0, 0, 0)
        ref2 = ReferenceFrame(0, 0, hcav, θ, 0)
        mirror1 = Mirror(R, n, n, ref1)
        mirror2 = Mirror(R, n, n, ref2)
        fp = [mirror1, mirror2]

	#Defines the maximum range of the spatial directions to consider
        srmax = √(4 * 16 / (ω^2 * (2π / λ)^2))
        sx = range(-srmax, srmax, length = 100)
        sy = range(-srmax, srmax, length = 100)

	#Gets the incident field
        field = Jolab.FieldAngularSpectrum_gaussian(sx, sy, ω, λ, n, 1, ref1)

	#Calculates the reflected and transmitted field
        if θ > 1E-15
                (fiedlr, fieldt) = lightinteraction_recursive(fp, field, thresold = 1E-4)
        else
                (fiedlr, fieldt) = lightinteraction(fp, field)
        end

	# Calculates the intensity of the transmitted field
        return intensity(fieldt)
end

# Condtion of the simulation
ω = 50E-6
λi = 1500E-9
R = .98
hcav = [25E-6, 50E-6, 100E-6, 200E-6]
n = 1.

itf = zeros(100, 100, length(hcav))
λ = zeros(size(itf))
p = fit([0, .5, 1], [0, .4E-3, 1.5E-3])
θ = p.(range(0, 1, length = size(itf,2)))
θ = repeat(θ', size(itf,1), 1)
for i in eachindex(hcav)
        p2 = fit([-1, -.5, .5, 1], [-1.5e-9, -.25E-9, .25E-9, 1.5E-9], 3)
	#Calculates the resonance wavelength
        N = round(2 * n * hcav[i] / λi)
        λaux = (4π * n * hcav[i]) / (N * 2π - π)
        λaux = λaux .+ p2.(range(-1, 1, length = size(itf,1)))
        λ[:,:,i] = repeat(λaux, 1, size(itf,2))
	Threads.@threads for i2 in 1:length(λ[:,1,i])
		for i3 in 1:length(λ[1,:,i])
			itf[i2,i3,i] = itfpoint(R, hcav[i], n, ω, θ[i2,i3], λ[i2,i3,i])
		end
	end
end
θ = θ[1,:]

#Calculates the visibility
v = reshape(maximum(itf, dims = 1), size(itf)[2:3])

#Calculates the sensitivity
sens = reshape(minimum((itf[2:end,:,:] .- itf[1:end-1,:,:]) ./ (λ[2:end,:,:] .- λ[1:end-1,:,:]), dims = 1), size(v))

@save "tilt_vs_H.jld2" v sens θ itf λ

plot(θ, v, ylims = (0, 1), xlims = (θ[1], θ[end]))
plot!(twinx(), θ, sens, ylims = (minimum(vec(sens)), 0), xlims = (θ[1], θ[end]))
savefig("theta_h.pdf")

plot(λ[:,1,1], itf[:,1,1])
plot!(λ[:,43,1], itf[:,43,1])
plot!(λ[:,66,1], itf[:,66,1])
plot!(λ[:,81,1], itf[:,81,1])
ylims!((0, 1))
xlims!((λ[1,1,1], λ[end,1,1]))
savefig("itf_h25.pdf")

plot(λ[:,1,2], itf[:,1,2])
plot!(λ[:,43,2], itf[:,43,2])
plot!(λ[:,66,2], itf[:,66,2])
plot!(λ[:,81,2], itf[:,81,2])
ylims!((0, 1))
xlims!((λ[1,1,2], λ[end,1,2]))
savefig("itf_h50.pdf")

plot(λ[:,1,3], itf[:,1,3])
plot!(λ[:,43,3], itf[:,43,3])
plot!(λ[:,66,3], itf[:,66,3])
plot!(λ[:,81,3], itf[:,81,3])
ylims!((0, 1))
xlims!((λ[1,1,3], λ[end,1,3]))
savefig("itf_h100.pdf")

plot(λ[:,1,4], itf[:,1,4])
plot!(λ[:,43,4], itf[:,43,4])
plot!(λ[:,66,4], itf[:,66,4])
plot!(λ[:,81,4], itf[:,81,4])
ylims!((0, 1))
xlims!((λ[1,1,4], λ[end,1,4]))
savefig("itf_h200.pdf")
