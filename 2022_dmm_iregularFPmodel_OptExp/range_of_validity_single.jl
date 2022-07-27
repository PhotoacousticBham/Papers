using Jolab, Plots, JLD2, MKL, Interpolations, ThreadTools
plotly()


# Vectors of the cavity offset
layer_offset = (-15:.25:15) .* 1E-9

function int(λ, fp)
    field = Jolab.FieldAngularSpectrumScalar_gaussian(sx, sx, 70E-6, λ, n_air(λ), 1, ReferenceFrame(0,0,0))
    (field_r, field_t) = lightinteraction_recursivegridded(fp, field, rtol = 1E-3)
    return intensity(field_t)
end

function int(λ, fp::MultilayerStructure)
    field = Jolab.FieldAngularSpectrumScalar_gaussian(sx, sx, 70E-6, λ, n_air(λ), 1, ReferenceFrame(0,0,0))
    (field_r, field_t) = lightinteraction(fp, field)
    return intensity(field_t)
end

# Sampling of kx and ky. kx = k sx and ky = k sy 
sx = range(-.2, .2, length = 128)

# Wavelenth points to calculate ITF
λ = range(1553E-9, 1555.5E-9, length = 50)

# Refractive index of all layers forming the mirror and cavity
n_air = Jolab.refractiveindex_air(printBool = false)
n_glass = Jolab.refractiveindex_quartz(printBool = false)
n_cryo = Jolab.refractiveindex_cryolite(printBool = false)
n_zns = Jolab.refractiveindex_zincsulfide(printBool = false)

λ_desing = 1500E-9
h_cav = 100E-6

h_cryo = λ_desing / 4 / real(n_cryo(λ_desing))
h_zns = λ_desing / 4 / real(n_zns(λ_desing))

n = [n_air; repeat([n_zns, n_cryo], 5); n_glass; repeat([n_zns, n_cryo], 5); n_air];
h = [0; repeat([h_zns, h_cryo], 5); h_cav; repeat([h_zns, h_cryo], 5)]
fp_ref = MultilayerStructure(n, h[2:end], ReferenceFrame(0,0,0))

mirror1 = MultilayerStructure(n[1:12], h[2:11], ReferenceFrame(0,0,0))
h_cum = cumsum(h)

itf = zeros(length(λ), length(layer_offset))
itf_ref = zeros(length(λ), length(layer_offset))
for i in eachindex(layer_offset)
    non_flat_profile(x,y) = layer_offset[i]
	h_cav = 100E-6 + layer_offset[i]
	h = [0; repeat([h_zns, h_cryo], 5); h_cav; repeat([h_zns, h_cryo], 5)]
	fp_ref = MultilayerStructure(n, h[2:end], ReferenceFrame(0,0,0))
	h_cum = cumsum(h)

	fp = [mirror1; map(i -> RoughMultilayerStructure([n[i], n[i+1]], zeros(0), [non_flat_profile], ReferenceFrame(0,0, h_cum[i])), 12:length(h_cum))]

    itf[:,i] .= map(λi -> int(λi, fp), λ)
	itf_ref[:,i] .= tmap(λi -> int(λi, fp_ref), λ)
end

# filename = homedir() * "/approx_glass_single.jld2"
# @save filename λ itf itf_ref layer_offset

function measurefwhm(λ, i)
    # Function calculating the FWHM of an ITF
    (maxi, argmax) = findmax(i);
	(mini, argmin) = findmin(i);
	ifwhm = maxi - (maxi - mini) / 2 ;
	argfirst = findfirst(i[1:argmax] .> ifwhm);
	argfirst === nothing && return 0

	arglast = findlast(i[argmax:end] .> ifwhm);
	arglast === nothing && return 0
	arglast = arglast + argmax-1;

	itp = LinearInterpolation(([i[argfirst-1]; i[argfirst]], ), [λ[argfirst-1]; λ[argfirst]])
	λfirst = itp(ifwhm);

	itp = LinearInterpolation(([i[arglast+1]; i[arglast]], ), [λ[arglast+1]; λ[arglast]])
	λlast = itp(ifwhm);
	return λlast - λfirst;
end

measurevisibility(i) = maximum(i)

# Function calculating the sensitivity
measuresensitivity(λ, i) = maximum(abs.((i[2:end] .- i[1:end-1]) ./ (λ[2:end] .- λ[1:end-1])))

# Calculation of the diferent metrics for the different ITFs
fwhmref = [measurefwhm(λ, itf_ref[:,i]) for i in 1:size(itf_ref,2)]
fwhm = [measurefwhm(λ, itf[:,i]) for i in 1:size(itf,2)]
vref = [measurevisibility(itf_ref[:,i]) for i in 1:size(itf_ref,2)]
v = [measurevisibility(itf[:,i]) for i in 1:size(itf,2)]
sensref = [measuresensitivity(λ, itf_ref[:,i]) for i in 1:size(itf_ref,2)]
sens = [measuresensitivity(λ, itf[:,i]) for i in 1:size(itf,2)]

# Plot figure 4c
fig4c = plot()
plot!(fig4c, λ, itf[:, 1:8:end], c = :blue)
plot!(fig4c, λ, itfref[:, 1:8:end], c = :red)

# Plot figure 4d
fig4d = plot()
# plot!(fig4d, layer_offset, abs.((sensref .- sens)) ./ sensref)
plot!(fig4d, layer_offset, abs.((fwhmref .- fwhm)))
plot!(fig4d, layer_offset, abs.((vref .- v)))
vline!(fig4d, [-8.5E-9; 8.5E-9])