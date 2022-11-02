using Jolab, Interpolations, Random, JLD2, MKL

seed = rand(UInt64)
Random.seed!(seed[1])
include(joinpath((@__DIR__), "roughProfileGenerator.jl"))

# k_correlation(1.3671933077726415e8, 1962.0215556503895, 0.787)
rough_model = k_correlation(0.06*(1E3)^4, 500*1E3, .5)
sx = range(-0.7, .7, length = 128*4 + 1); sx = sx[2:end]

n_glass = Jolab.refractiveindex_quartz(printBool = false)
n_air = Jolab.refractiveindex_air(printBool = false)
n_cryo = Jolab.refractiveindex_cryolite(printBool = false)
n_zns = Jolab.refractiveindex_zincsulfide(printBool = false)

λ_desing = 1500E-9
h_cav = 100E-6

h_cryo = λ_desing / 4 / real(n_cryo(λ_desing))
h_zns = λ_desing / 4 / real(n_zns(λ_desing))

n = [n_air; repeat([n_zns, n_cryo], 5); n_glass; repeat([n_zns, n_cryo], 5); n_air];
h = [0; repeat([h_zns, h_cryo], 5); h_cav; repeat([h_zns, h_cryo], 5)]
fp_mls = MultilayerStructure(n, h[2:end], ReferenceFrame(0,0,0))
h = cumsum(h)

roughs = map(i -> roughProfile(sx ./ 1550, rough_model, rfilter = 150E-6), 1:(length(h)))
roughs_f = map(i -> f(x,y) = roughs[i](x*1E9, y*1e9) *1E-9, 1:length(roughs))

fp = map(i -> RoughMultilayerStructure([n[i], n[i+1]], zeros(0), [roughs_f[i]], ReferenceFrame(0,0, h[i])), 1:length(h))
fp_mls2 = map(i -> MultilayerStructure([n[i], n[i+1]], zeros(0), ReferenceFrame(0,0, h[i])), 1:length(h))

function int(λ, fp)
    field = Jolab.FieldAngularSpectrumScalar_gaussian(sx, sx, 70E-6, λ, n_air(λ), 1, ReferenceFrame(0,0,0))
    (fieldr, fieldt) = lightinteraction_recursivegridded(fp, field, rtol = 1E-3)
    return intensity(fieldt)
end

function int(λ, fp::MultilayerStructure)
    field = Jolab.FieldAngularSpectrumScalar_gaussian(sx, sx, 70E-6, λ, n_air(λ), 1, ReferenceFrame(0,0,0))
    (fieldr, fieldt) = lightinteraction(fp, field)
    return intensity(fieldt)
end

λ = 1555.4E-9:.1e-9:1556e-9
λ = 1554.313E-9 .+ 1E-9 .* tan.(range(-1.35, 1.35, length = 50)) / (tan(1.35) - tan(-1.35))
itf = map(λi -> int(λi, fp), λ)
itf_mls = map(λi -> int(λi, fp_mls), λ)


@save "itfglass_$seed.jld2" itf itf_mls λ seed
