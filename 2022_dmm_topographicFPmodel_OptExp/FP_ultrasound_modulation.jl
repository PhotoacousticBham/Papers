using Jolab, Interpolations, Random, JLD2, MKL

sx = range(-0.7, .7, length = 128*4 + 1); sx = sx[2:end]

n_air = Jolab.refractiveindex_air(printBool = false)
n_quartz = Jolab.refractiveindex_quartz(printBool = false)
n_cryo = Jolab.refractiveindex_cryolite(printBool = false)
n_zns = Jolab.refractiveindex_zincsulfide(printBool = false)

λ_desing = 1500E-9
h_cav = 100E-6

h_cryo = λ_desing / 4 / real(n_cryo(λ_desing))
h_zns = λ_desing / 4 / real(n_zns(λ_desing))

n = [n_air; repeat([n_zns, n_cryo], 5); n_quartz; repeat([n_zns, n_cryo], 5); n_air];
h = [0; repeat([h_zns, h_cryo], 5); h_cav; repeat([h_zns, h_cryo], 5)]
fp_mls = MultilayerStructure(n, h[2:end], ReferenceFrame(0,0,0))
h = cumsum(h)

fp_mls_vec = map(i -> MultilayerStructure([n[i], n[i+1]], zeros(0), ReferenceFrame(0,0, h[i])), 1:length(h))

function int(λ, fp)
    field = Jolab.FieldAngularSpectrumScalar_gaussian(sx, sx, 70E-6, λ, n_air(λ), 1, ReferenceFrame(0,0,0))
    (fieldr, fieldt) = lightinteraction_recursivegridded(fp, field, rtol = 1E-5)
    return intensity(fieldt)
end

function int(λ, fp::MultilayerStructure)
    field = Jolab.FieldAngularSpectrumScalar_gaussian(sx, sx, 70E-6, λ, n_air(λ), 1, ReferenceFrame(0,0,0))
    (fieldr, fieldt) = lightinteraction(fp, field)
    return intensity(fieldt)
end

γ = 97.2E9 # Pa
P = 100E3 # Pa
c = 5570 # m/s
f = 1E6 # Hz
θ = π / 4 # rad

λ = 1554.3E-9:.005e-9:1554.4e-9
# itf_mls = map(λi -> int(λi, fp_mls), λ)

# (val, arg) = findmax(abs, diff(itf_mls))
# λ_bias = (λ[arg+1] + λ[arg]) / 2

# i_bias = int(λ_bias, fp_mls_vec)
# trunc_error = i_bias - int(λ_bias, fp_mls)

us_phase = range(0, 2π, length = 100)
modulation = zeros(length(us_phase))
m = zeros(length(us_phase))
for i_phase in eachindex(modulation)
    m_d(x,y) = P / γ * h_cav * cos(2π * cos(θ) * f / c * x + us_phase[i_phase]) * exp(-4(x^2 + y^2) / (70E-6)^2)
    norm_d(x,y) = exp(-4(x^2 + y^2) / (70E-6)^2)
    m[i_phase] = mapreduce(m_d, +, range(-140E-6, 140E-6, length = 5000), range(-140E-6, 140E-6, length = 5000)') / mapreduce(norm_d, +, range(-140E-6, 140E-6, length = 5000), range(-140E-6, 140E-6, length = 5000)')
    continue

    sound_modulation_mirror1(x,y) = P / γ * h_cav * cos(2π * cos(θ) * f / c * x + us_phase[i_phase]) - m
    sound_modulation_mirror2(x,y) = -P / γ * h_cav * cos(2π * cos(θ) * f / c * x + us_phase[i_phase]) + m

    h = [0; repeat([h_zns, h_cryo], 5); h_cav - 2m; repeat([h_zns, h_cryo], 5)]
    h = cumsum(h)

    fp_modulated = [map(i -> RoughMultilayerStructure([n[i], n[i+1]], zeros(0), [sound_modulation_mirror1], ReferenceFrame(0,0, h[i])), 1:11); map(i -> RoughMultilayerStructure([n[i], n[i+1]], zeros(0), [sound_modulation_mirror2], ReferenceFrame(0,0, h[i])), 12:length(h))]; 

    modulation[i_phase] = i_bias - int(λ_bias, fp_modulated)
end

plot(λ, itf_mls)
plot(us_phase, modulation .- i_bias)