using Jolab, ThreadTools

f_col = 10.9E-3

nsr = range(0, .3, length = 1500)

fibre = SingleModeFibre(10.4e-6, 1, 1, ReferenceFrame(0,0,0))
col = Lens(f_col, 1, ReferenceFrame(0,0,f_col))

field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(fibre, nsr, 1550E-9)
(tmp, field) = lightinteraction(col, field)
r_R = field.r_R
nsr2 = range(0, 0.10, length = 1000)

zns = Jolab.refractiveindex_zincsulfide()
naalf = Jolab.refractiveindex_cryolite()
sio = Jolab.refractiveindex_fusedsilica()
air = Jolab.refractiveindex_air()

λ_des = 1402E-9
hcav = 102E-6
n = cat(air, repeat([zns, naalf], 6), sio, repeat([naalf, zns], 6), air, dims = 1)
h = real.(λ_des ./ 4 ./ map(i -> i(λ_des), n[2:end-1]))
h[13] = hcav

fp = MultilayerStructure(n, h, ReferenceFrame(0,0, 2f_col))

function int_gauss(λ)
    field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(fibre, nsr, λ)
    (tmp, field) = lightinteraction(col, field)
    (tmp, field) = lightinteraction(obj, field)
    (rfield, tfield) = lightinteraction(fp_g, field) 
    (rfield, tmp) = lightinteraction(obj, rfield)
    (rfield, tmp) = lightinteraction(col, rfield)
    return (signal(fibre, rfield), intensity(tfield))
end

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

function findzmax(f, λ, z)
    i_z = map(i -> f(λ, i), z)
    return z[findmax(i_z)[2]]
end

# λ_m = range(1500E-9, 1630E-9, length = 10000)
λ_m = [i .+ range(-1.5E-9, .5E-9, length = 1000) for i in (1502.4, 1509.8, 1517.5, 1525.1, 1532.9, 1540.7, 1548.6, 1556.6, 1564.6, 1572.7, 1581.0, 1589.3, 1597.6, 1606.2, 1614.8, 1623.4) .* 1E-9]

nsr = range(0, .3, length = 1500)

f_obj = 150E-3 # must be 35E-3 for the focussed Gaussian beam
obj = Lens(f_obj, 1, ReferenceFrame(0,0,2f_col + f_obj))


axicon = AxiconFourier(r_R, r_R, nsr2, nsr2, 0.4π/180, ReferenceFrame(0,0,2f_col)) 
# for the focused Bessel beam:
# axicon = AxiconFourier(r_R, r_R, nsr2, nsr2, 2.4π/180, ReferenceFrame(0,0,2f_col))
fp_g = MultilayerStructure(n, h, ReferenceFrame(0,0, 2f_col + 2f_obj))

z_m = range(0, 30E-2, length = 100)
# for the focused Bessel beam:
# z_m = range(0, 3E-2, length = 100)

itf_bessel = map(i -> (data = int_bessel.(i, findzmax(int_bessel, i[1], z_m)); (first.(data), last.(data))), λ_m)
stop
itf_gauss = tmap(i -> (data = int_gauss.(i); (first.(data), last.(data))), λ_m)

t_bessel_fringes_m = map(i -> i[2], itf_bessel)
r_bessel_fringes_m = map(i -> i[1], itf_bessel)

t_gauss_fringes_m = map(i -> i[2], itf_gauss)
r_gauss_fringes_m = map(i -> i[1], itf_gauss)