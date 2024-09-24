using Jolab;
using Plots; plotly()

# This script calculates the ITFs used to validade the non-parallel multilayer structure model agains a parallel multilayer structure model, i.e. the ITFs of figure 4a

# Range of the spatial frequencies for numerical evaluation
sx_X = range(-.15, stop = .15, length = 100);

# Creates the non-parallel multilayer structure tilted by prad
# Creates the first FP mirror with reference frame perpendicular to the z axis
ref_1 = ReferenceFrame(0., 0., 0., 0., 0.)
mls1 = MultilayerStructure([1.; 2.; 3.; 1.], [100E-9; 100E-9], ref_1)

# Creates the second FP mirror with reference frame inclined by 1 prad to the z axis
ref_2 = ReferenceFrame(0., 0., 20E-6 + 200E-9, 1E-12, 0.)
mls2 = MultilayerStructure([1.; 2.; 3.; 1.], [100E-9; 100E-9], ref_2)

# Creates the FP etalon
mls = [mls1, mls2]

# The ITF data has conpared assuming a point detector
pointdet = PointDetector(ref_1);

# Creates the parallel multilayer structure for reference
mlsfull = MultilayerStructure([1.; 2; 3; 1; 2; 3; 1], [100E-9; 100E-9; 20E-6; 100E-9; 100E-9], ref_1)

# wavelength range
λ_M = range(1520E-9, stop = 1570E-9, length = 50)

# Preallocate the output
i_M = zeros(length(λ_M))
inorm_M = zeros(length(λ_M))
iref_M = zeros(length(λ_M))

for i in eachindex(λ_M)

    # Incident field
    angspe = FieldAngularSpectrum_gaussian(sx_X, sx_X, 28E-6 / sqrt(2), λ_M[i], 1, 1, ref_1);

    # Calculates the reflected and transmitted fields of the non-parallel multilayer structure
    (rangspe, tangspe) = lightinteraction_recursive(mls, angspe, thresold = 1E-4);

    # Calculates the intensity on the point detector
    i_M[i] = signal(pointdet, tangspe);

    # Calculates the transmitted and reflected fields of the parallel multilayer structure
    (rangsperef, tangsperef) = lightinteraction(mlsfull, angspe);

    # Calculates the transmitted field intensity at the point (0,0,0)
    iref_M[i] = signal(pointdet, tangsperef);

    # Reference signal for normalization
    inorm_M[i] = signal(pointdet, angspe)
end

plot(λ_M, i_M ./ inorm_M)
scatter!(λ_M[1:2:end], iref_M[1:2:end] ./ inorm_M[1:2:end])
ylims!((0,1))
