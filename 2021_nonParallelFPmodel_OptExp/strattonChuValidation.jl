using Jolab;
using Plots; plotly()

# Range of the spatial frequencies for numerical evaluation
sx_X = 2range(-.15, stop = .15, length = 100);

# Creates the non-parallel multilayer structure tilted by 20 mrad
ref_1 = ReferenceFrame(0., 0., 0., 0., 0.)
pointdet = PointDetector(ref_1);
mls1 = MultilayerStructure([1.; 2.; 3.; 1.], [100E-9; 100E-9], ref_1)

ref_2 = ReferenceFrame(0., 0., 20E-6 + 200E-9, 0.02, 0.)
mls2 = MultilayerStructure([1.; 2.; 3.; 1.], [100E-9; 100E-9], ref_2)

mls = [mls1, mls2]

#wavelength range
λ_M = collect(range(1520E-9, stop = 1570E-9, length = 50));

#Calculates incident field
angspe = FieldAngularSpectrum_gaussian.(Ref(sx_X), Ref(sx_X), Ref(28E-6 / sqrt(2)), λ_M, Ref(1), Ref(1), Ref(ref_1));
#Jolab.scalartovectorial!.(angspe);

# Calculates transmitted and reflected field
rtangspe = lightinteraction_recursive.(Ref(mls), angspe);
rangspe = first.(rtangspe);
tangspe = last.(rtangspe);

# Calculates the reflected intensity and the point detection
i2_M = signal.(Ref(pointdet), rangspe);

# Reference signal for normalization
iaM = signal.(Ref(pointdet), angspe)

plot(λ_M, i2_M ./ iaM)
