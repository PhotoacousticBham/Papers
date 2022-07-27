using Statistics, FFTW, DSP, Interpolations

struct k_correlation
    a::Float64
    b::Float64
    c::Float64
end

struct gauss_k_correlation
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
end

function psd_distribution(psd_model::k_correlation, sr)
    return .√(psd_model.a .* (1 .+ psd_model.b.^2 .* sr.^2).^(-psd_model.c - 1/2))
end

function psd_distribution(psd_model::gauss_k_correlation, sr)
    return @. √(psd_model.a * exp(-sr^2 * psd_model.b^2) + psd_model.c * (1 + psd_model.d^2 * sr^2)^(-psd_model.e - 1/2))
end

function roughProfile(sx, psd_model; filter_width = 0.01, sr_centre = .2, filter = false, rfilter = 100E-6)
    sizeX = length(sx)
    sizeY = sizeX
    x = fftshift(FFTW.fftfreq(sizeX, 1 / (sx[2] - sx[1])))
    Δx = x[2] - x[1]
    Δy = x[2] - x[1]
    area = (x[end] - x[1])^2
    ceil(Int, sizeX / 2) == 1 && error("length must be even")

    phase = rand(sizeX, sizeX) * (2π) .- π
    # return phase

    phase[1,1] = 0
    phase[1, Int(sizeX/2)+1] = 0
    phase[Int(sizeX/2)+1, Int(sizeX/2)+1] = 0
    phase[Int(sizeX/2)+1, 1] = 0
    phase[2:end, 2:Int(sizeX/2)] = .- rotl90(phase[2:end, Int(sizeX/2)+2:end],2)
    phase[1, 2:Int(sizeX/2)] = .- rotl90(reshape(phase[1, Int(sizeX/2)+2:end],1,:),2)
    phase[Int(sizeX/2) + 2:end,1] = .- rotl90(reshape(phase[2:Int(sizeX/2),1],:,1),2)
    phase[Int(sizeX/2) + 2:end,Int(sizeX/2)+1] = .- rotl90(reshape(phase[2:Int(sizeX/2),Int(sizeX/2)+1],:,1),2)
    # phase .= 0.

    sr = .√(sx.^2 .+ (sx').^2)
    if filter
        filter_function(sr_aux) = exp(-(sr_aux - sr_centre)^2 / filter_width^2 / 4)
        ampli = √(Δx * Δy / sizeX / sizeY) * psd_distribution(psd_model, sr) .* filter_function.(sr)
    else
        ampli = √area * psd_distribution(psd_model, sr)
    end
    # ampli = √area * psd_distribution(psd_model, sr)

    ampli[1,1] = 0
    ampli[1, Int(sizeX/2)+1] = 0
    ampli[Int(sizeX/2)+1, Int(sizeX/2)+1] = 0
    ampli[Int(sizeX/2)+1, 1] = 0
    ampli[2:end, 2:Int(sizeX/2)] = rotl90(ampli[2:end, Int(sizeX/2)+2:end],2)
    ampli[1, 2:Int(sizeX/2)] = rotl90(reshape(ampli[1, Int(sizeX/2)+2:end],1,:),2)
    ampli[Int(sizeX/2) + 2:end,1] = rotl90(reshape(ampli[2:Int(sizeX/2),1],:,1),2)
    ampli[Int(sizeX/2) + 2:end,Int(sizeX/2)+1] = rotl90(reshape(ampli[2:Int(sizeX/2),Int(sizeX/2)+1],:,1),2)

    rough =  ampli .* exp.(im .* phase)
    # return rough
    rough = real.(ifft(ifftshift(rough))) * sizeX * sizeY * (sx[2] - sx[1])^2
    @show √(sum(rough[:].^2) / length(rough))
    # return rough

    rough .+= -mean(rough)

    function filterPos(x, ω, x0)
        if x < -x0
            return exp(-(x + x0)^2 / ω^2 / 4)
        elseif x > x0
            return exp(-(x - x0)^2 / ω^2 / 4)
        else
            return 1.
        end
    end

    rough .*= filterPos.(.√(x.^2 .+ (x').^2), 10E3, rfilter * 1E9)
    # auxtukey = tukey(sizeX , .5)
    # auxtukey = LinearInterpolation(x, auxtukey, extrapolation_bc = 0.)
    # rough .*= auxtukey.(.√(x.^2 .+ (x').^2))
    # stop
    # rough .+= -mean(vec(rough))
    rough = LinearInterpolation((x,x), rough, extrapolation_bc = 0.)
    return rough
end

function radialAveragedPSD(x,y,z)
    sizeX = length(x)
    sizeY = length(y)
    Δx = x[2] - x[1]
    Δy = y[2] - y[1]
    kx = fftshift(FFTW.fftfreq(sizeX, 1 / (x[2] - x[1])))
    ky = fftshift(FFTW.fftfreq(sizeY, 1 / (y[2] - y[1])))

    ffz = fftshift(fft(z)) .* Δx .* Δy
    # return ffz

    kr = range(0, √(maximum(kx.^2) + maximum(ky.^2)), length = max(sizeX,sizeY))
    # kr = range(0, √(maximum(kx.^2) + maximum(ky.^2)), length = 80)
    sizeR = length(kr)
    fftr = zeros(length(kr)-1)

    for iR in 1:sizeR-1
        num = 0
        for iX in 1:sizeX
            for iY in 1:sizeY
                r = √(kx[iX]^2 + ky[iY]^2)
                if kr[iR] < r < kr[iR+1]
                    fftr[iR] = fftr[iR] * num / (num + 1) + abs2(ffz[iX,iY]) / (num + 1)
                    num += 1
                end
            end
        end
    end
    fftr ./= (x[end] - x[1]) * (y[end] - y[1])
    return (kr[1:end-1], fftr)
end
