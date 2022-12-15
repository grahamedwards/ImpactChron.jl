## Recalibrate ages based on decay constants
using ImpactChron: agerecal
using DelimitedFiles

cd(@__DIR__)

data_mat = readdlm("KArArages_compilation.csv",',')
data_headers = Symbol.(data_in_mat[1,:])
d = Dict{Symbol,Any}()
for i ∈ eachindex(data_headers)
    k = data_headers[i]
    dv = data_mat[2:end,i]
    k ∈ [:ageMa, :sigma] ? d[k] = float.(dv) : d[k] = string.(dv)
end

# Create copies of ages and uncertainties to overwrite with recalibrated data
ageMa = copy(d[:ageMa])
σ = copy(d[:sigma])


# Identify indices of Ar-Ar ages and K-Ar ages
indices_ArAr = .!contains.(d[:notes],"K-Ar")
indices_KAr = contains.(d[:notes],"K-Ar")


## Correct K-Ar ages for new decay constants. 
iKArλ = findall( @. indices_KAr * (d[:lambda] == "T69" || d[:lambda] == "H74"))
for i ∈ iKArλ
    ageMa[i],σ[i] = agerecal(d[:ageMa][i],d[:sigma][i],KAr=true)
end


## Papers using the Husain 1974 decay constants and monitors.
iH74 = findall(indices_ArAr .* (d[:lambda] .== "H74"))
for i ∈ iH74 
    ageMa[i],σ[i] = agerecal(d[:ageMa][i],d[:sigma][i],monitor_age=2.668e9)
end

## Papers using the Turner 1969 decay constants, unknown monitor.
iT69 = findall(indices_ArAr .* (d[:lambda] .== "T69"))
for i ∈ iT69
    ageMa[i],σ[i] = agerecal(d[:ageMa][i],d[:sigma][i])
end

outmat = hcat(data.name, data.group, data.type, ageMa,σ)

writedlm("KArArages.csv", outmat,',')
writedlm("ArArages.csv", outmat[indices_ArAr,:],',')
writedlm("KArages.csv", outmat[indices_KAr,:],',')