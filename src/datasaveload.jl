# Efficiently save and load data with DelimitedFiles functions.

    # nt2csv
    # dict2csv
    # data2csv
    # csv2nt
    # csv2dict
    # serial2dict

"""

```julia
nt2csv(filename::String,N::NamedTuple)
```

Save a `NamedTuple` to a .csv. Accounts for fields of any type.

For single-element entries, the rest of the array/table is filled with `NaN`s.

see also: [`dict2csv`](@ref), [`data2csv`](@ref)

"""

function nt2csv(filename::String,N::NamedTuple)
    kz = keys(N)
    nₖ = length(kz)
    rows = maximum( length(N[i]) for i ∈ kz)
    A = Array{Any}(undef,rows+1,nₖ) # preallocate output array
    for i ∈ 1:nₖ
        k = kz[i]
        if length(N[k]) == rows
            A[:,i] = vcat(k,N[k])
        elseif isone(length(N[k]))
            A[1,i] = k
            A[2,i] = N[k]
            A[3:end,i] .= NaN
        else
            println("key $k didn't save correctly")
        end
    end
    writedlm(filename,A,',')
end

"""

```julia
dict2csv(filename::String,D::Dict)
```
Save a `Dict` to a .csv. Accounts for fields of any type, keys must be `Symbol`s.

For single-element entries, the rest of the array/table is filled with `NaN`s.

see also: [`nt2csv`](@ref), [`data2csv`](@ref)
"""
dict2csv(filename::String,D::Dict) = nt2csv(filename,(; D...))


"""

```julia
data2csv(filename::String,data)
```

Save `data` to a .csv file.
`data` may be one of the following: `Dict`, `NamedTuple`, or an `Array`,
Accepts single-element or Vector fields.
For single-element entries, the rest of the array/table is filled with `NaN`s.

see also: `nt2csv`, `dict2csv`

"""
function data2csv(filename::String,data)
    if isa(data,Dict)
        dict2csv(filename,data)
    elseif isa(data,NamedTuple)
        nt2csv(filename,data)
    elseif isa(data,AbstractArray)
        writedlm(filename,data,',')
    else
        throw(ArgumentError("Input data must be an Array, Dict, or NamedTuple"))
    end
end

## Load data from a csv to a NamedTuple or Dict

"""

```julia
csv2nt(filename::String;symbol=true)
```

Read a .csv file into a `NamedTuple`.
Assumes columns reflect discrete entries, with the first row specifying that element's key.
For NaN-buffered columns, the `NaN`s are removed, returning a single-element entry.

A `true` value for `symbol` will convert data saved as arrays of `String`s into arrays of `Symbol`s.

see also: [`nt2csv`](@ref), [`csv2dict`](@ref)

"""
function csv2nt(filename::String;symbol::Bool=true)
    A = readdlm(filename,',') #Input Array
    kz = Tuple(Symbol.(A[1,:]))
    pN = Vector{Union{AbstractVector,Number}}(undef,length(kz)) #proto-NamedTuple Vector{Vector}
    for i = 1:length(kz)
        T = typeof(A[2,i])
        if T <: Number
            if isnan(A[3,i])
                pN[i] = float(A[2,i])
            else
                pN[i] = convert(Vector{T},A[2:end,i])
            end
        elseif T <: AbstractString
            symbol ? pN[i] = convert(Vector{Symbol},Symbol.(A[2:end,i])) : pN[i] = A[2:end,i]
        else
            k = kz[i]
            println("key $k did not match any criteria"); flush(stdout)
        end
    end
    return NamedTuple{kz}(pN)
end

"""

```julia
csv2dict(filename::String;symbol::Bool=true)
```

Read a .csv file into a `Dict`.
Assumes columns reflect discrete entries, with the first row specifying that entry's key.
For NaN-buffered columns, the `NaN`s are removed, returning a single-element entry.

A `true` value for `symbol` will convert data saved as arrays of `String`s into arrays of `Symbol`s.

see also: [`nt2csv`](@ref), [`csv2nt`](@ref)

"""
function csv2dict(filename::String;symbol::Bool=true)
    dnt = csv2nt(filename,symbol=symbol)
    return Dict(k => dnt[k] for k ∈ keys(dnt))
end

"""

```julia
serial2dict(file::String, vars::Tuple; n, ll=true, accept=true, perturbation=false)
```

Load a serialized archive file from [`thermochron_metropolis`](@ref) into a `Dict` for post-run analysis, only incorporating the variables in `vars` and covering `n` steps (all steps by default). 
Provide a `Bool` to incorporate log-likelihoods `ll`, acceptances `accept`, or the perturbed variables `perturbation` at each step.

"""
function serial2dict(file::String,vars::Tuple;n::Integer=0,ll::Bool=true,accept::Bool=true,perturbation::Bool=false)

    data = deserialize(file)
    iszero(n) && (n=length(data.llDist))


    out = Dict{Symbol,Any}((vars[i],data.pDist[1:n,i]) for i ∈ 1:length(vars))

    ll && (out[:ll] = data.llDist[1:n])
    accept && (out[:accept] = data.acceptanceDist[1:n])
    perturbation && (out[:prt] = data.prt[1:n])
    return out
end
