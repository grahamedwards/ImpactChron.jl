# Efficiently save and load data with DelimitedFiles functions.

## save data to a csv:

"""
```julia
nt2csv(filename::String,N::NamedTuple)
```

Save a `NamedTuple` to a .csv. Accounts for fields of any type.
For single-element entries, the rest of the array/table is filled with `NaN`s.

see also: `dict2csv`, `data2csv`
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
Save a `Dict` to a .csv. Accounts for fields of any type.
For single-element entries, the rest of the array/table is filled with `NaN`s.

see also: `nt2csv`, `data2csv`
"""
function dict2csv(filename::String,D::Dict)
    nt2csv(filename,(; D...))
end


"""
```julia
data2csv(filename::String,data)
```

Save `data` to a .csv file.
`data` may be a composite type of `Dict` or `NamedTuple` or an `Array`.
Accepts fields of any type element.
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
csv2nt(filename::String;symbol::Bool=true)
```

Read a .csv file into a `NamedTuple`.
Assumes columns reflect individual elements, with the first row specifying that element's key.
For NaN-buffered columns (see `nt2csv`), the `NaN`s are removed, returning a single-element entry.

A `true` value for `symbol` will convert data saved as arrays of `String`s into arrays of `Symbol`s.

see also: `nt2csv`, `csv2dict`
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

Read a .csv file into a `NamedTuple`.
Assumes columns reflect individual entries, with the first row specifying that entry's key.
For NaN-buffered columns (see `nt2csv`), the `NaN`s are removed, returning a single-element entry.

A `true` value for `symbol` will convert data saved as arrays of `String`s into arrays of `Symbol`s.

see also: `nt2csv`, `csv2dict`
"""
function csv2dict(filename::String;symbol::Bool=true)
    dnt = csv2nt(filename,symbol=symbol)
    return Dict(k => dnt[k] for k ∈ keys(dnt))
end


## Add NetCDF capability some day...
