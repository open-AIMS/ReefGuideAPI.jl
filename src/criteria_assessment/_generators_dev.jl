
abstract type DeploymentCriteria end

"""
    generate_criteria(name::Symbol, field_spec::Dict)
    generate_criteria(name::Symbol, fields::Vector{Symbol})
    generate_criteria(name::Symbol)

Generate structs for a given criteria
"""
function generate_criteria(name::Symbol, field_spec::Dict)
    fields = [:($(Symbol(k))::$(typeof(v)) = $(v)) for (k, v) in field_spec]
    @eval @kwdef struct $(name)
        $(fields...)
    end
end
function generate_criteria(name::Symbol, fields::Vector{Symbol})
    # fields = [:($(Symbol(field))) for field in field_spec]
    @eval struct $(name)
        $(fields...)
    end
end
function generate_criteria(name::Symbol)
    pascal_case_name = Symbol(name, :Criteria)

    @eval struct $(pascal_case_name){T<:Union{Int64,Float64}} <: DeploymentCriteria
        data::Raster
        lower_bound::T
        upper_bound::T
    end

    # Data-only constructor
    @eval function $(pascal_case_name)(data)
        return $(pascal_case_name)(data, minimum(data), maximum(data))
    end

    @eval export $(pascal_case_name)
end

function generate_criteria_structs()::Nothing
    patterns = criteria_data_map()

    # Create structs for each criteria
    for k in keys(patterns)
        generate_criteria(k)
    end
end

# function generate_criteria_state()::Nothing
# end
