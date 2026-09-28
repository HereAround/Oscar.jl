##############################################################
# 1: Closed subschemes presented by homogeneous Cox ideals
##############################################################

@attributes mutable struct ClosedSubschemeOfToricVariety
  toric_variety::NormalToricVarietyType
  defining_ideal::MPolyIdeal
  function ClosedSubschemeOfToricVariety(
    toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal
  )
    @req is_simplicial(toric_variety) "Closed subschemes are currently supported only for simplicial toric varieties"
    @req base_ring(defining_ideal) === cox_ring(toric_variety) "The defining ideal must belong to the ambient Cox ring"
    # OSCAR has no ideal-level homogeneity predicate; check the stored generators.
    @req all(is_homogeneous, gens(defining_ideal)) "The defining ideal must be homogeneous in the Cox grading"
    return new(toric_variety, defining_ideal)
  end
end

##############################################################
# 2: Generic constructor
##############################################################

@doc raw"""
    closed_subscheme_of_toric_variety(toric_variety::NormalToricVarietyType, defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}})

Construct a closed subscheme of a simplicial normal toric variety.
The polynomials must be homogeneous in the ambient Cox ring. Their ideal is
retained as presentation data; its associated ideal sheaf defines the subscheme.
No radical is taken, and reducedness, irreducibility, and smoothness are not
required. Use [`covered_scheme`](@ref) for the associated scheme.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> closed_subscheme_of_toric_variety(f2, [t1])
Closed subscheme of a normal toric variety
```
"""
closed_subscheme_of_toric_variety(
  toric_variety::NormalToricVarietyType,
  defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem,QQMPolyRingElem}},
) = ClosedSubschemeOfToricVariety(
  toric_variety, ideal(defining_polynomials)
)

@doc raw"""
    closed_subscheme_of_toric_variety(toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal)

Construct a closed subscheme of a simplicial normal toric variety.
The ideal must belong to the ambient Cox ring and have homogeneous generators;
the current generator list is checked at construction (OSCAR has no
ideal-level homogeneity predicate). The original ideal is retained,
without taking its radical or saturation. Its associated ideal sheaf defines
the subscheme, which can be singular, reducible, or nonreduced. Different Cox
ideals can define the same associated scheme. Use [`covered_scheme`](@ref) for
scheme properties and [`defining_ideal`](@ref) for the stored presentation.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> closed_subscheme_of_toric_variety(f2, ideal([t1]))
Closed subscheme of a normal toric variety
```
"""
closed_subscheme_of_toric_variety(
  toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal
) = ClosedSubschemeOfToricVariety(
  toric_variety, defining_ideal
)

######################
# 3: Display
######################s

function Base.show(io::IO, c::ClosedSubschemeOfToricVariety)
  properties_string = ["Closed"]
  push_attribute_if_exists!(properties_string, c, :is_empty, "empty")
  push!(properties_string, "subscheme of a normal toric variety")
  join(io, properties_string, ", ", " ")
end
