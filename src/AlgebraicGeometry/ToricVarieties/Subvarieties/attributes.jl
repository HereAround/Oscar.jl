############################
# Attributes
############################

@doc raw"""
    toric_variety(c::ClosedSubschemeOfToricVariety)

Return the ambient normal toric variety supplied at construction.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subscheme_of_toric_variety(f2, [t1])
Closed subscheme of a normal toric variety

julia> toric_variety(c) == f2
true
```
"""
@attr NormalToricVarietyType toric_variety(c::ClosedSubschemeOfToricVariety) =
  c.toric_variety

@doc raw"""
    defining_ideal(c::ClosedSubschemeOfToricVariety)

Return the original homogeneous Cox-ring ideal supplied at construction.
This is presentation data: different ideals can define the same ideal sheaf
and hence the same closed subscheme. No radical or saturation is taken.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subscheme_of_toric_variety(f2, [t1])
Closed subscheme of a normal toric variety

julia> defining_ideal(c) == ideal([t1])
true
```
"""
@attr MPolyIdeal defining_ideal(c::ClosedSubschemeOfToricVariety) = c.defining_ideal

@doc raw"""
    radical(c::ClosedSubschemeOfToricVariety)

Return the radical of the stored Cox-ring ideal. This is neither a canonical
ideal of the associated subscheme nor its scheme-theoretic reduction, and it
does not test reducedness of that scheme. Use [`is_reduced`](@ref) on the
subscheme for the latter. No saturation is performed.

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subscheme_of_toric_variety(f2, [t1])
Closed subscheme of a normal toric variety

julia> radical(c) == ideal([t1])
true
```
"""
@attr MPolyIdeal function radical(c::ClosedSubschemeOfToricVariety)
  return radical(defining_ideal(c))
end
