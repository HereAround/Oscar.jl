######################
# Properties
######################

@doc raw"""
    is_empty(c::ClosedSubschemeOfToricVariety)

Check whether the associated subscheme is empty by saturating its Cox ideal
by the irrelevant ideal. For simplicial toric varieties without torus factors
over the complex numbers, see Proposition 5.2.6 in [CLS11](@cite).

# Examples
```jldoctest
julia> f2 = hirzebruch_surface(NormalToricVariety, 2);

julia> (t1, x1, t2, x2) = gens(cox_ring(f2));

julia> c = closed_subscheme_of_toric_variety(f2, [t1])
Closed subscheme of a normal toric variety

julia> is_empty(c)
false

julia> c2 = closed_subscheme_of_toric_variety(f2, [x1,x2])
Closed subscheme of a normal toric variety

julia> is_empty(c2)
true
```
"""
@attr Bool function is_empty(c::ClosedSubschemeOfToricVariety)
  B = irrelevant_ideal(toric_variety(c))
  return (is_one(defining_ideal(c)) || is_one(saturation(defining_ideal(c), B)))
end

@attr AbsCoveredScheme function _toric_subscheme(c::ClosedSubschemeOfToricVariety)
  I = defining_ideal(c)
  if is_empty(c)
    return CoveredScheme(coefficient_ring(toric_variety(c)))
  end
  return subscheme(ideal_sheaf(toric_variety(c), I))
end

@doc raw"""
    covered_scheme(c::ClosedSubschemeOfToricVariety)

Return the closed subscheme defined by `c` as an OSCAR covered scheme.
The construction uses the ideal sheaf on the ambient toric charts, excludes
the irrelevant locus, and preserves nonreduced structure. Empty subschemes
give an empty covered scheme over the coefficient ring of the ambient variety.

The result is cached and reused by dimension and smoothness computations.
Conversion itself does not check smoothness or expected codimension.
"""
covered_scheme(c::ClosedSubschemeOfToricVariety) = _toric_subscheme(c)

@doc raw"""
    is_smooth(c::ClosedSubschemeOfToricVariety)

Test smoothness of the associated scheme via its toric ideal sheaf. This uses
the toric charts (thus excludes the irrelevant locus) and retains nonreduced
scheme structure. The result is cached only after a successful computation;
unsupported computations and other errors propagate rather than returning false.
"""
@attr Bool is_smooth(c::ClosedSubschemeOfToricVariety) = is_smooth(covered_scheme(c))

@doc raw"""
    dim(c::ClosedSubschemeOfToricVariety)

Return the dimension of the associated covered scheme, or `-inf` if it is empty.
"""
dim(c::ClosedSubschemeOfToricVariety) = dim(covered_scheme(c))

@doc raw"""
    is_reduced(c::ClosedSubschemeOfToricVariety)

Test whether the associated covered scheme has no nonzero nilpotents.
This tests the scheme after sheafification, not radicality of the stored Cox ideal.
"""
is_reduced(c::ClosedSubschemeOfToricVariety) = is_reduced(covered_scheme(c))

@doc raw"""
    is_irreducible(c::ClosedSubschemeOfToricVariety)

Test irreducibility of the associated covered scheme. The empty scheme is not
irreducible. Unsupported scheme computations propagate errors.
"""
is_irreducible(c::ClosedSubschemeOfToricVariety) =
  !is_empty(c) && is_irreducible(covered_scheme(c))

@doc raw"""
    is_integral(c::ClosedSubschemeOfToricVariety)

Test whether the associated covered scheme is nonempty, reduced, and irreducible.
This does not assert geometric integrality, which OSCAR requires for varieties.
"""
is_integral(c::ClosedSubschemeOfToricVariety) = is_integral(covered_scheme(c))
