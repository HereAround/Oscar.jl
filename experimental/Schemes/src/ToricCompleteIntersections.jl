# Covered-scheme integration depends on the experimental Schemes machinery.
# The scheme conversion belongs to experimental/Schemes, which owns the
# Cox-ideal to toric ideal-sheaf conversion and its affine-chart methods.

@doc raw"""
    covered_scheme(Y::ToricCompleteIntersection)

Return the covered scheme determined by the homogeneous Cox ideal of `Y`.
"""
@attr AbsCoveredScheme function covered_scheme(Y::ToricCompleteIntersection)
  return subscheme(ideal_sheaf(toric_variety(Y), defining_ideal(Y)))
end
covered_scheme(Y::ToricHypersurface) = covered_scheme(_toric_presentation(Y))

dim(Y::AbsToricCompleteIntersection) = dim(covered_scheme(Y))
@attr Bool is_smooth(Y::ToricCompleteIntersection) = is_smooth(covered_scheme(Y))
is_smooth(Y::ToricHypersurface) = is_smooth(_toric_presentation(Y))
is_reduced(Y::AbsToricCompleteIntersection) = is_reduced(covered_scheme(Y))
is_irreducible(Y::AbsToricCompleteIntersection) = is_irreducible(covered_scheme(Y))
is_integral(Y::AbsToricCompleteIntersection) = is_integral(covered_scheme(Y))
