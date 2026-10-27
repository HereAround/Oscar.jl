# Properties and validation of toric complete-intersection presentations.

@doc raw"""
    is_empty(Y::ToricCompleteIntersection)

Test whether the saturated defining ideal is the unit ideal.
"""
@attr Bool is_empty(Y::ToricCompleteIntersection) = is_one(_saturated_toric_ideal(Y))
is_empty(Y::ToricHypersurface) = is_empty(_toric_presentation(Y))

@doc raw"""
    has_expected_codimension(Y::ToricCompleteIntersection)

Test whether the chosen equations cut out a nonempty subscheme of the expected
codimension, using saturation by the irrelevant ideal.
"""
@attr Bool function has_expected_codimension(Y::ToricCompleteIntersection)
  I = _saturated_toric_ideal(Y)
  return !is_one(I) &&
         dim(I) == ngens(cox_ring(toric_variety(Y))) - length(defining_equations(Y))
end
has_expected_codimension(Y::ToricHypersurface) =
  has_expected_codimension(_toric_presentation(Y))

function _check_toric_presentation(
  Y::AbsToricCompleteIntersection;
  check_smoothness::Bool=false, check_codimension::Bool=false,
)
  if check_codimension
    X = toric_variety(Y)
    @req is_smooth(X) ||
         all(D -> is_cartier(toric_divisor(D)), defining_divisor_classes(Y)) "Certified complete intersections require Cartier defining classes"
    @req has_expected_codimension(Y) "The equations do not have the expected codimension"
  end
  if check_smoothness
    @req is_smooth(Y) "The defined subscheme is singular"
  end
  return nothing
end
