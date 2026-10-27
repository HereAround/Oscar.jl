# Chow classes and intersections of complete-intersection presentations.

@doc raw"""
    rational_equivalence_class(Y::AbsToricCompleteIntersection; check_codimension=true)

Return the toric Chow class of a chosen complete intersection from the divisor
classes of its defining equations. This is not a general conversion of Cox
ideals to scheme-theoretic cycles.
"""
function rational_equivalence_class(
  Y::AbsToricCompleteIntersection; check_codimension::Bool=true,
)
  X = toric_variety(Y)
  @req is_simplicial(X) && is_complete(X) "Chow classes require a complete simplicial toric variety"
  if check_codimension
    _check_toric_presentation(Y; check_codimension=true)
  else
    P = _toric_presentation(Y)
    @req !has_attribute(P, :has_expected_codimension) ||
         get_attribute(P, :has_expected_codimension) "The equations are known not to have the expected codimension"
  end
  result = RationalEquivalenceClass(X, one(chow_ring(X)))
  for D in defining_divisor_classes(Y)
    result *= rational_equivalence_class(D)
  end
  return result
end

function Base.:*(ac::RationalEquivalenceClass, Y::AbsToricCompleteIntersection)
  @req toric_variety(ac) === toric_variety(Y) "Both objects must have the same ambient variety"
  return ac * rational_equivalence_class(Y)
end
Base.:*(Y::AbsToricCompleteIntersection, ac::RationalEquivalenceClass) = ac * Y
