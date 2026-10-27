# Chow classes and intersections of complete-intersection presentations.

@doc raw"""
    fundamental_class(Y::AbsToricCompleteIntersection; check_codimension=true)

Return the class of the fundamental cycle of a chosen complete-intersection
presentation in a simplicial (not necessarily complete) toric variety.
The product of the defining divisor classes represents this cycle subject to
the codimension and Cartier conditions.

The expensive codimension check can be disabled with `check_codimension=false`
only when the equations are known to cut out a nonempty complete intersection.
"""
function fundamental_class(
  Y::AbsToricCompleteIntersection; check_codimension::Bool=true,
)
  X = toric_variety(Y)
  @req is_simplicial(X) "Chow classes require a simplicial toric variety"
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

@doc raw"""
    rational_equivalence_class(Y::AbsToricCompleteIntersection; check_codimension=true)

Forward to [`fundamental_class`](@ref) with the same keyword argument.
"""
function rational_equivalence_class(
  Y::AbsToricCompleteIntersection; check_codimension::Bool=true,
)
  return fundamental_class(Y; check_codimension)
end

function Base.:*(ac::RationalEquivalenceClass, Y::AbsToricCompleteIntersection)
  @req toric_variety(ac) === toric_variety(Y) "Both objects must have the same ambient variety"
  return ac * fundamental_class(Y)
end
Base.:*(Y::AbsToricCompleteIntersection, ac::RationalEquivalenceClass) = ac * Y
