################################################
# 1: Intersections involving complete-intersection presentations
################################################

function Base.:*(ac::RationalEquivalenceClass, Y::AbsToricCompleteIntersection)
  @req toric_variety(ac) === toric_variety(Y) "The rational equivalence class and the complete intersection must be defined on the same toric variety"
  return ac * rational_equivalence_class(Y)
end

function Base.:*(Y::AbsToricCompleteIntersection, ac::RationalEquivalenceClass)
  @req toric_variety(ac) === toric_variety(Y) "The rational equivalence class and the complete intersection must be defined on the same toric variety"
  return ac * rational_equivalence_class(Y)
end
