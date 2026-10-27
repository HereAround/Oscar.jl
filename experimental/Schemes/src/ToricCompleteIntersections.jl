# A chosen regular-sequence presentation, not an arbitrary Cox ideal.
abstract type AbsToricCompleteIntersection end

@attributes mutable struct ToricCompleteIntersection <: AbsToricCompleteIntersection
  ambient::NormalToricVariety
  equations::Vector{MPolyDecRingElem{QQFieldElem,QQMPolyRingElem}}

  function ToricCompleteIntersection(
    X::NormalToricVariety, equations::AbstractVector;
    check_codimension::Bool=true, check_smoothness::Bool=false,
  )
    fs = _toric_equations(X, equations)
    Y = new(X, fs)
    _check_toric_presentation(Y; check_codimension, check_smoothness)
    return Y
  end
end

@attributes mutable struct ToricHypersurface <: AbsToricCompleteIntersection
  presentation::ToricCompleteIntersection

  function ToricHypersurface(presentation::ToricCompleteIntersection)
    @req length(defining_equations(presentation)) == 1 "A hypersurface requires one equation"
    return new(presentation)
  end
end

_toric_presentation(Y::ToricCompleteIntersection) = Y
_toric_presentation(Y::ToricHypersurface) = Y.presentation
toric_variety(Y::AbsToricCompleteIntersection) = _toric_presentation(Y).ambient
defining_equations(Y::AbsToricCompleteIntersection) = copy(_toric_presentation(Y).equations)
defining_equation(Y::ToricHypersurface) = only(defining_equations(Y))

# Divisor classes are derived from the equations; inconsistent class data cannot
# be supplied through an alternate constructor.
defining_divisor_classes(Y::AbsToricCompleteIntersection) = [
  toric_divisor_class(toric_variety(Y), degree(f)) for f in defining_equations(Y)
]

@attr MPolyIdeal defining_ideal(Y::ToricCompleteIntersection) =
  ideal(cox_ring(toric_variety(Y)), defining_equations(Y))
defining_ideal(Y::ToricHypersurface) = defining_ideal(_toric_presentation(Y))

@attr MPolyIdeal function _saturated_toric_ideal(Y::ToricCompleteIntersection)
  return saturation(defining_ideal(Y), irrelevant_ideal(toric_variety(Y)))
end

@attr Bool is_empty(Y::ToricCompleteIntersection) = is_one(_saturated_toric_ideal(Y))
is_empty(Y::ToricHypersurface) = is_empty(_toric_presentation(Y))

# The scheme conversion belongs to experimental/Schemes, which owns the
# Cox-ideal to toric ideal-sheaf conversion and its affine-chart methods.
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

function _toric_equations(X::NormalToricVariety, equations::AbstractVector)
  @req 1 <= length(equations) <= dim(X) "Provide between one and dim(X) equations"
  R = cox_ring(X)
  @req all(f -> f isa MPolyDecRingElem && parent(f) === R, equations) "Equations must belong to the ambient Cox ring"
  @req all(f -> !iszero(f) && !is_constant(f), equations) "Defining equations must be nonzero and nonconstant"
  @req all(is_homogeneous, equations) "Defining equations must be homogeneous in the Cox grading"
  return elem_type(R)[copy(f) for f in equations]
end

@doc raw"""
    toric_complete_intersection(X::NormalToricVariety, equations::AbstractVector;
                                check_codimension=true, check_smoothness=false)

Represent a chosen complete-intersection presentation by homogeneous Cox-ring
equations and their divisor classes in a normal toric ambient variety. The
chosen equations, not an arbitrary ideal generating set, define the presentation.

By default, check the expected codimension using saturation by the irrelevant
ideal and an ideal dimension computation. This can be expensive; pass
check_codimension=false to assume the condition. On a singular ambient variety,
certification additionally requires Cartier defining classes. Smoothness is
checked separately on request, using the intrinsic covered scheme.
"""
function toric_complete_intersection(
  X::NormalToricVariety, equations::AbstractVector;
  check_codimension::Bool=true, check_smoothness::Bool=false,
)
  return ToricCompleteIntersection(
    X, equations; check_codimension, check_smoothness
  )
end

@doc raw"""
    toric_hypersurface(X::NormalToricVariety, equation::MPolyDecRingElem;
                       check_codimension=true, check_smoothness=false)

Construct a chosen toric hypersurface presentation from one homogeneous Cox
equation. The scheme may be singular or nonreduced. The expected-codimension
check is enabled by default; smoothness is checked only on request.
"""
function toric_hypersurface(
  X::NormalToricVariety, equation::MPolyDecRingElem;
  check_codimension::Bool=true, check_smoothness::Bool=false,
)
  return ToricHypersurface(toric_complete_intersection(
    X, [equation]; check_codimension, check_smoothness
  ))
end

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

Base.show(io::IO, ::ToricHypersurface) =
  print(io, "Hypersurface in a normal toric variety")
Base.show(io::IO, Y::ToricCompleteIntersection) =
  print(io, "Complete-intersection presentation with ",
        length(defining_equations(Y)), " equations in a normal toric variety")
