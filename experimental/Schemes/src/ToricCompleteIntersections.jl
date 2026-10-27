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


function _check_toric_euler_ambient(
  X::NormalToricVariety, classes::AbstractVector{ToricDivisorClass};
  completeness_check::Bool=true,
)
  @req is_smooth(X) "Intersection theory requires a smooth ambient toric variety"
  if completeness_check
    @req is_complete(X) "Intersection theory requires a complete ambient toric variety"
  end
  @req 1 <= length(classes) <= dim(X) "Provide between one and dim(X) divisor classes"
  @req all(D -> toric_variety(D) === X, classes) "Divisor classes must belong to the ambient toric variety"
  return nothing
end

function _check_toric_euler_presentation(
  Y::AbsToricCompleteIntersection;
  completeness_check::Bool=true, check_smoothness::Bool=false,
  check_codimension::Bool=false,
)
  _check_toric_euler_ambient(toric_variety(Y), defining_divisor_classes(Y); completeness_check)
  _check_toric_presentation(Y; check_smoothness, check_codimension)
  P = _toric_presentation(Y)
  # A known failure is never silently overridden by an unchecked Euler call.
  @req !has_attribute(P, :is_smooth) ||
       get_attribute(P, :is_smooth) "The hypersurface is known to be singular"
  if has_attribute(P, :covered_scheme)
    S = get_attribute(P, :covered_scheme)::AbsCoveredScheme
    @req !has_attribute(S, :is_smooth) ||
         get_attribute(S, :is_smooth) "The hypersurface is known to be singular"
  end
  @req !has_attribute(P, :has_expected_codimension) ||
       get_attribute(P, :has_expected_codimension) "The equations are known not to have the expected codimension"
  return nothing
end

# Formal adjunction, truncated at the desired degree in the Chow ring.
function _toric_ci_chern_polynomial(
  X::NormalToricVariety, classes::AbstractVector{ToricDivisorClass}, k::Int;
  completeness_check::Bool=true,
)
  cs = [polynomial(chern_class(X, j; completeness_check)) for j in 0:k]
  for D in classes
    d = polynomial(cohomology_class(D; completeness_check))
    for j in 1:k
      cs[j + 1] -= d * cs[j]
    end
  end
  return cs[k + 1]
end

@doc raw"""
    chern_class(Y::AbsToricCompleteIntersection, k::Int;
                completeness_check=true, check_smoothness=false, check_codimension=false)

Return an ambient cohomology class restricting to the k-th Chern class of a
smooth complete intersection. The ambient toric variety must be smooth and
complete, and the equations must form a regular sequence. Unless explicitly
checked, hypersurface smoothness and expected codimension are assumptions.
"""
function chern_class(
  Y::AbsToricCompleteIntersection, k::Int;
  completeness_check::Bool=true, check_smoothness::Bool=false,
  check_codimension::Bool=false,
)
  X = toric_variety(Y)
  classes = defining_divisor_classes(Y)
  @req 0 <= k <= dim(X) - length(classes) "Chern class index is outside the expected dimension"
  _check_toric_euler_presentation(Y; completeness_check, check_smoothness, check_codimension)
  p = _toric_ci_chern_polynomial(X, classes, k; completeness_check)
  return cohomology_class(X, p; completeness_check)
end

@doc raw"""
    euler_characteristic(X::NormalToricVariety, classes::AbstractVector{ToricDivisorClass};
                         completeness_check=true, algorithm=:intersection_theory)

Compute the Euler characteristic of a smooth member of the expected
complete-intersection codimension with the given divisor classes in a smooth,
complete toric ambient variety. No equations are supplied or checked.
This is separate from the existing Euler characteristic of the toric variety X.
"""
function euler_characteristic(
  X::NormalToricVariety, classes::AbstractVector{ToricDivisorClass};
  algorithm::Symbol=:intersection_theory, completeness_check::Bool=true,
)
  @req algorithm == :intersection_theory "Unknown Euler characteristic algorithm"
  _check_toric_euler_ambient(X, classes; completeness_check)
  k = dim(X) - length(classes)
  p = _toric_ci_chern_polynomial(X, classes, k; completeness_check)
  for D in classes
    p *= polynomial(cohomology_class(D; completeness_check))
  end
  value = integrate(cohomology_class(X, p; completeness_check); completeness_check)
  @req isone(denominator(value)) "The computed Euler characteristic is not integral"
  return ZZ(value)
end

@doc raw"""
    euler_characteristic(Y::AbsToricCompleteIntersection;
                         algorithm=:intersection_theory, completeness_check=true,
                         check_smoothness=false, check_codimension=false)

Compute the Euler characteristic of a smooth complete-intersection presentation
using adjunction and ambient intersection theory. The ambient variety must be
smooth and complete. Smoothness and regularity of the equations are assumptions
unless checked. Previously known failures are rejected even when checks are
disabled. Arbitrary or singular closed subschemes are not supported.
"""
function euler_characteristic(
  Y::AbsToricCompleteIntersection;
  algorithm::Symbol=:intersection_theory, completeness_check::Bool=true,
  check_smoothness::Bool=false, check_codimension::Bool=false,
)
  @req algorithm == :intersection_theory "Unknown Euler characteristic algorithm"
  _check_toric_euler_presentation(Y; completeness_check, check_smoothness, check_codimension)
  return euler_characteristic(
    toric_variety(Y), defining_divisor_classes(Y); algorithm, completeness_check
  )
end
