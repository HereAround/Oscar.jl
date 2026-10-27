# Toric hypersurfaces and chosen complete-intersection presentations.

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


Base.show(io::IO, ::ToricHypersurface) =
  print(io, "Hypersurface in a normal toric variety")
Base.show(io::IO, Y::ToricCompleteIntersection) =
  print(io, "Complete-intersection presentation with ",
        length(defining_equations(Y)), " equations in a normal toric variety")
