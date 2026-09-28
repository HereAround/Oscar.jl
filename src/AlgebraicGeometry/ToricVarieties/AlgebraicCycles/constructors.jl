##############################################
# 1: The Julia type for toric algebraic cycles
##############################################

@attributes mutable struct RationalEquivalenceClass
  v::NormalToricVarietyType
  p::MPolyQuoRingElem
  RationalEquivalenceClass(v::NormalToricVarietyType, p::MPolyQuoRingElem) = new(v, p)
end

####################################################
# 2: Generic constructors
####################################################

@doc raw"""
    rational_equivalence_class(v::NormalToricVarietyType, p::MPolyQuoRingElem)

Construct the rational equivalence class of algebraic cycles corresponding to a linear combination of cones.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> chow_ring(P2)
Quotient
  of multivariate polynomial ring in 3 variables x1, x2, x3
    over rational field
  by ideal (x1 - x3, x2 - x3, x1*x2*x3)

julia> (x1, x2, x3) = gens(chow_ring(P2))
3-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x1
 x2
 x3

julia> rational_equivalence_class(P2, x1)
Rational equivalence classon a normal toric variety represented by V(x3)
```
"""
function rational_equivalence_class(v::NormalToricVarietyType, p::MPolyQuoRingElem)
  @req (is_simplicial(v) && is_complete(v)) "Currently, algebraic cycles are only supported for toric varieties that are simplicial and complete"
  @req parent(p) == chow_ring(v) "The polynomial must reside in the Chow ring of the toric variety"
  return RationalEquivalenceClass(v, p)
end

@doc raw"""
    rational_equivalence_class(v::NormalToricVarietyType, coefficients::Vector{T}) where {T <: IntegerUnion}

Construct the rational equivalence class of algebraic cycles corresponding to a linear combination of cones.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> rational_equivalence_class(P2, [6, 5, 4, 3, 2, 1])
Rational equivalence class on a normal toric variety represented by 15V(x1,x3)+6V(x3)
```
"""
function rational_equivalence_class(
  v::NormalToricVarietyType, coefficients::Vector{T}
) where {T<:IntegerUnion}
  @req (is_simplicial(v) && is_complete(v)) "Currently, algebraic cycles are only supported for toric varieties that are simplicial and complete"
  @req length(coefficients) == n_cones(v; trivial=false) "The number of coefficients must match the number of all cones (but the trivial one) in the fan of the toric variety"
  mons = gens_of_rational_equivalence_classes(v)
  return RationalEquivalenceClass(
    v, sum(coefficients[i] * mons[i] for i in 1:length(coefficients))
  )
end

####################################################
# 3: Special constructors
####################################################

@doc raw"""
    rational_equivalence_class(d::ToricDivisor)

Construct the rational equivalence class of algebraic cycles corresponding to the toric divisor `d`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> d = toric_divisor(P2, [1, 2, 3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)
```
"""
function rational_equivalence_class(d::ToricDivisor)
  v = toric_variety(d)
  if is_trivial(d)
    return RationalEquivalenceClass(v, zero(chow_ring(v)))
  end
  coeffs = coefficients(d)
  indets = gens(chow_ring(v))
  return RationalEquivalenceClass(v, sum(coeffs[k] * indets[k] for k in 1:length(coeffs)))
end

@doc raw"""
    rational_equivalence_class(c::ToricDivisorClass)

Construct the algebraic cycle corresponding to the toric divisor class `c`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> tdc = toric_divisor_class(P2, [2])
Divisor class on a normal toric variety

julia> rational_equivalence_class(tdc)
Rational equivalence class on a normal toric variety represented by 2V(x3)
```
"""
rational_equivalence_class(c::ToricDivisorClass) = rational_equivalence_class(
  toric_divisor(c)
)

@doc raw"""
    RationalEquivalenceClass(l::ToricLineBundle)

Construct the toric algebraic cycle corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(P2, [2])
Toric line bundle on a normal toric variety

julia> polynomial(rational_equivalence_class(l))
2*x3
```
"""
rational_equivalence_class(l::ToricLineBundle) = rational_equivalence_class(
  toric_divisor(l)
)

@doc raw"""
    rational_equivalence_class(cc::CohomologyClass)

Construct the toric algebraic cycle corresponding to the cohomology class `cc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> (x1, x2, x3) = gens(cohomology_ring(P2))
3-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 x1
 x2
 x3

julia> cc = cohomology_class(P2, x1+x2)
Cohomology class on a normal toric variety given by x1 + x2

julia> rational_equivalence_class(cc)
Rational equivalence class on a normal toric variety represented by 2V(x3)
```
"""
rational_equivalence_class(cc::CohomologyClass) = RationalEquivalenceClass(
  toric_variety(cc), polynomial(chow_ring(toric_variety(cc)), cc)
)

@doc raw"""
    rational_equivalence_class(Y::AbsToricCompleteIntersection;
                               check_codimension=true)

Construct the toric Chow class of the chosen complete-intersection
presentation of `Y` by multiplying the classes of its defining Cartier
divisors. The ambient toric variety must be complete and simplicial. By
default, the expected-codimension/Cartier check is performed; set
`check_codimension=false` only when those complete-intersection hypotheses
are being supplied as an explicit assumption.

This method is intentionally not defined for an arbitrary
`ClosedSubschemeOfToricVariety`. Its stored Cox ideal is presentation data,
and OSCAR currently has no general conversion from an intrinsic
scheme-theoretic cycle of an arbitrary covered scheme to a toric
`RationalEquivalenceClass`.
"""
function rational_equivalence_class(
  Y::AbsToricCompleteIntersection; check_codimension::Bool=true,
)
  X = toric_variety(Y)
  @req is_simplicial(X) && is_complete(X) "Currently, algebraic cycles are only supported for complete simplicial toric varieties"
  if check_codimension
    _check_toric_presentation(Y; check_codimension=true)
  else
    for Z in (Y, _toric_presentation(Y))
      @req !has_attribute(Z, :has_expected_codimension) || get_attribute(Z, :has_expected_codimension) "The equations are known not to have the expected codimension"
    end
  end
  result = RationalEquivalenceClass(X, one(chow_ring(X)))
  for D in defining_divisor_classes(Y)
    result = result * rational_equivalence_class(D)
  end
  return result
end

####################################################
# 4: Addition, subtraction and scalar multiplication
####################################################

function Base.:+(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
  @req toric_variety(ac1) === toric_variety(ac2) "The rational equivalence classes must be defined on the same toric variety"
  return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) + polynomial(ac2))
end

function Base.:-(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
  @req toric_variety(ac1) === toric_variety(ac2) "The rational equivalence classes must be defined on the same toric variety"
  return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) - polynomial(ac2))
end

Base.:*(c::QQFieldElem, ac::RationalEquivalenceClass) = RationalEquivalenceClass(
  toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac)
)
Base.:*(c::Rational{Int64}, ac::RationalEquivalenceClass) = RationalEquivalenceClass(
  toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac)
)
Base.:*(c::T, ac::RationalEquivalenceClass) where {T<:IntegerUnion} =
  RationalEquivalenceClass(
    toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac)
  )

####################################################
# 5: Intersection product
####################################################

function Base.:*(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
  @req toric_variety(ac1) === toric_variety(ac2) "The rational equivalence classes must be defined on the same toric variety"
  return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) * polynomial(ac2))
end

Base.:^(ac::RationalEquivalenceClass, p::T) where {T<:IntegerUnion} =
  RationalEquivalenceClass(
    toric_variety(ac), polynomial(ac)^p
  )

function Base.:*(ac::RationalEquivalenceClass, Y::AbsToricCompleteIntersection)
  @req toric_variety(ac) === toric_variety(Y) "The rational equivalence class and the complete intersection must be defined on the same toric variety"
  return ac * rational_equivalence_class(Y)
end

function Base.:*(Y::AbsToricCompleteIntersection, ac::RationalEquivalenceClass)
  @req toric_variety(ac) === toric_variety(Y) "The rational equivalence class and the complete intersection must be defined on the same toric variety"
  return ac * rational_equivalence_class(Y)
end

function Base.:*(Y1::AbsToricCompleteIntersection, Y2::AbsToricCompleteIntersection)
  @req toric_variety(Y1) === toric_variety(Y2) "The complete intersections must be defined on the same toric variety"
  return rational_equivalence_class(Y1) * rational_equivalence_class(Y2)
end

####################################################
# 6: Equality and hash
####################################################

function Base.:(==)(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
  return toric_variety(ac1) === toric_variety(ac2) && polynomial(ac1) == polynomial(ac2)
end

function Base.hash(ac::RationalEquivalenceClass, h::UInt)
  b = 0xb5d4ac6b9084eb6e % UInt
  h = hash(toric_variety(ac), h)
  h = hash(polynomial(ac), h)
  return xor(h, b)
end

####################################################
# 7: Display
####################################################

function Base.show(io::IO, ac::RationalEquivalenceClass)
  if is_trivial(ac)
    join(io, "Trivial rational equivalence class on a normal toric variety")
  else
    # otherwise, extract properties to represent the rational equivalence class
    r = representative(ac)
    coeffs = [c for c in AbstractAlgebra.coefficients(r)]
    expos = [
      matrix(ZZ, [k for k in AbstractAlgebra.exponent_vectors(m)]) for
      m in AbstractAlgebra.monomials(r)
    ]
    indets = gens(chow_ring(toric_variety(ac)))

    # form string to be printed
    properties_string = String[]
    for i in 1:length(coeffs)
      m = String[]
      for j in 1:ncols(expos[i])
        for k in 1:expos[i][1, j]
          push!(m, string(indets[j]))
        end
      end
      tmp = join(m, ",")
      if i == 1 && coeffs[i] == 1
        push!(
          properties_string,
          "Rational equivalence classon a normal toric variety represented by V($tmp)",
        )
      elseif i == 1 && coeffs[i] != 1
        push!(
          properties_string,
          "Rational equivalence class on a normal toric variety represented by $(coeffs[i])V($tmp)",
        )
      elseif i > 1 && coeffs[i] == 1
        push!(properties_string, "+V($tmp)")
      elseif i > 1 && coeffs[i] > 0 && coeffs[i] != 1
        push!(properties_string, "+$(coeffs[i])V($tmp)")
      elseif i > 1 && coeffs[i] < 0
        push!(properties_string, "$(coeffs[i])V($tmp)")
      end
    end
    # print information
    join(io, properties_string)
  end
end
