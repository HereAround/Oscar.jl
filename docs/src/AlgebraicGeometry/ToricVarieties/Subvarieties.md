# Closed Subschemes of Toric Varieties

## Cox ideals and associated schemes

A homogeneous ideal in the Cox ring of a simplicial normal toric variety
presents a closed subscheme, which may be singular, reducible, or nonreduced.
`closed_subscheme_of_toric_variety(X, I)` retains the ambient variety and the
original homogeneous ideal `I`. Its intrinsic scheme is `covered_scheme(c)`,
constructed from `ideal_sheaf(X, I)`. No radical or saturation is taken when
storing the presentation. Homogeneity is checked at construction.

The set-theoretic and scheme-theoretic statements in [CLS11](@cite) differ.
Over the complex numbers, for a simplicial toric variety without torus factors,
Proposition 5.2.4 describes closed vanishing loci of homogeneous Cox ideals.
These loci depend only on radical information. Proposition 5.2.7 gives their
bijection with radical homogeneous ideals contained in the irrelevant ideal
``B(\Sigma)``. These statements do not describe nilpotent structure.

In the Chapter 6 appendix, again assuming no torus factors, sheafification of a
homogeneous ideal defines a closed subscheme; Proposition 6.A.6 says every
closed subscheme arises this way. The stored Cox ideal is not unique.
For a **smooth** toric variety without torus factors, Proposition 6.A.7 identifies
equality of the associated ideal sheaves with equality of the saturations
``I:B(\Sigma)^\infty``. This equivalence is not asserted here for singular
simplicial varieties. OSCAR's scheme computations use the associated ideal
sheaf over the ambient coefficient field.

OSCAR uses *variety* for a geometrically integral scheme. These subschemes
have no such requirement. The legacy names `ClosedSubvarietyOfToricVariety`
and `closed_subvariety_of_toric_variety` are deprecated aliases for the
canonical subscheme API; they introduce no geometric-integrality check.

## General Closed Subschemes

### Constructors

```@docs
closed_subscheme_of_toric_variety(toric_variety::NormalToricVarietyType, defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}})
closed_subscheme_of_toric_variety(toric_variety::NormalToricVarietyType, defining_ideal::MPolyIdeal)
```

### Conversion to a scheme

Use `covered_scheme(Y)` for the associated OSCAR scheme of a closed subscheme,
hypersurface, or complete intersection. The conversion respects the irrelevant
locus and retains any nonreduced structure after sheafification.

```@docs
covered_scheme(::ClosedSubschemeOfToricVariety)
covered_scheme(::AbsToricCompleteIntersection)
```

### Properties

Dimension, smoothness, reducedness, irreducibility, and integrality refer to the
associated scheme. In particular, a nonradical Cox ideal can present a reduced
scheme if its nonreduced structure is confined to the irrelevant locus.
Emptiness is tested by saturation with the irrelevant ideal. Unsupported scheme
computations propagate errors. Scheme checks can be computationally expensive.

```@docs
dim(::ClosedSubschemeOfToricVariety)
dim(::AbsToricCompleteIntersection)
is_empty(::ClosedSubschemeOfToricVariety)
is_empty(::AbsToricCompleteIntersection)
is_smooth(::ClosedSubschemeOfToricVariety)
is_smooth(::AbsToricCompleteIntersection)
is_reduced(::ClosedSubschemeOfToricVariety)
is_reduced(::AbsToricCompleteIntersection)
is_irreducible(::ClosedSubschemeOfToricVariety)
is_irreducible(::AbsToricCompleteIntersection)
is_integral(::ClosedSubschemeOfToricVariety)
is_integral(::AbsToricCompleteIntersection)
```

### Presentation data

`defining_ideal(Y)` returns the stored Cox ideal. `radical(Y)` returns its
ordinary ideal-theoretic radical, without saturation. Neither accessor supplies
an intrinsic defining ideal or computes a scheme-theoretic reduction.
For reductions of affine charts, OSCAR provides `reduced_scheme(U)`; reducedness
of the entire associated scheme is tested by `is_reduced(Y)`.

```@docs
toric_variety(::ClosedSubschemeOfToricVariety)
toric_variety(::AbsToricCompleteIntersection)
defining_ideal(::ClosedSubschemeOfToricVariety)
defining_ideal(::AbsToricCompleteIntersection)
radical(::ClosedSubschemeOfToricVariety)
radical(::AbsToricCompleteIntersection)
```

## Hypersurfaces and Complete Intersections

A toric hypersurface or complete intersection retains its chosen homogeneous
Cox equations and divisor classes together with the underlying closed subscheme.
Here *toric* describes the ambient variety. The subscheme itself need not be
toric, smooth, reduced, or irreducible.

The constructors check the complete-intersection conditions by default:
Cartier defining classes and a nonempty subscheme of codimension equal to the
number of equations. The corresponding local equations then form a regular
sequence, since a normal toric variety is Cohen--Macaulay. On a smooth ambient
variety every divisor class is Cartier. Non-Cartier presentations on singular
ambient varieties are not certified by this check.

`check_codimension=false` explicitly assumes these conditions and leaves them
unverified. Such an unchecked object may fail to be a complete intersection;
use `closed_subscheme_of_toric_variety` for arbitrary homogeneous ideals.
Smoothness is a separate, optional check, requested with `check_smoothness=true`.
Skipping it never caches a positive smoothness result.

### Constructors

```@docs
toric_hypersurface(::Union{NormalToricVariety,ClosedSubschemeOfToricVariety}, ::MPolyDecRingElem)
toric_complete_intersection(::NormalToricVariety, ::AbstractVector)
```

When given an existing closed subscheme, the supplied equations must generate
its stored Cox ideal exactly. A hypersurface therefore requires an explicit
single generator; the constructor does not infer one from arbitrary ideal
generators. Equality of ideal sheaves alone is not used to replace the stored
presentation. `AbstractVector` permits vector views as input; every entry must
be a homogeneous polynomial in the ambient Cox ring.

### Attributes

`closed_subscheme_of_toric_variety(Y)` returns the stored subscheme wrapper,
sharing its scheme computations. `defining_equations(Y)` and
`defining_divisor_classes(Y)` also work for hypersurfaces, returning one-element
vectors. The hypersurface wrapper supports the single-equation accessor
`defining_equation(Y)` and shares the complete-intersection algorithms.

```@docs
closed_subscheme_of_toric_variety(::AbsToricCompleteIntersection)
defining_equation(::ToricHypersurface)
defining_equations(::AbsToricCompleteIntersection)
defining_divisor_classes(::AbsToricCompleteIntersection)
```

### Codimension

`has_expected_codimension(Y)` tests nonemptiness and the dimension condition.
For non-Cartier defining classes on a singular ambient variety, this predicate
alone does not establish that the presentation is a regular sequence.

```@docs
has_expected_codimension
```

### Euler characteristics by intersection theory

For a smooth complete intersection of codimension ``r`` in a smooth, complete
toric variety ``X``, the engine computes
```math
\chi(Y) = \int_X \left[c(TX)\prod_{i=1}^r \frac{D_i}{1+D_i}\right]_{\dim X}.
```
The inverse factors are truncated to the required degree. All arithmetic is
exact, and the result is an integer in `ZZ`.

For example, a Fermat quartic surface has Euler characteristic 24:
```julia
X = projective_space(NormalToricVariety, 3)
f = sum(x^4 for x in gens(cox_ring(X)))
Y = toric_hypersurface(X, f)
euler_characteristic(Y) # 24, assuming smoothness
```

The constructors check the complete-intersection conditions by default. Euler
computations assume smoothness unless `check_smoothness=true` is requested.
For an object constructed with `check_codimension=false`, pass
`check_codimension=true` to the Euler method to verify that condition too.
Skipped checks remain unverified. Ambient smoothness is
always required; `completeness_check=false` skips only completeness verification.
Checks and algorithm selection are performed even when cached data exist.

The class-only method `euler_characteristic(X, classes)` instead computes a
numerical invariant of a smooth member of the expected codimension, **assuming
such a member exists**. It does not create or certify any actual subvariety,
and must not be interpreted as the Euler characteristic of a singular member.

```@docs
chern_class(::AbsToricCompleteIntersection, ::Int)
euler_characteristic(::AbsToricCompleteIntersection)
euler_characteristic(::NormalToricVariety, ::AbstractVector{ToricDivisorClass})
```
