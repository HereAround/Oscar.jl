# Complete Intersections

## Introduction

A toric hypersurface or complete intersection is described by a chosen list of
homogeneous equations in the Cox ring of a normal toric variety. The types
`ToricHypersurface` and `ToricCompleteIntersection` retain these equations,
their ambient toric variety, and the resulting divisor classes.

Unlike a general `ClosedSubvarietyOfToricVariety`, these objects record a
specific complete-intersection presentation. By default, the constructors
check that the chosen equations have the expected codimension. These checks
may require an expensive saturation by the irrelevant ideal. The optional
`check_smoothness` keyword requests a separate smoothness test.

## Constructors

```@docs
toric_hypersurface(X::NormalToricVariety, equation::MPolyDecRingElem)
toric_complete_intersection(X::NormalToricVariety, equations::AbstractVector)
```

## Attributes

```@docs
toric_variety(Y::AbsToricCompleteIntersection)
defining_equations(Y::AbsToricCompleteIntersection)
defining_equation(Y::ToricHypersurface)
defining_divisor_classes(Y::AbsToricCompleteIntersection)
defining_ideal(Y::ToricCompleteIntersection)
```

The associated covered scheme, its dimension, and its geometric properties
are accessible through the experimental Schemes functionality.

## Properties

```@docs
is_empty(Y::ToricCompleteIntersection)
has_expected_codimension(Y::ToricCompleteIntersection)
```

## Rational equivalence class

For a simplicial (not necessarily complete) ambient toric variety, the product
of the defining divisor classes represents the fundamental cycle under the
complete-intersection assumptions.

```@docs
fundamental_class(Y::AbsToricCompleteIntersection)
rational_equivalence_class(Y::AbsToricCompleteIntersection)
```
