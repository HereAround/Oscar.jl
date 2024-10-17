# This file contains code related to G-sets
# The idea of the implementation is that available GAP functionality
# for computing orbits, permutations, actions, etc. can be used,
# but that local replacements by pure Julia code (independent of GAP)
# are welcome.

import Hecke.orbit


# G-sets are "sets" (in a very general sense, these do not need to be objects of type `Set`)
# with an action by a group G::T.
# Alternatively one can define G-sets as a union of G-orbits.
# Potential examples include:
# - orbits of integers under a permutation
# - conjugacy classes of group elements
# - conjugacy classes of subgroups
# - block system


# TODO: add lots of concrete subtypes constructors, e.g. for
# - regular action of a group on itself
# - action of a perm group on its moved points
# - ...


"""
    GSetByElements{T,S} <: GSet{T,S}

Objects of this type represent G-sets that are willing to write down
orbits and elements lists as vectors.
These G-sets are created by default by [`gset`](@ref).

The fields are
- the group that acts, of type `T`,
- the Julia function (for example `on_tuples`) that describes the action,
- the seeds (something iterable of eltype `S`) whose closure under the action is the G-set
- the dictionary used to store attributes (orbits, elements, ...).
"""
@attributes mutable struct GSetByElements{T,S} <: GSet{T,S}
    group::T
    action_function::Function
    seeds

    function GSetByElements(G::T, fun::Function, seeds; closed::Bool = false) where {T<:Union{GAPGroup, FinGenAbGroup}}
        @assert !isempty(seeds)
        Omega = new{T,eltype(seeds)}(G, fun, seeds, Dict{Symbol,Any}())
        closed && set_attribute!(Omega, :elements => unique!(collect(seeds)))
        return Omega
    end
end
#TODO: How can I specify that `seeds` should be an iterable object?

function Base.show(io::IO, ::MIME"text/plain", x::GSetByElements)
  println(io, "G-set of")
  io = pretty(io)
  print(io, Indent())
  println(io, Lowercase(), x.group)
  print(io, "with seeds ", x.seeds)
  print(io, Dedent())
end

function Base.show(io::IO, x::GSetByElements)
  if is_terse(io)
    print(io, "G-set")
  else
    print(io, "G-set of ")
    io = pretty(io)
    print(terse(io), Lowercase(), x.group, " with seeds ", x.seeds)
  end
end

"""
    acting_group(Omega::GSetByElements)

Return the group `G` acting on `Omega`.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> acting_group(gset(G, [1])) == G
true
```
"""
acting_group(Omega::GSetByElements) = Omega.group

@doc raw"""
    action_function(Omega::GSetByElements)

Return the function $f: \Omega \times G \to \Omega$ that defines the G-set.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> action_function(gset(G, [1])) == ^
true

julia> action_function(gset(G, [[1, 2]])) == on_tuples
true

julia> action_function(gset(G, on_sets, [[1, 2]])) == on_sets
true
```
"""
action_function(Omega::GSetByElements) = Omega.action_function

# The following works for all G-set types that support attributes
# and for which the number of elements is an `Int`.
@attr PermGroup function action_range(Omega::GSet)
  return symmetric_group(length(Int, Omega))
end


#############################################################################
##
##  general method with explicit action function

"""
    gset(G::Union{GAPGroup, FinGenAbGroup}[, fun::Function], seeds, closed::Bool = false)

Return the G-set `Omega` that consists of the closure of the seeds `seeds`
under the action of `G` defined by `fun`.

This means that `Omega` contains all elements `fun(omega, g)`
for `omega` in `seeds` and `g` in `G`.

`fun` can be omitted if the element type of `seeds` implies
a reasonable default,
for example, if `G` is a `PermGroup` and `seeds` is a `Vector{T}`
where `T` is one of `Int`, `Set{Int}`, `Vector{Int}`.

If `closed` is set to `true` then `seeds` is assumed to be closed
under the action of `G`.
In this case, `collect(Omega)` is guaranteed to be equal to `collect(seeds)`;
in particular, the ordering of points in `seeds` (if applicable) is kept.
Note that the indexing of points in `Omega` is used by
[`action_homomorphism`](@ref).

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> length(gset(G, [1]))  # natural action
4

julia> length(gset(G, [[1, 2]]))  # action on ordered pairs
12

julia> length(gset(G, on_sets, [[1, 2]]))  # action on unordered pairs
6
```
"""
function gset(G::Union{GAPGroup, FinGenAbGroup}, fun::Function, seeds; closed::Bool = false)
  return GSetByElements(G, fun, seeds; closed = closed)
end


#############################################################################
##
##  G-sets where the action function can be omitted
##
##  (We use an indirection via `gset_by_type`, in order to admit specifying
##  a default action depending on the element type of `seeds` (which can be
##  any iterable collection.)

gset(G::T, seeds; closed::Bool = false) where T<:GAPGroup = gset_by_type(G, seeds, eltype(seeds); closed = closed)


## natural action of permutations on positive integers
function gset_by_type(G::PermGroup, Omega, ::Type{T}; closed::Bool = false) where T<:IntegerUnion
  return GSetByElements(G, ^, Omega; closed = closed)
end

## action of permutations on sets of positive integers
function gset_by_type(G::PermGroup, Omega, ::Type{T}; closed::Bool = false) where T<:Set{T2} where T2<:IntegerUnion
  return GSetByElements(G, on_sets, Omega; closed = closed)
end

## action of permutations on vectors of positive integers
function gset_by_type(G::PermGroup, Omega, ::Type{T}; closed::Bool = false) where T<:Vector{T2} where T2<:IntegerUnion
  return GSetByElements(G, on_tuples, Omega; closed = closed)
end

## action of permutations on tuples of positive integers
function gset_by_type(G::PermGroup, Omega, ::Type{T}; closed::Bool = false) where T<:Tuple{T2,Vararg{T2}} where T2<:IntegerUnion
  return GSetByElements(G, on_tuples, Omega; closed = closed)
end

## action of matrices on vectors via right multiplication
function gset_by_type(G::MatrixGroup{E, M}, Omega, ::Type{AbstractAlgebra.Generic.FreeModuleElem{E}}; closed::Bool = false) where E where M
  return GSetByElements(G, *, Omega; closed = closed)
end

## action of matrices on sets of vectors via right multiplication
function gset_by_type(G::MatrixGroup{E, M}, Omega, ::Type{T}; closed::Bool = false) where T <: Set{AbstractAlgebra.Generic.FreeModuleElem{E}} where E where M
  return GSetByElements(G, on_sets, Omega; closed = closed)
end

## action of matrices on vectors of vectors via right multiplication
function gset_by_type(G::MatrixGroup{E, M}, Omega, ::Type{T}; closed::Bool = false) where T <: Vector{AbstractAlgebra.Generic.FreeModuleElem{E}} where E where M
  return GSetByElements(G, on_tuples, Omega; closed = closed)
end

## action of matrices on subspaces via right multiplication
function gset_by_type(G::MatrixGroup{E, M}, Omega, ::Type{T}; closed::Bool = false) where T <: AbstractAlgebra.Generic.Submodule{E} where E where M
  return GSetByElements(G, ^, Omega; closed = closed)
end

## action of matrices on polynomials via `on_indeterminates`
function gset_by_type(G::MatrixGroup{E, M}, Omega, ::Type{T}; closed::Bool = false) where T <: MPolyRingElem{E} where E where M
  return GSetByElements(G, on_indeterminates, Omega; closed = closed)
end

## (add more such actions: on sets of sets, on sets of tuples, ...)

## natural action of a permutation group on the integers 1, ..., degree
gset(G::PermGroup) = gset(G, 1:G.deg; closed = true)

## natural action of a matrix group over a finite field on vectors
function gset(G::MatrixGroup{T, MT}) where T <: FinFieldElem where MT
    V = free_module(base_ring(G), degree(G))
    return gset(G, collect(V); closed = true)
end


#############################################################################
##
#TODO: Compute membership without writing down all elements,
#      using what is called `RepresentativeAction` in GAP.

function Base.in(omega::S, Omega::GSetByElements{T,S}) where {T,S}
    omega in Omega.seeds && return true
    return omega in elements(Omega)
end


#############################################################################
##
##  G-sets given by the complete set

function as_gset(G::T, fun::Function, Omega) where T<:Union{GAPGroup,FinGenAbGroup}
    return GSetByElements(G, fun, Omega; closed = true)
end

as_gset(G::T, Omega) where T<:Union{GAPGroup,FinGenAbGroup} = as_gset(G, ^, Omega)


#############################################################################
##
##  wrapper objects for elements of G-sets,
##  with fields `gset` (the G-set) and `objects` (the unwrapped object)
##
##  These objects are optional ("syntactic sugar"), they can be used to
##  - apply group elements via `^`,
##    not via the action function stored in the G-set,
##  - write something like `orbit(omega)`, `stabilizer(omega)`.

struct ElementOfGSet{T, S, G <: GSet{T, S}}
    gset::G
    obj::S
end

function (Omega::GSet{T, S})(obj::S) where {T, S}
    return ElementOfGSet(Omega, obj)
end

function ^(omega::ElementOfGSet, g::T) where {T<:AbstractAlgebra.GroupElem}
    Omega = omega.gset
    fun = action_function(Omega)
    return ElementOfGSet(Omega, fun(omega.obj, g))
end

==(omega1::ElementOfGSet, omega2::ElementOfGSet) =
  ((omega1.gset == omega2.gset) && (omega1.obj == omega2.obj))

function Base.hash(omega::ElementOfGSet, h::UInt)
  b = 0x4dd1b3e65edeab89 % UInt
  h = hash(omega.gset, h)
  h = hash(omega.obj, h)
  return xor(h, b)
end

Base.in(omega::ElementOfGSet, Omega::GSet) = Base.in(omega.obj, Omega)

Base.in(omega::ElementOfGSet, Omega::GSetByElements) = Base.in(omega.obj, Omega)

orbit(omega::ElementOfGSet) = orbit(omega.gset, omega.obj)


unwrap(omega::Any) = omega

unwrap(omega::ElementOfGSet) = omega.obj


#############################################################################
##
##  `:orbit`

"""
    orbit(G::Union{GAPGroup, FinGenAbGroup}[, fun::Function], omega)

Return the G-set that consists of the images of `omega`
under the action of `G` defined by `fun`.

This means that the result contains all elements `fun(omega, g)`
for `g` in `G`.

`fun` can be omitted if the type of `Omega` implies a reasonable default,
for example, if `G` is a `PermGroup` and `omega` is
one of `Int`, `Set{Int}`, `Vector{Int}`.

# Examples
```jldoctest
julia> G = symmetric_group(4);

julia> length(orbit(G, 1))
4

julia> length(orbit(G, [1, 2]))
12

julia> length(orbit(G, on_sets, [1, 2]))
6
```
"""
orbit(G::GAPGroup, omega) = gset_by_type(G, [omega], typeof(omega))

orbit(G::Union{GAPGroup, FinGenAbGroup}, fun::Function, omega) = GSetByElements(G, fun, [omega])

"""
    orbit(Omega::GSet, omega)

Return the G-set that consists of the elements `fun(omega, g)` where
`g` is in the group of `Omega` and `fun` is the underlying action of `Omega`.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Permutation group of degree 6 and order 16

julia> Omega = gset(G, [1, 5]);

julia> length(orbit(Omega, 1))
4
```
"""
function orbit(Omega::GSetByElements{<:GAPGroup}, omega::T) where T
    G = acting_group(Omega)
    acts = GapObj(gens(G))
    gfun = GapObj(action_function(Omega))

    # The following works only because GAP does not check
    # whether the given (dummy) group 'GapObj(G)' fits to the given generators,
    # or whether the elements of 'acts' are group elements.
    orb = Vector{T}(GAP.Globals.Orbit(GapObj(G), omega, acts, acts, gfun)::GapObj)

    res = as_gset(acting_group(Omega), action_function(Omega), orb)
    # We know that this G-set is transitive.
    set_attribute!(res, :orbits => [orb])
    return res
end
#T check whether omega lies in Omega?

function orbit(Omega::GSetByElements{FinGenAbGroup}, omega::T) where T
    return orbit_via_Julia(Omega, omega)
end

# simpleminded alternative directly in Julia
function orbit_via_Julia(Omega::GSet, omega)
    acts = gens(acting_group(Omega))
    orbarray = [omega]
    orb = Set(orbarray)
    fun = action_function(Omega)
    for p in orbarray
      for g in acts
        img = fun(p, g)
        if !(img in orb)
          push!(orbarray, img)
          push!(orb, img)
        end
      end
    end

    res = as_gset(acting_group(Omega), action_function(Omega), orbarray)
    # We know that this G-set is transitive.
    set_attribute!(res, :orbits => [orbarray])
    return res
end


#############################################################################
##
##  `:orbits` a vector of G-sets

"""
    orbits(Omega::GSet)

Return the vector of transitive G-sets in `Omega`.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Permutation group of degree 6 and order 16

julia> orbs = orbits(gset(G));

julia> map(collect, orbs)
2-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [5, 6]
```
"""
@attr Vector{GSetByElements{T,S}} function orbits(Omega::GSetByElements{T,S}) where {T <: Union{GAPGroup, FinGenAbGroup},S}
  orbs = GSetByElements{T,S}[]
  for p in Omega.seeds
    if all(o -> !(p in o), orbs)
      push!(orbs, orbit(Omega, p))
    end
  end
  return orbs
end

"""
    orbits(G::PermGroup)

Return the orbits of the natural G-set of `G`.

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Permutation group of degree 6 and order 16

julia> orbs = orbits(G);

julia> map(length, orbs)
2-element Vector{Int64}:
 4
 2
```
"""
@attr Vector{GSetByElements{PermGroup, Int}} orbits(G::PermGroup) = orbits(gset(G))


"""
    stabilizer(Omega::GSet{T,S})
    stabilizer(Omega::GSet{T,S}, omega::S = representative(Omega); check::Bool = true) where {T,S}

Return the subgroup of `G = acting_group(Omega)` that fixes `omega`,
together with the embedding of this subgroup into `G`.
If `check` is `false` then it is not checked whether `omega` is in `Omega`.

# Examples
```jldoctest
julia> Omega = gset(symmetric_group(3));

julia> stabilizer(Omega)
(Permutation group of degree 3 and order 2, Hom: permutation group -> Sym(3))
```
"""
@attr Tuple{sub_type(T), Map{sub_type(T), T}} function stabilizer(Omega::GSet{T,S}) where {T,S}
    return stabilizer(Omega, representative(Omega), check = false)
end

function stabilizer(Omega::GSet{T,S}, omega::S; check::Bool = true) where {T,S}
    check && @req omega in Omega "omega must be an element of Omega"
    G = acting_group(Omega)
    gfun = action_function(Omega)
    return stabilizer(G, omega, gfun)
end


#############################################################################
##
##  `:elements` a vector of points;
##  if `:seeds` is known to be closed under the action then
##  keep its ordering of points

@attr Vector{S} function elements(Omega::GSetByElements{T,S}) where {T,S}
  orbs = orbits(Omega)
  return union(map(collect, orbs)...)
end


#############################################################################
##

# In fact, '<:GAPGroup' is not used at all in this function.
"""
    permutation(Omega::GSetByElements{T}, g::BasicGAPGroupElem{T}) where T<:GAPGroup

Return the element of the permutation group that describes the action
of `g` on `Omega`, where `g` is an element of `acting_group(Omega)`.

# Examples

```jldoctest
julia> G = symmetric_group(4);

julia> Omega = gset(G, [[1, 2]]);

julia> x = gen(G, 1)
(1,2,3,4)

julia> permutation(Omega, x)
(1,2,4,7)(3,6,9,12)(5,8,10,11)
```
"""
function permutation(Omega::GSetByElements{T}, g::Union{GAPGroupElem, FinGenAbGroupElem}) where T<:Union{GAPGroup, FinGenAbGroup}
    omega_list = GAP.Obj(elements(Omega))
    gfun = GAP.Obj(action_function(Omega))

    # The following works only because GAP does not check
    # whether the given group element 'g' is a group element.
    pi = GAP.Globals.PermutationOp(g, omega_list, gfun)

    return group_element(action_range(Omega), pi)
end


@doc raw"""
    GSetBySubgroupTransversal{T, S, E} <: GSet{T}

Objects of this type represent G-sets that describe the left or right cosets
of a subgroup $H$ in a group $G$.
The group $G$ acts on the G-set by multiplication from the right or (after
taking inverses) from the left.
These G-sets store just transversals,
see [`right_transversal`](@ref) and [`left_transversal`](@ref).
The construction of explicit right or left cosets is not necessary in order
to compute the permutation action of elements of $G$ on the cosets.

The fields are
- the group that acts, of type `T`, with elements of type `E`,
- the subgroup whose cosets are the elements, of type `S`,
- the side from which the group acts (`:right` or `:left`),
- the (left or right) transversal, of type `SubgroupTransversal{T, S, E}`,
- the dictionary used to store attributes (orbits, elements, ...).
"""
@attributes mutable struct GSetBySubgroupTransversal{T, S, E} <: GSet{T,GroupCoset{T, E}}
    group::T
    subgroup::S
    side::Symbol
    transversal::SubgroupTransversal{T, S, E}

    function GSetBySubgroupTransversal(G::T, H::S, side::Symbol; check::Bool = true) where {T<:GAPGroup, S<:GAPGroup}
        check && @req is_subgroup(H, G)[1] "H must be a subgroup of G"
        E = eltype(G)
        if side == :right
          tr = right_transversal(G, H)
        elseif side == :left
          tr = left_transversal(G, H)
        else
          throw(ArgumentError("side must be :right or :left"))
        end
        return new{T, S, E}(G, H, side, tr, Dict{Symbol,Any}())
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::GSetBySubgroupTransversal)
  side = (x.side == :right ? "Right" : "Left")
  println(io, "$side cosets of")
  io = pretty(io)
  print(io, Indent())
  println(io, Lowercase(), x.subgroup, " in")
  print(io, Lowercase(), x.group)
  print(io, Dedent())
end

function Base.show(io::IO, x::GSetBySubgroupTransversal)
  side = (x.side == :right ? "Right" : "Left")
  if is_terse(io)
    print(io, "$side cosets of groups")
  else
    print(io, "$side cosets of ")
    io = pretty(io)
    print(terse(io), Lowercase(), x.subgroup, " in ", Lowercase(), x.group)
  end
end

acting_group(Omega::GSetBySubgroupTransversal) = Omega.group
action_function(Omega::GSetBySubgroupTransversal) = ((Omega.side == :right) ? (Base.:*) : function(omega, g) return inv(g)*omega; end)

function Base.in(omega::GroupCoset, Omega::GSetBySubgroupTransversal)
    return omega.side == Omega.side &&
           omega.G == Omega.group && omega.H == Omega.subgroup
end

Base.length(Omega::GSetBySubgroupTransversal) = index(Int, Omega.group, Omega.subgroup)
Base.length(::Type{T}, Omega::GSetBySubgroupTransversal) where T <: IntegerUnion = index(T, Omega.group, Omega.subgroup)

Base.lastindex(Omega::GSetBySubgroupTransversal) = length(Omega)

Base.keys(Omega::GSetBySubgroupTransversal) = keys(1:length(Omega))

function representative(Omega::GSetBySubgroupTransversal)
  if Omega.side == :right
    return right_coset(Omega.subgroup, one(Omega.group))
  else
    return left_coset(Omega.subgroup, one(Omega.group))
  end
end

function Base.iterate(Omega::GSetBySubgroupTransversal, state = 1)
  T = Omega.transversal
  state > length(T) && return nothing
  if Omega.side == :right
    return (right_coset(Omega.subgroup, T[state]), state+1)
  else
    return (left_coset(Omega.subgroup, T[state]), state+1)
  end
end

Base.eltype(::Type{GSetBySubgroupTransversal{T, S, E}}) where {S, T, E} = GroupCoset{T, E}

function Base.getindex(Omega::GSetBySubgroupTransversal, i::Int)
  if Omega.side == :right
    return right_coset(Omega.subgroup, Omega.transversal[i])
  else
    return left_coset(Omega.subgroup, Omega.transversal[i])
  end
end

is_transitive(Omega::GSetBySubgroupTransversal) = true

function orbit(G::T, omega::GroupCoset{T, S}) where T <: GAPGroup where S
    @req G == omega.G "omega must be a left or right coset in G"
    return GSetBySubgroupTransversal(G, omega.H, omega.side, check = false)
end
# We could admit the more general `is_subset(G, omega.G)`.
# One problem would be that `omega` would not be a point in the orbit,
# according to the definition of equality for cosets.

function orbit(Omega::GSetBySubgroupTransversal{T, S, E}, omega::GroupCoset{T, E}) where T <: GAPGroup where S <: GAPGroup where E
  @req (Omega.group == omega.G && Omega.subgroup == omega.H && Omega.side == omega.side) "omega is not in Omega"
  return Omega
end

orbits(Omega::GSetBySubgroupTransversal) = [Omega]

function permutation(Omega::GSetBySubgroupTransversal{T, S, E}, g::E) where T <: GAPGroup where S <: GAPGroup where E
  # The following works because GAP uses its `PositionCanonical`.
  # Note that we use `GAP.Globals.OnRight` also for the case of
  # a left transversal, since a right transversal is used on the GAP side.
  pi = GAP.Globals.PermutationOp(GapObj(g), Omega.transversal.X, GAP.Globals.OnRight)::GapObj
  return group_element(action_range(Omega), pi)
end

@attr GAPGroupHomomorphism{T, PermGroup} function action_homomorphism(Omega::GSetBySubgroupTransversal{T, S, E}) where T <: GAPGroup where S <: GAPGroup where E
  G = Omega.group

  # The following works because GAP uses its `PositionCanonical`.
  # Note that we use `GAP.Globals.OnRight` also for the case of
  # a left transversal, since a right transversal is used on the GAP side.
  acthom = GAP.Globals.ActionHomomorphism(GapObj(G), Omega.transversal.X, GAP.Globals.OnRight)::GapObj

  # See the comment about `SetJuliaData` in the `action_homomorphism` method
  # for `GSetByElements`.
  GAP.Globals.SetJuliaData(acthom, GAP.Obj([Omega, G]))

  return GAPGroupHomomorphism(G, action_range(Omega), acthom)
end


############################################################################
##
##  action homomorphisms

"""
    action_homomorphism(Omega::GSetByElements{T}) where T<:GAPGroup

Return the group homomorphism `act` with domain `G = acting_group(Omega)`
and codomain `symmetric_group(n)` that describes the permutation action
of `G` on `Omega`, where `Omega` has `n` elements.

This means that if an element `g` in `G` maps `collect(Omega)[i]` to
`collect(Omega)[j]` then `act(g)` maps `i` to `j`.

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> Omega = gset(G, [Set([1, 2])]);  # action on unordered pairs

julia> acthom = action_homomorphism(Omega)
Group homomorphism
  from Sym(6)
  to Sym(15)

julia> g = gen(G, 1)
(1,2,3,4,5,6)

julia> elms = collect(Omega);

julia> actg = acthom(g)
(1,2,3,5,7,10)(4,6,8,11,14,13)(9,12,15)

julia> elms[1]^g == elms[2]
true

julia> 1^actg == 2
true
```
"""
@attr GAPGroupHomomorphism{T, PermGroup} function action_homomorphism(Omega::GSetByElements{T}) where T<:GAPGroup
  G = acting_group(Omega)
  omega_list = GAP.Obj(collect(Omega))
  gap_gens = GapObj(gens(G); recursive = true)
  gfun = GAP.Obj(action_function(Omega))

  # The following works only because GAP does not check
  # whether the given generators in GAP and Julia fit together.
  acthom = GAP.Globals.ActionHomomorphism(GapObj(G), omega_list, gap_gens, GAP.Obj(gens(G)), gfun)::GapObj

  # The first difficulty on the GAP side is `ImagesRepresentative`
  # (which is the easy direction of the action homomorphism):
  # In the method in question, GAP does not really know how to compute
  # the group element that actually acts from the given group element;
  # there is only a rudimentary `FunctionAction` inside the
  # `UnderlyingExternalSet` of the GAP homomorphism object `acthom`.
  # We could replace this function here,
  # but this would introduce overhead for mapping each point.
  # Thus we install a special `ImagesRepresentative` method in GAP;
  # note that we know how to get the Julia "actor" from the GAP group
  # element, by wrapping it into the corresponding Julia group element.
  # (Yes, this is also overhead.
  # The alternative would be to create a new type of Oscar homomorphism,
  # which uses `permutation` or something better for mapping elements.)
  GAP.Globals.SetJuliaData(acthom, GAP.Obj([Omega, G]))

  return GAPGroupHomomorphism(G, action_range(Omega), acthom)
end

# for convenience: create the G-set on the fly
# (Here we assume that `Omega` is closed, this is dangerous.)
function action_homomorphism(G::PermGroup, Omega)
  return action_homomorphism(gset_by_type(G, Omega, eltype(Omega); closed = true))
end

function action_homomorphism(G::PermGroup, fun::Function, Omega)
  return action_homomorphism(GSetByElements(G, fun, Omega, closed = true))
end


"""
    is_conjugate(Omega::GSet, omega1, omega2)

Return `true` if `omega1`, `omega2` are in the same orbit of `Omega`,
and `false` otherwise.
To also obtain a conjugating element use [`is_conjugate_with_data`](@ref).

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Permutation group of degree 6 and order 16

julia> Omega = gset(G);

julia> is_conjugate(Omega, 1, 2)
true

julia> is_conjugate(Omega, 1, 5)
false
```
"""
is_conjugate(Omega::GSet, omega1, omega2) = omega2 in orbit(Omega, omega1)


"""
    is_conjugate_with_data(Omega::GSet, omega1, omega2)

Determine whether `omega1`, `omega2` are in the same orbit of `Omega`.
If yes, return `(true, g)` where `g` is an element in the group `G` of
`Omega` that maps `omega1` to `omega2`.
If not, return `(false, nothing)`.
If the conjugating element `g` is not needed, use [`is_conjugate`](@ref).

# Examples
```jldoctest
julia> G = sylow_subgroup(symmetric_group(6), 2)[1]
Permutation group of degree 6 and order 16

julia> Omega = gset(G);

julia> is_conjugate_with_data(Omega, 1, 2)
(true, (1,2))

julia> is_conjugate_with_data(Omega, 1, 5)
(false, ())
```
"""
function is_conjugate_with_data(Omega::GSet, omega1, omega2)
    # We do not call GAP's 'RepresentativeAction' with points, generators,
    # and actors.
    # The method in question would create a new 'ExternalSet' object
    # with a useless 'FunctionAction' value.
    # Instead, we delegate to the image of the action homomorphism.
    # (For that, we write down the elements of the G-set.
    # Computing the orbit of `omega1` or `omega2` would in principle suffice.)
    G = acting_group(Omega)
    acthom = action_homomorphism(Omega)
    elms = collect(Omega)
    pos1 = findfirst(isequal(omega1), elms)
    pos1 === nothing && return false, one(G)
    pos2 = findfirst(isequal(omega2), elms)
    pos2 === nothing && return false, one(G)
    img = GAP.Globals.RepresentativeAction(GapObj(image(acthom)[1]), pos1, pos2)
    img == GAP.Globals.fail && return false, one(G)
    pre = has_preimage_with_preimage(acthom, group_element(image(acthom)[1], img))
    @assert(pre[1])
    return true, pre[2]
end

############################################################################

acting_domain(Omega::GSet) = acting_group(Omega)

Base.length(Omega::GSetByElements) = length(elements(Omega))
Base.length(::Type{T}, Omega::GSetByElements) where T <: IntegerUnion = T(length(elements(Omega)))

representative(Omega::GSetByElements) = first(Omega.seeds)

function Base.iterate(Omega::GSetByElements, state = 1)
  elms = elements(Omega)
  state > length(elms) && return nothing
  return (elms[state], state+1)
end

Base.eltype(::Type{GSetByElements{T,S}}) where {T,S} = S

Base.getindex(Omega::GSetByElements, i::Int) = elements(Omega)[i]

blocks(G::GSet) = error("not implemented")
maximal_blocks(G::GSet) = error("not implemented")
minimal_block_reps(G::GSet) = error("not implemented")
all_blocks(G::GSet) = error("not implemented")

function is_transitive(Omega::GSetByElements)
    length(Omega.seeds) == 1 && return true
    return length(orbits(Omega)) == 1
end

is_regular(Omega::GSet) = is_transitive(Omega) && length(Omega) == order(acting_group(Omega))

function is_semiregular(Omega::GSet)
    ord = order(acting_group(Omega))
    return all(orb -> length(orb) == ord, orbits(Omega))
end

is_primitive(G::GSet) = error("not implemented")


"""
    blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))

Return a G-set that is a block system for the action of `G` on `L`,
i.e., a non-trivial partition of `L` preserved by the action of `G`.

Here, `L` must be a subvector of `1:degree(G)` on which `G` acts transitively.
`G` may move points outside `L`, in this case the restriction of the action
of the set stabilizer of `L` in `G` to `L` is considered.

An exception is thrown if this action is not transitive.

# Examples
```jldoctest
julia> g = sylow_subgroup(symmetric_group(4), 2)[1]
Permutation group of degree 4 and order 8

julia> collect(blocks(g))
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3, 4]
```
"""
function blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
   @assert is_transitive(G, L) "The group action is not transitive"
   bl = Vector{Vector{Int}}(GAP.Globals.Blocks(GapObj(G), GapObj(L))::GapObj)
   return gset(G, on_sets, bl; closed = true)
end

"""
    maximal_blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))

Return a G-set that is a maximal block system for the action of `G` on `L`,
i.e., a maximal non-trivial partition of `L` preserved by the action of `G`.

Here, `L` must be a subvector of `1:degree(G)` on which `G` acts transitively.
`G` may move points outside `L`, in this case the restriction of the action
of the set stabilizer of `L` in `G` to `L` is considered.

An exception is thrown if this action is not transitive.

# Examples
```jldoctest
julia> G = transitive_group(8, 2)
Permutation group of degree 8

julia> collect(maximal_blocks(G))
2-element Vector{Vector{Int64}}:
 [1, 2, 3, 8]
 [4, 5, 6, 7]
```
"""
function maximal_blocks(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
   @assert is_transitive(G, L) "The group action is not transitive"
   bl = Vector{Vector{Int}}(GAP.Globals.MaximalBlocks(GapObj(G), GapObj(L))::GapObj)
   return gset(G, bl; closed = true)
end


"""
    minimal_block_reps(G::PermGroup, L::AbstractVector{Int} = moved_points(G))

Return a vector of block representatives for all minimal non-trivial block
systems for the action of `G` on `L`.

Here, `L` must be a subvector of `1:degree(G)` on which `G` acts transitively.
`G` may move points outside `L`, in this case the restriction of the action
of the set stabilizer of `L` in `G` to `L` is considered.

An exception is thrown if this action is not transitive.

# Examples
```jldoctest
julia> G = transitive_group(8, 2)
Permutation group of degree 8

julia> minimal_block_reps(G)
3-element Vector{Vector{Int64}}:
 [1, 3]
 [1, 5]
 [1, 7]
```
"""
function minimal_block_reps(G::PermGroup, L::AbstractVector{Int} = moved_points(G))
   @assert is_transitive(G, L) "The group action is not transitive"
   return Vector{Vector{Int}}(GAP.Globals.RepresentativesMinimalBlocks(GapObj(G), GapObj(L))::GapObj)
end


"""
    all_blocks(G::PermGroup)

Return a vector of smallest representatives of all block systems
for the action of `G` on the set of moved points of `G`.

# Examples
```jldoctest
julia> G = transitive_group(8, 2)
Permutation group of degree 8

julia> all_blocks(G)
6-element Vector{Vector{Int64}}:
 [1, 2, 3, 8]
 [1, 5]
 [1, 3, 5, 7]
 [1, 3]
 [1, 3, 4, 6]
 [1, 7]
```
"""
all_blocks(G::PermGroup) = Vector{Vector{Int}}(GAP.Globals.AllBlocks(GapObj(G)))
#TODO: Do we really want to act on the set of moved points?


"""
    rank_action(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return the rank of the transitive action of `G` on `L`.
This is defined as the number of `G`-orbits in the action on ordered pairs
of points in `L`,
and is equal to the number of orbits of the stabilizer of a point in `L`
on `L`, see [Cam99](@cite) Section 1.11.

An exception is thrown if `G` is not transitive on `L`.

# Examples
```jldoctest
julia> G = symmetric_group(4); rank_action(G)  # 4-transitive
2

julia> H = sylow_subgroup(G, 2)[1]
Permutation group of degree 4 and order 8

julia> rank_action(H)  # not 2-transitive
3

julia> K = stabilizer(G, 1)[1]
Permutation group of degree 4 and order 6

julia> rank_action(K, 2:4)  # 2-transitive
2

julia> rank_action(K, 3:5)
ERROR: ArgumentError: the group is not transitive
```
"""
function rank_action(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
   @req is_transitive(G, L) "the group is not transitive"
   @req length(L) != 0 "the action domain is empty"
   H = stabilizer(G, L[1])[1]
   return length(orbits(gset(H, L, closed = true)))
end

"""
    transitivity(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return the maximum `k` such that `G` acts `k`-transitively on `L`,
that is, every `k`-tuple of points in `L` can be mapped simultaneously
to every other `k`-tuple by an element of `G`.

The output is `0` if `G` acts intransitively on `L`,
and an exception is thrown if `G` does not act on `L`.

# Examples
```jldoctest
julia> transitivity(mathieu_group(24))
5

julia> transitivity(symmetric_group(6))
6

julia> transitivity(symmetric_group(6), 1:7)
0

julia> transitivity(symmetric_group(6), 1:5)
ERROR: ArgumentError: the group does not act
```
"""
function transitivity(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))
  gL = GapObj(L)
  res = GAP.Globals.Transitivity(GapObj(G), gL)::Int
  @req res !== GAP.Globals.fail "the group does not act"
  # If the result is `0` then it may be that `G` does not act on `L`,
  # and in this case we want to throw an exception.
  if res == 0 && length(L) > 0
    lens = GAP.Globals.OrbitLengths(GapObj(G), gL)
#TODO: Compute the orbit lengths more efficiently than GAP does.
    @req sum(lens) == length(L) "the group does not act"
  end
  return res
end

"""
    is_transitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether `G` acts transitively on `L`, that is,
`L` is an orbit of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> is_transitive(G)
true

julia> is_transitive(sylow_subgroup(G, 2)[1])
false

julia> is_transitive(stabilizer(G, 1)[1])
false
```
"""
is_transitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsTransitive(GapObj(G), GapObj(L))
# Note that this definition does not coincide with that of the
# property `GAP.Globals.IsTransitive`, for which the default domain
# of the action is the set of moved points.


"""
    is_primitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether the action of `G` on `L` is primitive, that is,
the action is transitive and the point stabilizers are maximal in `G`.

# Examples
```jldoctest
julia> G = alternating_group(6);

julia> mx = filter(is_transitive, map(representative, maximal_subgroup_classes(G)))
3-element Vector{PermGroup}:
 Permutation group of degree 6 and order 24
 Permutation group of degree 6 and order 36
 Permutation group of degree 6 and order 60

julia> [(order(H), is_primitive(H)) for H in mx]
3-element Vector{Tuple{ZZRingElem, Bool}}:
 (24, 0)
 (36, 0)
 (60, 1)
```
"""
is_primitive(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsPrimitive(GapObj(G), GapObj(L))


"""
    is_regular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether the action of `G` on `L` is regular
(i.e., transitive and semiregular).

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> H = sub(G, [G([2, 3, 4, 5, 6, 1])])[1]
Permutation group of degree 6

julia> is_regular(H)
true

julia> is_regular(G)
false
```
"""
is_regular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsRegular(GapObj(G), GapObj(L))


"""
    is_semiregular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G))

Return whether the action of `G` on `L` is semiregular
(i.e., the stabilizer of each point is the identity).

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> H = sub(G, [G([2, 3, 1, 5, 6, 4])])[1]
Permutation group of degree 6

julia> is_semiregular(H)
true

julia> is_regular(H)
false
```
"""
is_semiregular(G::PermGroup, L::AbstractVector{Int} = 1:degree(G)) = GAPWrap.IsSemiRegular(GapObj(G), GapObj(L))

"""
    orbit_representatives_and_stabilizers(G::MatrixGroup{E}, k::Int) where E <: FinFieldElem

Return a vector of pairs `(orb, stab)` where `orb` runs over representatives
of orbits of `G` on the `k`-dimensional subspaces of `F^n`,
where `G` is a subgroup of `general_linear_group(F, n)`.

# Examples
```jldoctest
julia> G = orthogonal_group(1, 4, GF(3))
GO+(4,3)

julia> res = orbit_representatives_and_stabilizers(G, 1);

julia> length(res)
3

julia> print(sort([index(G, stab) for (U, stab) in res]))
ZZRingElem[12, 12, 16]
```
"""
function orbit_representatives_and_stabilizers(G::MatrixGroup{E}, k::Int) where E <: FinFieldElem
  F = base_ring(G)
  n = degree(G)
  q = GAP.Obj(order(F))
  V = vector_space(F, n)
  orbs = GAP.Globals.Orbits(GapObj(G), GAP.Globals.Subspaces(GAPWrap.GF(q)^n, k))
  orbreps = [GAP.Globals.BasisVectors(GAPWrap.Basis(orb[1])) for orb in orbs]
  stabs = [Oscar._as_subgroup_bare(G, GAP.Globals.Stabilizer(GapObj(G), v, GAP.Globals.OnSubspacesByCanonicalBasis)) for v in orbreps]::Vector{typeof(G)}
  orbreps1 = [[[F(x) for x in v] for v in bas] for bas in orbreps]::Vector{Vector{Vector{elem_type(F)}}}
  orbreps2 = [sub(V, [V(v) for v in bas])[1] for bas in orbreps1]
  return [(orbreps2[i], stabs[i]) for i in 1:length(stabs)]
end
