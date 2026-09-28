@doc raw"""
    chern_class(m::AbstractFTheoryModel, k::Int)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
this very cohomology class ``h``, otherwise it raises an error.

Smoothness and the expected codimension of the hypersurface are assumed unless
`check_smoothness=true` and `check_codimension=true`, respectively. Computation
delegates to the generic toric complete-intersection Chern-class method; existing
cached classes are retained, with requested checks performed before cache access.

!!! note "Completeness check"
    The implemented algorithm is guaranteed to work only for toric ambient spaces
    that are smooth and **complete**. Verifying completeness can be very time 
    consuming. To skip this check, pass the optional keyword argument 
    `completeness_check=false`.

!!! warning
    This method works ONLY for F-theory models which are hypersurfaces in a toric ambient space.
  
!!! warning
    This method represents the Chern classes of said hypersurface by cohomology classes on the toric ambient space.
    These classes counterparts must be restricted to the hypersurface to truly represent the Chern class in question.
    Internally, we integrate those ambient space classes against the class of the hypersurface, which automatically
    executes the restriction to the hypersurface.

# Examples
```jldoctest; setup = :(Oscar.ensure_qsmdb_installed())
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> h = chern_class(qsm_model, 4; completeness_check = false);

julia> is_trivial(h)
false
```
"""
function chern_class(
  m::AbstractFTheoryModel, k::Int;
  completeness_check::Bool=true, check_smoothness::Bool=false,
  check_codimension::Bool=false,
)
  Y = toric_hypersurface(m; check_codimension=false)
  @req 0 <= k <= dim(ambient_space(m)) - 1 "Chern class index must lie between zero and the expected dimension"
  _check_toric_euler_presentation(Y; completeness_check, check_smoothness, check_codimension)
  if has_attribute(m, :chern_classes)
    cs = get_attribute(m, :chern_classes)::Dict{Int64,CohomologyClass}
    haskey(cs, k) && return cs[k]
  else
    cs = Dict{Int64,CohomologyClass}()
  end
  cs[k] = chern_class(Y, k; completeness_check)
  set_attribute!(m, :chern_classes, cs)
  return cs[k]
end

@doc raw"""
    chern_classes(m::AbstractFTheoryModel)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
a vector with the cohomology classes corresponding to all non-trivial Chern classes ``c_k`` of
``Y_n``. Otherwise, this methods raises an error.

As of right now, this method is computationally expensive for involved toric ambient spaces,
such as in the example below.

!!! note "Completeness check"
    The implemented algorithm is guaranteed to work only for toric ambient spaces
    that are simplicial and **complete**. Verifying completeness can be very time 
    consuming. To skip this check, pass the optional keyword argument 
    `completeness_check=false`.

# Examples
```jldoctest; setup = :(Oscar.ensure_qsmdb_installed())
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> h = chern_classes(qsm_model; completeness_check = false);

julia> is_one(polynomial(h[0]))
true

julia> is_trivial(h[1])
true

julia> is_trivial(h[2])
false
```
"""
function chern_classes(m::AbstractFTheoryModel; completeness_check::Bool=true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric ambient space"
  for k in 0:(dim(ambient_space(m)) - 1)
    chern_class(m, k; completeness_check)
  end
  return get_attribute(m, :chern_classes)::Dict{Int64,CohomologyClass}
end

@doc raw"""
    euler_characteristic(m::AbstractFTheoryModel)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a smooth hypersurface of the expected codimension in a smooth, complete toric
ambient space, compute its Euler characteristic in any dimension by delegation
to [`toric_hypersurface`](@ref) and the generic intersection-theory engine.
The return type remains `Int` for compatibility. Set `check_smoothness=true`
and `check_codimension=true` to verify the subvariety assumptions; by default
these checks are skipped and smoothness remains unknown. A known singular
subvariety is rejected. `algorithm=:intersection_theory` selects the engine.
Cached or serialized Euler values do not bypass checks or algorithm selection;
the computed result refreshes the existing Euler attribute.

!!! note "Completeness check"
    The implemented algorithm is guaranteed to work only for toric ambient spaces
    that are smooth and **complete**. Verifying completeness can be very time 
    consuming. To skip this check, pass the optional keyword argument 
    `completeness_check=false`.

# Examples
```jldoctest; setup = :(Oscar.ensure_qsmdb_installed())
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> h = euler_characteristic(qsm_model; completeness_check = false)
378
```
"""
function euler_characteristic(
  m::AbstractFTheoryModel;
  algorithm::Symbol=:intersection_theory, completeness_check::Bool=true,
  check_smoothness::Bool=false, check_codimension::Bool=false,
)
  Y = toric_hypersurface(m; check_codimension=false)
  value = Int(euler_characteristic(Y; algorithm, completeness_check, check_smoothness, check_codimension))
  set_attribute!(m, :euler_characteristic, value)
  return value
end

@doc raw"""
    euler_characteristic(m::CompleteIntersectionModel, equations::AbstractVector;
                         algorithm=:intersection_theory, completeness_check=true,
                         check_smoothness=false, check_codimension=false)

Compute the Euler characteristic using the supplied complete-intersection
presentation and the generic toric intersection-theory engine. The equations
must generate the stored ideal. The ambient variety must be smooth and complete;
subvariety smoothness and expected codimension are assumed unless checked.
The existing model representation and serialized fields are unchanged.
"""
function euler_characteristic(
  m::CompleteIntersectionModel, equations::AbstractVector;
  algorithm::Symbol=:intersection_theory, completeness_check::Bool=true,
  check_smoothness::Bool=false, check_codimension::Bool=false,
)
  Y = toric_complete_intersection(m, equations; check_codimension=false)
  return Int(euler_characteristic(Y; algorithm, completeness_check, check_smoothness, check_codimension))
end
