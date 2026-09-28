@doc raw"""
    base_space(m::AbstractFTheoryModel)

Return the base space of the F-theory model.

# Examples
```jldoctest
julia> using Random;

julia> w = weierstrass_model_over_projective_space(2, rng = Random.Xoshiro(1234));

julia> dim(base_space(w))
2
```
"""
function base_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) ||
    @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY base space.\n"
  return m.base_space
end

@doc raw"""
    ambient_space(m::AbstractFTheoryModel)

Return the ambient space of the F-theory model.

# Examples
```jldoctest
julia> using Random;

julia> w = weierstrass_model_over_projective_space(2, rng = Random.Xoshiro(1234));

julia> dim(ambient_space(w))
4
```
"""
function ambient_space(m::AbstractFTheoryModel)
  is_base_space_fully_specified(m) ||
    @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning AUXILIARY ambient space.\n"
  return m.ambient_space
end

@doc raw"""
    fiber_ambient_space(m::AbstractFTheoryModel)

Return the fiber ambient space of an F-theory model.

# Examples
```jldoctest
julia> using Random;

julia> w = weierstrass_model_over_projective_space(2, rng = Random.Xoshiro(1234));

julia> dim(fiber_ambient_space(w))
2
```
"""
function fiber_ambient_space(m::AbstractFTheoryModel)
  @req hasfield(typeof(m), :fiber_ambient_space) "fiber_ambient_space not supported for this F-theory model"
  return m.fiber_ambient_space
end

@doc raw"""
    model_index(m::AbstractFTheoryModel)
Return database index of a literature model. This index is a unique identifier that can be used to more conveniently construct the model. 
All models have a model_index and these will not change in the future.

# Examples
```jldoctest
julia> using Random

julia> t = literature_model(31, rng = Random.Xoshiro(1234))
Weierstrass model over a not fully specified base -- F-theory weierstrass model dual to hypersurface model with fiber ambient space F_10 based on arXiv paper 1408.4808 Eq. (3.130)

julia> model_index(t)
31
```
"""
function model_index(m::AbstractFTheoryModel)
  model_indices = _model_indices()
  return parse(Int, model_indices["model" * literature_identifier(m) * ".json"])
end

@doc raw"""
    calabi_yau_hypersurface(m::AbstractFTheoryModel)

Return the Calabi–Yau hypersurface that defines the hypersurface model
as a closed subscheme of its toric ambient space.

# Examples
```jldoctest
julia> using Random;

julia> w =  weierstrass_model_over_projective_space(3, rng = Random.Xoshiro(1234))
Weierstrass model over a concrete base

julia> calabi_yau_hypersurface(w)
Closed subscheme of a normal toric variety

julia> t = global_tate_model_over_projective_space(2, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base

julia> calabi_yau_hypersurface(t)
Closed subscheme of a normal toric variety

julia> b = projective_space(NormalToricVariety, 2);

julia> kb = anticanonical_divisor_class(b);

julia> fiber_ambient = weighted_projective_space(NormalToricVariety, [2, 3, 1]);

julia> set_coordinate_names(fiber_ambient, ["x", "y", "z"]);

julia> p = "x^3 - y^2 + x1^12 * x * z^4 + x2^18 * z^6 + 13 * x3^3*x*y*z";

julia> h = hypersurface_model(b, fiber_ambient, [2 * kb, 3 * kb], p; completeness_check=false)
Hypersurface model over a concrete base

julia> calabi_yau_hypersurface(h)
Closed subscheme of a normal toric variety
```
"""
@attr ClosedSubschemeOfToricVariety function calabi_yau_hypersurface(
  m::AbstractFTheoryModel
)
  @req (m isa HypersurfaceModel || m isa WeierstrassModel || m isa GlobalTateModel) "Calabi-Yau hypersurface currently not supported for this model"
  @req base_space(m) isa NormalToricVariety "Calabi-Yau hypersurface currently only supported for toric varieties as base space"
  is_base_space_fully_specified(m) ||
    @vprint :FTheoryModelPrinter 1 "Base space was not fully specified. Returning hypersurface in AUXILIARY ambient space.\n"
  return closed_subscheme_of_toric_variety(ambient_space(m), [hypersurface_equation(m)])
end

@doc raw"""
    toric_hypersurface(m::AbstractFTheoryModel;
                      check_smoothness=false, check_codimension=true)

Adapt a concrete toric hypersurface, Weierstrass, or global Tate model to a
`ToricHypersurface`, retaining its actual equation. The underlying closed
subscheme is `calabi_yau_hypersurface(m)`. Expected codimension and Cartier
classes are checked by default; smoothness remains unknown unless requested.
Set `check_codimension=false` to assume the complete-intersection conditions.
Models over auxiliary families or nontoric ambient spaces are unsupported.
"""
function toric_hypersurface(
  m::AbstractFTheoryModel;
  check_smoothness::Bool=false, check_codimension::Bool=true,
)
  @req m isa Union{HypersurfaceModel,WeierstrassModel,GlobalTateModel} "A hypersurface model is required; for a CompleteIntersectionModel supply its chosen equations to toric_complete_intersection"
  @req base_space(m) isa NormalToricVariety "The model must have a concrete toric base"
  @req ambient_space(m) isa NormalToricVariety "The model must have a concrete toric ambient space"
  return toric_hypersurface(calabi_yau_hypersurface(m), hypersurface_equation(m); check_smoothness, check_codimension)
end

@doc raw"""
    defining_ideal(m::CompleteIntersectionModel)

Return the stored ideal of a complete-intersection model. Its generators are
not assumed to be a chosen regular sequence.
"""
defining_ideal(m::CompleteIntersectionModel) = m.defining_ideal

@doc raw"""
    toric_complete_intersection(m::CompleteIntersectionModel, equations::AbstractVector;
                               check_smoothness=false, check_codimension=true)

Adapt a concrete toric model using an explicitly supplied equation list that
generates its stored ideal. The model stores only an ideal, so a chosen
complete-intersection presentation cannot be recovered automatically.
The smoothness and codimension options are independent, as for the generic
toric constructor. This adapter does not change the model or its serialization.
"""
function toric_complete_intersection(
  m::CompleteIntersectionModel, equations::AbstractVector;
  check_smoothness::Bool=false, check_codimension::Bool=true,
)
  @req base_space(m) isa NormalToricVariety "The model must have a concrete toric base"
  @req ambient_space(m) isa NormalToricVariety "The model must have a concrete toric ambient space"
  c = closed_subscheme_of_toric_variety(ambient_space(m), defining_ideal(m))
  return toric_complete_intersection(c, equations; check_smoothness, check_codimension)
end

######################################################################################
### Attributes for flux families (not exported, rather for serialization overhaul)
######################################################################################

function _normalize_flux_family_algorithm(
  algorithm::Union{Symbol,String}, caller::Symbol
)
  normalized_algorithm = algorithm isa String ? Symbol(algorithm) : algorithm
  normalized_algorithm in (:default, :special) ||
    throw(ArgumentError("Unknown value for keyword `algorithm`: $algorithm"))
  if algorithm isa String
    Base.depwarn(
      "Passing `algorithm=\"$algorithm\"` is deprecated; use `algorithm=:$algorithm` instead",
      caller,
    )
  end
  return normalized_algorithm
end

@attr QQMatrix function matrix_integral_quant_transverse(
  m::AbstractFTheoryModel;
  completeness_check::Bool=true,
  algorithm::Union{Symbol,String}=:default,
  rng::AbstractRNG=Random.default_rng(),
)
  return matrix_integral(
    special_flux_family(m; completeness_check, algorithm, rng)
  )
end

@attr QQMatrix function matrix_rational_quant_transverse(
  m::AbstractFTheoryModel;
  completeness_check::Bool=true,
  algorithm::Union{Symbol,String}=:default,
  rng::AbstractRNG=Random.default_rng(),
)
  return matrix_rational(
    special_flux_family(m; completeness_check, algorithm, rng)
  )
end

@attr Vector{QQFieldElem} function offset_quant_transverse(
  m::AbstractFTheoryModel;
  completeness_check::Bool=true,
  algorithm::Union{Symbol,String}=:default,
  rng::AbstractRNG=Random.default_rng(),
)
  return offset(special_flux_family(m; completeness_check, algorithm, rng))
end

@attr QQMatrix function matrix_integral_quant_transverse_nobreak(
  m::AbstractFTheoryModel;
  completeness_check::Bool=true,
  algorithm::Union{Symbol,String}=:default,
  rng::AbstractRNG=Random.default_rng(),
)
  return matrix_integral(
    special_flux_family(
      m; not_breaking=true, completeness_check, algorithm, rng
    ),
  )
end

@attr QQMatrix function matrix_rational_quant_transverse_nobreak(
  m::AbstractFTheoryModel;
  completeness_check::Bool=true,
  algorithm::Union{Symbol,String}=:default,
  rng::AbstractRNG=Random.default_rng(),
)
  return matrix_rational(
    special_flux_family(
      m; not_breaking=true, completeness_check, algorithm, rng
    ),
  )
end

@attr Vector{QQFieldElem} function offset_quant_transverse_nobreak(
  m::AbstractFTheoryModel;
  completeness_check::Bool=true,
  algorithm::Union{Symbol,String}=:default,
  rng::AbstractRNG=Random.default_rng(),
)
  return offset(
    special_flux_family(
      m; not_breaking=true, completeness_check, algorithm, rng
    ),
  )
end
