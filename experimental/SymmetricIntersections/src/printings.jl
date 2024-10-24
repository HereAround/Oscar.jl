###############################################################################
#
#  Elevators
#
###############################################################################

function Base.show(io::IO, EC::ElevCtx{T, U}) where {T, U}
  d = degree_of_elevations(EC)
  if get(io, :supercompact, false)
    print(io, "Elevator")
  else
    print(io, "Degree $d elevator of a list with entries of type $T")
  end
end

###############################################################################
#
#  Representations
#
###############################################################################

# Linear representations

function Base.show(io::IO, ::MIME"text/plain", LR::LinRep)
  println(io, "Linear representation")
  io = AbstractAlgebra.pretty(io)
  println(io, AbstractAlgebra.Indent(), "of ", underlying_group(LR))
  println(io, "over ", AbstractAlgebra.Lowercase(), base_field(representation_ring(LR)))
  print(io, AbstractAlgebra.Dedent(), "of dimension ", dimension_representation(LR))
end 

function Base.show(io::IO, LR::LinRep)
  if get(io, :supercompact, false)
    print(io, "Linear representation")
  else
    print(io, "Linear representation of finite group of dimension $(dimension_representation(LR))")
  end
end

# Projective representations

function Base.show(io::IO, ::MIME"text/plain", PR::ProjRep)
  println(io, "Projective representation")
  io = AbstractAlgebra.pretty(io)
  println(io, AbstractAlgebra.Indent(), "of ", underlying_group(PR))
  println(io, "over ", AbstractAlgebra.Lowercase(), base_field(representation_ring_linear_lift(PR)))
  print(io, AbstractAlgebra.Dedent(), "of dimension ", dimension_representation(PR))
end

function Base.show(io::IO, PR::ProjRep)
  if get(io, :supercompact, false)
    print(io, "Projective representation")
  else
    print(io, "Projective representation of finite group of dimension ", dimension_representation(PR))
  end
end

# Representation rings

function Base.show(io::IO, ::MIME"text/plain", RR::RepRing)
  println(io, "Representation ring")
  io = AbstractAlgebra.pretty(io)
  println(io, AbstractAlgebra.Indent(), "of ", underlying_group(RR))
  print(io, "over ", AbstractAlgebra.Lowercase(), base_field(RR))
  print(io, AbstractAlgebra.Dedent())
end

function Base.show(io::IO, RR::RepRing)
  if get(io, :supercompact, false)
    print(io, "Representation ring")
  else
    print(io, "Representation ring of finite group over a field of characteristic 0")
  end
end

###############################################################################
#
#  Symmetric Grassmannians
#
###############################################################################

# Isotypical Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::IsotGrass)
  chi = submodule_character(M)
  println(io, "Symmetric Grassmannian of $(degree(chi))-dimensional submodules")
  io = AbstractAlgebra.pretty(io)
  println(io, AbstractAlgebra.Indent(), "of ", AbstractAlgebra.Lowercase(), module_representation(M))
  println(io, AbstractAlgebra.Dedent(), "with isotypical character")
  print(io, chi)
end

function Base.show(io::IO, M::IsotGrass)
  if get(io, :supercompact, false)
    print(io, "Isotypical Grassmannian")
  else
    print(io, "Isotypical Grassmannian of dimension $(projective_dimension(M))")
  end
end

# Character Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::CharGrass)
  chi = submodule_character(M)
  println(io, "Symmetric Grassmannian of $(degree(chi))-dimensional submodules")
  io = AbstractAlgebra.pretty(io)
  print(io, AbstractAlgebra.Indent(), "of ", AbstractAlgebra.Lowercase(), module_representation(M))
  println(io)
  println(io, AbstractAlgebra.Dedent(), "with character")
  print(io, chi)
end

function Base.show(io::IO, M::CharGrass)
  if get(io, :supercompact, false)
    print(io, "Character Grassmannian")
  else
    print(io, "Character Grassmannian of dimension $(projective_dimension(M))")
  end
end

# Determinant Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::DetGrass)
  chi = submodule_determinant_character(M)
  println(io, "Symmetric Grassmannian of $(M.d)-dimensional submodules")
  io = AbstractAlgebra.pretty(io)
  println(io, AbstractAlgebra.Indent(), "of ", AbstractAlgebra.Lowercase(), module_representation(M))
  println(io, AbstractAlgebra.Dedent(), "with determinant character")
  print(io, chi)
end

function Base.show(io::IO, M::DetGrass)
  if get(io, :supercompact, false)
    print(io, "Determinant Grassmannian")
  else
    print(io, "Determinant Grassmannian of dimension $(projective_dimension(M))")
  end
end

# Invariant Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::InvGrass)
  println(io, "Symmetric Grassmannian of $(M.d)-dimensional submodules")
  io = AbstractAlgebra.pretty(io)
  print(io, AbstractAlgebra.Indent(), "of ", AbstractAlgebra.Lowercase(), module_representation(M))
  print(io, AbstractAlgebra.Dedent())
end

function Base.show(io::IO, M::InvGrass)
  if get(io, :supercompact, false)
    print(io, "Invariant Grassmannian")
  else
    print(io, "Invariant Grassmannian of dimension $(projective_dimension(M))")
  end
end

###############################################################################
#
#  Symmetric intersections
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", symci::SymInter)
  M = symci.para
  j = symci.j
  RR = representation_ring_linear_lift(symci.prep)
  F = base_field(RR)
  G = underlying_group(symci.prep)
  n = ngens(codomain(j))
  t = submodule_dimension(M)
  d = total_degree(j(gens(domain(j))[1]))
  if t == 1
    ty = "($d)"
  elseif t == 2
    ty = "($d, $d)"
  else
    ty = "($d, "
    for i in 2:t-1
      ty *= "$d, "
    end
    ty *= "$d)"
  end
  io = AbstractAlgebra.pretty(io)
  println(io, "Parameter space for intersections")
  println(io, AbstractAlgebra.Indent(), "of type $(ty)")
  println(io, "in projective $(n-1)-space")
  print(io, AbstractAlgebra.Indent(), "over ", AbstractAlgebra.Lowercase())
  Base.show(IOContext(io, :supercompact => true), F)
  println(io)
  print(io, AbstractAlgebra.Dedent(), AbstractAlgebra.Dedent(), "preserved under the action of ", AbstractAlgebra.Lowercase(), G)
end

function Base.show(io::IO, symci::SymInter)
  if get(io, :supercompact, false)
    print(io, "Symmetric intersections")
  else
    print(io, "Parameter space for symmetric intersections")
  end
end

