# Check Cox presentations, constructor invariants, and legacy compatibility.
@testset "Closed subschemes" begin
  antv = normal_toric_variety(
    incidence_matrix([[1, 2, 3, 4]]), [[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0]]
  )

  ntv = normal_toric_variety(cube(2))
  (x1, x2, y1, y2) = gens(cox_ring(ntv))
  sv1 = closed_subscheme_of_toric_variety(ntv, [x1])
  sv2 = closed_subscheme_of_toric_variety(ntv, [x1^2 + x1 * x2 + x2^2, y2])

  P3 = projective_space(NormalToricVariety, 3)
  sv3 = closed_subscheme_of_toric_variety(P3, [gens(cox_ring(P3))[1]^2])

  @testset "Invalid closed subscheme presentations" begin
    @test_throws ArgumentError closed_subscheme_of_toric_variety(ntv, [x1 - y1])
    @test_throws ArgumentError closed_subscheme_of_toric_variety(ntv, ideal([x1 - y1]))
    @test_throws ArgumentError ClosedSubschemeOfToricVariety(ntv, ideal([x1 - y1]))
    @test_throws ArgumentError closed_subscheme_of_toric_variety(ntv, defining_ideal(sv3))
    @test_throws ArgumentError closed_subscheme_of_toric_variety(
      antv, [gens(cox_ring(antv))[1]]
    )
  end

  @testset "Basic properties" begin
    @test is_empty(sv1) == false
    @test is_empty(sv2) == false
    @test is_empty(sv3) == false
    @test radical(sv1) == defining_ideal(sv1)
    @test dim(toric_variety(sv1)) == 2
  end

  @test getfield(Oscar, :ClosedSubvarietyOfToricVariety) === ClosedSubschemeOfToricVariety
  legacy = @test_deprecated closed_subvariety_of_toric_variety(ntv, defining_ideal(sv1))
  @test legacy isa ClosedSubschemeOfToricVariety
  @test defining_ideal(legacy) === defining_ideal(sv1)
end

# Test intrinsic properties on P¹, including Cox nilpotents and presentation-independent scheme cycles.
@testset "Associated toric schemes" begin
  P1 = projective_space(NormalToricVariety, 1)
  x, y = gens(cox_ring(P1))
  I = ideal([x^2, x*y])
  c = closed_subscheme_of_toric_variety(P1, I)
  @test defining_ideal(c) === I
  @test radical(c) != I
  @test radical(c) == ideal([x])
  @test dim(c) == 0
  @test is_reduced(c)
  @test is_irreducible(c)
  @test is_integral(c)
  @test is_smooth(c)

  A = closed_subscheme_of_toric_variety(P1, ideal([x]))
  B = closed_subscheme_of_toric_variety(P1, ideal([x^2, x * y]))
  IA = ideal_sheaf(P1, defining_ideal(A))
  IB = ideal_sheaf(P1, defining_ideal(B))
  @test covered_scheme(A) !== covered_scheme(B)
  @test IA == IB
  @test algebraic_cycle(IA) == algebraic_cycle(IB)
  @test_throws MethodError rational_equivalence_class(A)
  @test_throws MethodError rational_equivalence_class(B)

  doubled = toric_hypersurface(P1, x^2)
  @test has_expected_codimension(doubled)
  @test !has_attribute(doubled, :is_smooth)
  @test !is_reduced(doubled)
  @test is_irreducible(doubled)
  @test !is_integral(doubled)
  @test !is_smooth(doubled)
  ci = toric_complete_intersection(closed_subscheme_of_toric_variety(doubled), [x^2])
  @test !is_reduced(ci)
  @test is_irreducible(ci)
  @test !is_integral(ci)
  @test covered_scheme(ci) === covered_scheme(doubled)
  @test (@test_deprecated closed_subvariety_of_toric_variety(ci)) === closed_subscheme_of_toric_variety(ci)

  two_points = closed_subscheme_of_toric_variety(P1, [x*y])
  @test is_reduced(two_points)
  @test !is_irreducible(two_points)
  @test !is_integral(two_points)
  @test is_smooth(two_points)

  empty = closed_subscheme_of_toric_variety(P1, [x, y])
  @test is_empty(empty)
  @test dim(empty) == -Oscar.inf
  @test is_reduced(empty)
  @test !is_irreducible(empty)
  @test !is_integral(empty)
  @test is_smooth(empty)
end

# Test chosen presentations, scheme conversion, independent checks, and adjunction.
@testset "Toric complete intersections" begin
  P2 = projective_space(NormalToricVariety, 2)
  x, y, z = gens(cox_ring(P2))
  smooth = closed_subscheme_of_toric_variety(P2, [x^2 - y*z])
  singular = closed_subscheme_of_toric_variety(P2, [x*y])
  doubled = closed_subscheme_of_toric_variety(P2, [x^2])
  smooth_scheme = covered_scheme(smooth)
  @test smooth_scheme isa CoveredScheme
  @test smooth_scheme === covered_scheme(smooth)
  @test !has_attribute(smooth, :is_smooth)
  @test !is_reduced(covered_scheme(doubled))
  @test is_smooth(smooth)
  @test !is_smooth(singular)
  @test !is_smooth(doubled)
  @test radical(toric_hypersurface(doubled, x^2)) == ideal([x])
  @test radical(toric_complete_intersection(P2, [x^2, y^2])) == ideal([x, y])

  H = toric_hypersurface(P2, x^3 + y^3 + z^3)
  @test H isa AbsToricCompleteIntersection
  @test !has_attribute(H, :is_smooth)
  @test !has_attribute(closed_subscheme_of_toric_variety(H), :is_smooth)
  @test get_attribute(H, :has_expected_codimension)
  @test toric_variety(H) === P2
  @test defining_equation(H) == x^3 + y^3 + z^3
  @test defining_ideal(H) == ideal([defining_equation(H)])
  @test covered_scheme(H) === covered_scheme(closed_subscheme_of_toric_variety(H))
  @test divisor_class(only(defining_divisor_classes(H))) == degree(defining_equation(H))
  @test euler_characteristic(H) == ZZ(0)
  @test euler_characteristic(H) isa ZZRingElem
  @test !has_attribute(H, :is_smooth)
  @test !has_attribute(closed_subscheme_of_toric_variety(H), :is_smooth)
  @test euler_characteristic(H; check_smoothness=true, check_codimension=true) == 0
  @test get_attribute(H, :is_smooth)
  @test has_expected_codimension(H)
  @test dim(H) == 1
  @test !is_empty(H)
  @test closed_subscheme_of_toric_variety(toric_hypersurface(smooth, x^2 - y*z)) === smooth
  @test_throws ArgumentError toric_hypersurface(smooth, x*y)
  @test_throws ArgumentError toric_hypersurface(P2, x*y; check_smoothness=true)
  @test_throws ArgumentError toric_hypersurface(P2, x^2; check_smoothness=true)
  @test_throws ArgumentError euler_characteristic(toric_hypersurface(singular, x*y))
  unchecked = toric_hypersurface(P2, x*y; check_codimension=false)
  @test !has_attribute(unchecked, :has_expected_codimension)
  @test !has_attribute(closed_subscheme_of_toric_variety(unchecked), :_toric_subscheme)
  @test euler_characteristic(unchecked) == 2
  @test !has_attribute(unchecked, :is_smooth)
  @test_throws ArgumentError euler_characteristic(unchecked; check_smoothness=true)
  @test_throws ArgumentError euler_characteristic(unchecked)
  @test_throws ArgumentError euler_characteristic(H; algorithm=:unknown)
  @test_throws ArgumentError toric_hypersurface(P2, x + y^2)
  @test_throws ArgumentError toric_hypersurface(P2, zero(cox_ring(P2)))
  @test_throws ArgumentError toric_hypersurface(P2, one(cox_ring(P2)))
  @test_throws ArgumentError toric_complete_intersection(P2, [])
  @test_throws ArgumentError toric_complete_intersection(P2, [x, y, z])

  P3 = projective_space(NormalToricVariety, 3)
  w = gens(cox_ring(P3))
  @test_throws ArgumentError toric_hypersurface(P2, w[1])
  quartic = toric_hypersurface(P3, sum(t^4 for t in w); check_codimension=false)
  @test euler_characteristic(quartic) == 24
  @test euler_characteristic(P3, defining_divisor_classes(quartic)) == 24
  @test_throws ArgumentError euler_characteristic(P2, defining_divisor_classes(quartic))
  @test_throws ArgumentError chern_class(quartic, 3)
  fs = [sum(t^2 for t in w), sum(i*w[i]^2 for i in eachindex(w))]
  C = toric_complete_intersection(P3, fs; check_smoothness=true, check_codimension=true)
  @test is_smooth(C)
  @test dim(C) == 1
  @test has_expected_codimension(C)
  @test covered_scheme(C) === covered_scheme(closed_subscheme_of_toric_variety(C))
  @test defining_equations(C) == fs
  fs[1] = w[1]
  @test defining_equations(C)[1] != fs[1]
  @test euler_characteristic(C) == 0
  @test iszero(polynomial(chern_class(C, 1)))
  redundant = toric_complete_intersection(P3, [w[1], w[1]]; check_smoothness=true, check_codimension=false)
  @test is_smooth(redundant)
  @test !has_expected_codimension(redundant)
  @test_throws ArgumentError euler_characteristic(redundant)
  @test_throws ArgumentError toric_complete_intersection(P3, [w[1], w[1]])
  unchecked_codim = toric_complete_intersection(P3, [w[1], w[1]]; check_codimension=false)
  @test !has_attribute(unchecked_codim, :has_expected_codimension)
  euler_characteristic(unchecked_codim)
  @test_throws ArgumentError euler_characteristic(unchecked_codim; check_codimension=true)
  @test_throws ArgumentError euler_characteristic(unchecked_codim)

  P4 = projective_space(NormalToricVariety, 4)
  quintic = toric_hypersurface(P4, sum(t^5 for t in gens(cox_ring(P4))); check_codimension=false)
  @test euler_characteristic(quintic) == -200
  P1 = projective_space(NormalToricVariety, 1)
  @test euler_characteristic(toric_hypersurface(P1, sum(t^3 for t in gens(cox_ring(P1))))) == 3
  product = P1 * P1
  a, b, c, d = gens(cox_ring(product))
  diagonal = toric_hypersurface(product, a*c + b*d; check_smoothness=true, check_codimension=true)
  @test euler_characteristic(diagonal) == 2
  @test euler_characteristic(product, [anticanonical_divisor_class(product)]) == 0
  @test_throws ArgumentError toric_complete_intersection(product, [a, b])
  empty = toric_complete_intersection(product, [a, b]; check_codimension=false)
  @test is_empty(empty)
  # Check the actual empty covered scheme, including its base and dimension.
  empty_scheme = covered_scheme(empty)
  @test empty_scheme === covered_scheme(closed_subscheme_of_toric_variety(empty))
  @test empty_scheme isa CoveredScheme
  @test base_ring(empty_scheme) === coefficient_ring(product)
  @test isempty(affine_charts(empty_scheme))
  @test is_empty(empty_scheme)
  @test dim(empty_scheme) == -Oscar.inf
  @test is_smooth(empty_scheme)
  @test !has_expected_codimension(empty)
  @test_throws ArgumentError euler_characteristic(empty)

  weighted = weighted_projective_space(NormalToricVariety, [1, 1, 2])
  t = gens(cox_ring(weighted))
  @test !is_smooth(weighted)
  @test_throws ArgumentError toric_hypersurface(weighted, t[1])
  W = toric_hypersurface(weighted, t[1]^4 + t[2]^4 + t[3]^2)
  @test has_expected_codimension(W)
  @test_throws ArgumentError euler_characteristic(W; completeness_check=false)
  @test_throws ArgumentError euler_characteristic(weighted, defining_divisor_classes(W))
end

# Check completeness independently of cached rings, Chern classes, and prior unchecked calls.
@testset "Toric Euler check forwarding" begin
  X = normal_toric_variety(
    incidence_matrix([[1, 2], [2, 3], [1, 3]]), [[1, 0], [0, 1], [-1, -1]]
  )
  @test !has_attribute(X, :is_complete)
  x = gens(cox_ring(X))
  Y = toric_hypersurface(X, sum(t^3 for t in x))
  @test euler_characteristic(Y; completeness_check=false) == 0
  @test !has_attribute(X, :is_complete)
  @test euler_characteristic(Y; completeness_check=true) == 0
  @test get_attribute(X, :is_complete)

  A = normal_toric_variety(incidence_matrix([[1, 2]]), [[1, 0], [0, 1]])
  Z = toric_hypersurface(A, gens(cox_ring(A))[1])
  @test euler_characteristic(Z; completeness_check=false) == 0
  @test has_attribute(A, :cohomology_ring)
  @test has_attribute(A, :chern_classes)
  @test_throws ArgumentError euler_characteristic(Z)
  @test_throws ArgumentError chern_class(Z, 0)
end
