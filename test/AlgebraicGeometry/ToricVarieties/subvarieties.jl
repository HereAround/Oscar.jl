@testset "Closed subvarieties" begin
  antv = normal_toric_variety(
    incidence_matrix([[1, 2, 3, 4]]), [[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0]]
  )

  ntv = normal_toric_variety(cube(2))
  (x1, x2, y1, y2) = gens(cox_ring(ntv))
  sv1 = closed_subvariety_of_toric_variety(ntv, [x1])
  sv2 = closed_subvariety_of_toric_variety(ntv, [x1^2 + x1 * x2 + x2^2, y2])

  P3 = projective_space(NormalToricVariety, 3)
  sv3 = closed_subvariety_of_toric_variety(P3, [gens(cox_ring(P3))[1]^2])

  @testset "Should fail due to bad arguments (closed subvarieties of toric varieties)" begin
    @test_throws ArgumentError closed_subvariety_of_toric_variety(ntv, [x1 - y1])
    @test_throws ArgumentError closed_subvariety_of_toric_variety(
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
end


# Check that chosen Cox equations yield certified geometric presentations.
@testset "Toric hypersurface and complete-intersection presentations" begin
  X = projective_space(NormalToricVariety, 2)
  x, y, z = gens(cox_ring(X))
  H = toric_hypersurface(X, x^2 - y*z)
  @test H isa ToricHypersurface
  @test has_expected_codimension(H)
  @test is_smooth(H)
  @test defining_equation(H) == x^2 - y*z
  @test defining_divisor_classes(H) == [toric_divisor_class(X, degree(x^2))]
  @test fundamental_class(H) ==
    rational_equivalence_class(only(defining_divisor_classes(H)))
  @test fundamental_class(H) == rational_equivalence_class(H)
  @test covered_scheme(H) === covered_scheme(H)

  C = toric_complete_intersection(X, [x, y])
  @test C isa ToricCompleteIntersection
  @test has_expected_codimension(C)
  @test dim(C) == 0
  @test_throws ArgumentError toric_complete_intersection(X, [x, x])
  @test_throws MethodError ToricCompleteIntersection(
    X, [x], [toric_divisor_class(X, degree(x^3))]
  )
  @test_throws ArgumentError toric_hypersurface(X, x + y^2)
  @test_throws ArgumentError toric_hypersurface(X, zero(cox_ring(X)))

  singular = toric_hypersurface(X, x*y)
  @test !is_smooth(singular)
  @test is_smooth(covered_scheme(H))
end
