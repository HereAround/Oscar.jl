using Random
our_rng = Random.Xoshiro(1234)

qsm_model = literature_model(;
  arxiv_id="1903.00009", model_parameters=Dict("k" => 283), rng=our_rng
)
as = ambient_space(qsm_model)
x1, x2, x3, x4, x5, x6, x7, v, e3, e2, u, e4, e1, w = gens(cohomology_ring(as))
X1 = cohomology_class(as, x1)
X2 = cohomology_class(as, x2)
X3 = cohomology_class(as, x3)
X4 = cohomology_class(as, x4)
X5 = cohomology_class(as, x5)
X6 = cohomology_class(as, x6)
X7 = cohomology_class(as, x7)
E3 = cohomology_class(as, e3)
E4 = cohomology_class(as, e4)
E1 = cohomology_class(as, e1) + E4
E2 = cohomology_class(as, e2) + E3
U = cohomology_class(as, u)
V = cohomology_class(as, v)
W = cohomology_class(as, w)
H = W + E1
cy = cohomology_class(toric_divisor_class(as, degree(hypersurface_equation(qsm_model))))
c2_B = cohomology_class(as, 4 * x4 * x6 + 8 * x5 * x6 + 14 * x6 * x7 + 5 * x7^2)
Kbar = cohomology_class(ambient_space(qsm_model), x1 + x2 + x3 + x4 + x5 + x6 + x7)
Kbar3 = kbar3(qsm_model)
fg_not_breaking = special_flux_family(
  qsm_model; not_breaking=true, completeness_check=false, rng=our_rng
)
g4_exp = flux_instance(fg_not_breaking, [3], [])

@testset "Test properties of the QSM" begin
  # Legacy QSM database attributes are canonicalized when the model is loaded.
  @test has_attribute(qsm_model, :kbar3)
  @test !has_attribute(qsm_model, :Kbar3)
  @test has_attribute(qsm_model, :has_quick_triangulation)
  @test !has_attribute(qsm_model, :triang_quick)
  @test H == U + E1 + E2 + E4 == V + E2 + E3
  @test chern_classes(qsm_model)[2] ==
    c2_B - 7 * E4^2 - E1 * Kbar - E2 * Kbar - E3 * Kbar - E4 * Kbar + 3 * H * Kbar
  @test euler_characteristic(qsm_model) ==
    integrate((3 * Kbar * (4 * c2_B + 5 * Kbar^2)) * V * cy)
  @test d3_tadpole_constraint(g4_exp) == 12 + 5//8 * Kbar3 - 45//(2 * Kbar3)
end

@testset "Advanced intersection theory and QSM-fluxes" begin
  qsm_model = literature_model(;
    arxiv_id="1903.00009", model_parameters=Dict("k" => 4), rng=our_rng
  )
  h22_converter_dict = converter_dict_h22_ambient(qsm_model; completeness_check=false)
  coh_ring = cohomology_ring(ambient_space(qsm_model); completeness_check=false)
  coh_ring_gens = gens(coh_ring)
  for (key, value) in h22_converter_dict
    obj1 = coh_ring_gens[key[1]] * coh_ring_gens[key[2]]
    obj2 = sum(
      value[k][1] * coh_ring_gens[value[k][2][1]] * coh_ring_gens[value[k][2][2]] for
      k in 1:length(value)
    )
    @test obj1 == obj2
  end
  qsm_g4_flux = qsm_flux(qsm_model)

  # Check that string-valued algorithm names are deprecated.
  @test_deprecated Oscar._normalize_flux_family_algorithm(
    "special", :special_flux_family
  ) ===
    :special

  h22_basis = gens_of_h22_hypersurface_indices(qsm_model; completeness_check=false)
  flux_poly_str = string(polynomial(cohomology_class(qsm_g4_flux)))
  ring = base_ring(parent(polynomial(cohomology_class(qsm_g4_flux))))
  flux_poly = Oscar.eval_poly(flux_poly_str, ring)
  coeffs = collect(coefficients(flux_poly))
  raw_exps = collect(exponents(flux_poly))
  massaged_exps = [
    length(pos) == 1 ? (pos[1], pos[1]) : Tuple(sort(pos)) for
    pos in (findall(!=(0), mon) for mon in raw_exps)
  ]
  flux_vector = fill(QQ(0), length(h22_basis))
  for (i, exp_pair) in enumerate(massaged_exps)
    idx = findfirst(==(exp_pair), h22_basis)
    flux_vector[idx] = coeffs[i]
  end
  flux_vector = transpose(matrix(QQ, [flux_vector]))
  fg = special_flux_family(
    qsm_model; not_breaking=true, completeness_check=false, algorithm=:special, rng=our_rng
  )
  @test ncols(matrix_integral(fg)) == 1
  @test nrows(matrix_integral(fg)) == nrows(matrix_rational(fg))
  @test unique(offset(fg)) == [0]
  M1 = matrix_integral(fg)
  M2 = matrix_rational(fg)
  large_M = hcat(M1, M2)
  @test rank(large_M) == minimum(size(large_M))
  solution = solve(large_M, flux_vector; side=:right)
  @test is_integer(solution[1])
  reconstructed_flux = flux_instance(fg, matrix(ZZ, [[solution[1]]]), solution[2:end, :])
  @test cohomology_class(qsm_g4_flux) == cohomology_class(reconstructed_flux)

  # Empty coefficient lists must retain the affine offset of a flux family.
  shifted_offset = copy(offset(fg))
  shifted_offset[1] = 1
  shifted_fg = family_of_g4_fluxes(
    qsm_model, M1, M2, shifted_offset; completeness_check=false
  )
  @test fg == fg
  @test hash(fg) == hash(fg)
  @test shifted_fg != fg
  empty_coeff_flux = flux_instance(
    shifted_fg, Int[], Rational{Int}[]; completeness_check=false, consistency_check=false
  )
  zero_matrix_flux = flux_instance(
    shifted_fg, zero_matrix(ZZ, ncols(M1), 1), zero_matrix(QQ, ncols(M2), 1);
    completeness_check=false, consistency_check=false,
  )
  @test cohomology_class(empty_coeff_flux) == cohomology_class(zero_matrix_flux)
end

# Compare intersection theory with independent legacy Hodge attributes of the local k=4 QSM.
# The guard prevents downloads when this testset is selected in isolation.
# Running this entire file already requires QSMDB in the setup above.
@testset "QSM k=4 Euler comparison" begin
  qsm_hash = Oscar.artifact_hash("QSMDB", joinpath(Oscar.oscardir, "Artifacts.toml"))
  if isnothing(qsm_hash) || !Oscar.artifact_exists(qsm_hash)
    @info "Skipping k=4 Euler testset: QSMDB is not installed locally"
    @test_skip false
  else
    m = literature_model(; arxiv_id="1903.00009", model_parameters=Dict("k" => 4), rng=Random.Xoshiro(1234))
    hs = [get_attribute(m, s) for s in (:hodge_h11, :hodge_h12, :hodge_h13, :hodge_h22)]
    @test hs == [31, 10, 34, 284]
    h11, h12, h13, h22 = hs
    expected = 4 + 2*h11 - 4*h12 + 2*h13 + h22
    @test expected == 378
    Y = toric_hypersurface(m; check_codimension=false)
    @test !has_attribute(Y, :is_smooth)
    @test euler_characteristic(Y; completeness_check=false) == expected
    @test euler_characteristic(m; completeness_check=false) == expected
    D = cohomology_class(only(defining_divisor_classes(Y)); completeness_check=false)
    @test integrate(chern_class(m, 4; completeness_check=false) * D; completeness_check=false) == expected
    @test verify_euler_characteristic_from_hodge_numbers(m; completeness_check=false)
    @test_throws ArgumentError verify_euler_characteristic_from_hodge_numbers(m; completeness_check=false, algorithm=:unknown)
    @test !has_attribute(Y, :is_smooth)
    @test !has_attribute(calabi_yau_hypersurface(m), :is_smooth)
  end
end
