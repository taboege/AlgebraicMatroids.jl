module AlgebraicMatroids

export algebraic_matroid, ground_set, prime_ideal, base_ring, coefficient_ring,
  coordinate_ideal, rank, is_independent, is_basis, is_circuit, fundamental_circuit,
  circuit_polynomial, base_degree

using Oscar
using Combinatorics: powerset

struct AlgebraicMatroid
  ground_set::Set
  ideal::Dict{Set, MPolyIdeal}
  rank::Dict{Set, Int}
end

"""
  algebraic_matroid(P::MPolyIdeal)

Constructs the algebraic matroid of the prime ideal `P`.
"""
function algebraic_matroid(P::MPolyIdeal)
  is_prime(P) || throw(DomainError(P, "ideal of an algebraic matroid must be prime"))

  E = Set(gens(base_ring(P)))
  AlgebraicMatroid(E,
    Dict{Set, MPolyIdeal}([(E, P)]),
    Dict{Set, Int}(),
  )
end

"""
  ground_set(M)

Returns the ground set of the matroid `M`.
"""
function ground_set(M::AlgebraicMatroid)
  M.ground_set
end

"""
  prime_ideal(M)

Returns the prime ideal associated to the algebraic matroid `M`.
"""
function Oscar.prime_ideal(M::AlgebraicMatroid)
  M.ideal[M.ground_set]
end

"""
  polynomial_ring(M)

Returns the polynomial ring over which the algebraic matroid `M` is defined.
"""
function Oscar.polynomial_ring(M::AlgebraicMatroid)
  base_ring(prime_ideal(M))
end

Oscar.base_ring(M::AlgebraicMatroid) = polynomial_ring(M)

"""
  coefficient_ring(M)

Returns the coefficient ring of `base_ring(M)`, e.g., `QQ`.
"""
function coefficient_ring(M::AlgebraicMatroid)
  base_ring(polynomial_ring(M))
end

"""
  coordinate_ideal(M, S)

Returns the elimination ideal of the algebraic matroid `M` in which only
coordinates in `S` are present.
"""
function coordinate_ideal(M::AlgebraicMatroid, S::Set)
  haskey(M.ideal, S) && return M.ideal[S]
  M.ideal[S] = eliminate(prime_ideal(M), collect(setdiff(ground_set(M), S)))
end

Oscar.ideal(M::AlgebraicMatroid) = prime_ideal(M)
Oscar.ideal(M::AlgebraicMatroid, S::Set) = coordinate_ideal(M, S)

"""
  rank(M, S)

Returns the rank of `S` in `M`.
"""
function Oscar.rank(M::AlgebraicMatroid, S::Set)
  length(S) - codim(coordinate_ideal(M, S))
end

"""
  rank(M)

Returns the rank of the matroid `M`, which is the rank of its `ground_set`.
"""
Oscar.rank(M::AlgebraicMatroid) = rank(M, ground_set(M))

"""
  is_independent(M, S)

Checks if `S` is an independent set in `M`.
"""
function is_independent(M::AlgebraicMatroid, S::Set)
  rank(M, S) == length(S)
end

"""
  is_basis(M, S)

Checks if `S` is a basis in `M`.
"""
function is_basis(M::AlgebraicMatroid, S::Set)
  rank(M) == length(S) && is_independent(M, S)
end

"""
  is_circuit(M, S)

Checks if `S` is a circuit in `M`. The implementation rests on the
fact that `S` is a circuit in an algebraic matroid `M` if and only
if `coordinate_ideal(M, S)` is principal with generator `f` and all
variables in `S` occur in `f`.
"""
function is_circuit(M::AlgebraicMatroid, S::Set)
  I = coordinate_ideal(M, S)
  length(gens(I)) == 1 && Set(vars(gens(I)[1])) == S
end

"""
  fundamental_circuit(M, B, x)

Compute the fundamental circuit of element `x` with respect to basis `B`
in the matroid `M`.
"""
function fundamental_circuit(M::AlgebraicMatroid, B::Set, x)
  first([S for S in powerset([x, B...]) if is_circuit(M, S)])
end

"""
  circuit_polynomial(M, C)

Compute the circuit polynomial of the circuit `C` in the matroid `M`.
This is the up to scale unique generator of the principal ideal
"""
function circuit_polynomial(M::AlgebraicMatroid, C::Set)
  @assert is_circuit(M, C)
  I = coordinate_ideal(M, C)
  gens(I)[1]
end

"""
  circuit_polynomial(M, B, x)

Compute the `circuit_polynomial` of the `fundamental_circuit(M, B, x)`.
"""
circuit_polynomial(M::AlgebraicMatroid, B::Set, x) =
  circuit_polynomial(M, fundamental_circuit(M, B, x))

"""
  base_degree(M, B)
"""
function base_degree(M::AlgebraicMatroid, B::Set)
  @assert is_basis(M, B)
  a, b = BigInt(1), BigInt(2)^100
  K = coefficient_ring(M)
  I = ideal(polynomial_ring(M), [x - K(rand(a:b)) for x in B])
  degree(I + prime_ideal(M))
end

end
