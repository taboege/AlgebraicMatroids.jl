push!(LOAD_PATH, ".")
using AlgebraicMatroids, Oscar
using Combinatorics: powerset

R, A,B,C,D = polynomial_ring(QQ, "A" => 1:2, "B" => 1:2, "C" => 1:2, "D" => 1:3)
I = ideal(R, [
  C[1]*A[1] + C[2] - A[2],
  C[1]*B[1] + C[2] - B[2],
  D[1]*A[1]^2 + D[2]*A[1] + D[3] - A[2],
  D[1]*B[1]^2 + D[2]*B[1] + D[3] - B[2],
])

abscomp = absolute_primary_decomposition(I)
P = abscomp[1][1]
M = algebraic_matroid(P)
for K in powerset([A, B, C, D])
  @show K, rank(M, Set(Iterators.flatten(K)))
end
