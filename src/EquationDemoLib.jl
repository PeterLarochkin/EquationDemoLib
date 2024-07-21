module EquationDemoLib
using Printf
export q, u, F, psi, apply_A, get_B, scalar_product, proximity_search

"""
    q(x::Float64, y::Float64)::Float64

This function calculates the sum of two floating-point numbers `x` and `y`.

# Arguments
- `x::Float64`: The first number.
- `y::Float64`: The second number.

# Returns
- `Float64`: The sum of `x` and `y`.
"""
function q(x::Float64,y::Float64)::Float64
    return x + y
end 

"""
    u(x::Float64, y::Float64) -> Float64

Calculates the analytical solution for a given experiment based on the input values `x` and `y`. 

The function computes the square root of the sum of 4.0 and the product of `x` and `y`.

# Arguments
- `x::Float64`: A floating-point number representing the first input parameter.
- `y::Float64`: A floating-point number representing the second input parameter.

# Returns
- `Float64`: The result of the analytical calculation, which is the square root of the sum of 4.0 and `x*y`.
"""
function u(x::Float64,y::Float64)::Float64
    return sqrt(4.0+x*y)
end 

"""
    F(x::Float64, y::Float64) -> Float64

Computes the right part of the equation for the experiment. The function calculates a specific value based on the input parameters `x` and `y`, which involves an analytical expression.

The function performs the following steps:
1. Computes the intermediate value `u` as the square root of `4.0 + x * y`.
2. Uses `u` to compute the final result, which is the sum of two terms:
   - The first term is \$(\\frac{x^2 + y^2}{4u^3})\$.
   - The second term is \$((x + y) \\cdot u)\$.

# Arguments
- `x::Float64`: A floating-point number representing the first input parameter.
- `y::Float64`: A floating-point number representing the second input parameter.

# Returns
- `Float64`: The computed result of the right part of the equation, which combines the contributions from `x` and `y` and the intermediate value `u`.
"""
function F(x::Float64,y::Float64)::Float64
    u::Float64 = sqrt(4.0+x*y)
    return 1/(4*u*u*u)*(x*x + y*y) + (x + y)*u
end 
"""
    psi(x::Float64, y::Float64, A1::Float64, A2::Float64, B1::Float64, B2::Float64, h1::Float64, h2::Float64) -> Float64

Defines the boundary values for a specific rectangular area. The function calculates the value at the boundary of the rectangular domain based on the input parameters. 

The function handles various edge cases of the rectangle defined by corners (A1, B1), (A2, B1), (A2, B2), and (A1, B2), and computes values based on different conditions:

1. **Vertical edges:**
   - For (x = A2 and (B1 < y < B2, it uses the function `u(x, y)` to compute the boundary value.
   - For (x = A1 and (B1 < y < B2, it returns -y/4 + 2.

2. **Horizontal edges:**
   - For y = B2 and A1 < x < A2, it uses the function `u(x, y)` to compute the boundary value.
   - For y = B1 and A1 < x < A2, it returns -x/4 + 2.

3. **Corners:**
   - For (x, y) = (A1, B1), it returns a weighted average of boundary values at (A1, B1).
   - For (x, y) = (A1, B2), it returns a weighted average of boundary values at (A1, B2).
   - For (x, y) = (A2, B1), it returns a weighted average of boundary values at (A2, B1).
   - For (x, y) = (A2, B2), it returns a weighted average of boundary values at (A2, B2).

# Arguments
- `x::Float64`: The x-coordinate of the point.
- `y::Float64`: The y-coordinate of the point.
- `A1::Float64`: The x-coordinate of the left boundary of the rectangle.
- `A2::Float64`: The x-coordinate of the right boundary of the rectangle.
- `B1::Float64`: The y-coordinate of the bottom boundary of the rectangle.
- `B2::Float64`: The y-coordinate of the top boundary of the rectangle.
- `h1::Float64`: A weighting factor for the boundary values.
- `h2::Float64`: A weighting factor for the boundary values.

# Returns
- `Float64`: The computed boundary value at the point (x, y).

# Throws
- Throws an `ErrorException` if the point (x, y) is not on the boundary of the rectangle defined by [A1, A2] \times [B1, B2].
"""
function psi(x::Float64, y::Float64, A1::Float64, A2::Float64, B1::Float64, B2::Float64, h1::Float64, h2::Float64)::Float64
    if x==A2 && B1<y && y<B2
        u_ = u(x,y)
        return y/(2*u_)+u_
    elseif x==A1 && B1<y && y<B2
        return -y/4+2
    elseif y==B2 && A1<x && x<A2
        u_ = u(x,y)
        return x/(2*u_)+u_
    elseif y==B1 && A1<x && x<A2
        return -x/4+2
    elseif x==A1 && y==B1
        return (h1*(-x/4+2) + h2*(-y/4+2))/(h1+h2)
    elseif x==A1 && y==B2
        u_ = u(x,y)
        return (h1*(x/(2*u_)+u_)+h2*(-y/4+2))/(h1+h2)
    elseif x==A2 && y==B1
        u_ = u(x,y)
        return (h1*(-x/4+2) + h2*(y/(2*u_)+u_)) / (h1 + h2)
    elseif x==A2 && y==B2
        u_ = u(x,y)
        return (h1*(x/(2*u_)+u_)+h2*(y/(2*u_)+u_))/(h1+h2)
    end
    error("($x, $y) is not in square [$A1,$A2]x[$B1,$B2]")
end
"""
    apply_A(matrix_apply_to::Matrix{Float64}, A1::Float64, A2::Float64, B1::Float64, B2::Float64, h1::Float64, h2::Float64) -> Matrix{Float64}

Applies a specific operator `A` to a given matrix. The operator is defined for the boundary conditions of a rectangular domain. This function returns a new matrix that is the result of applying the operator `A` to the input matrix `matrix_apply_to`, based on the specified boundaries and weighting factors.

The function iterates over each element of the input matrix and applies the operator `A` based on the position of the element relative to the boundaries of the rectangular domain.

# Arguments
- `matrix_apply_to::Matrix{Float64}`: The matrix to which the operator `A` will be applied. The function returns a new matrix with the results.
- `A1::Float64`: The x-coordinate of the left boundary of the rectangular domain.
- `A2::Float64`: The x-coordinate of the right boundary of the rectangular domain.
- `B1::Float64`: The y-coordinate of the bottom boundary of the rectangular domain.
- `B2::Float64`: The y-coordinate of the top boundary of the rectangular domain.
- `h1::Float64`: A weighting factor used in the operator `A`.
- `h2::Float64`: Another weighting factor used in the operator `A`.

# Returns
- `Matrix{Float64}`: A new matrix with the operator `A` applied, based on the boundary conditions.

"""
function apply_A(matrix_apply_to::Matrix{Float64}, A1::Float64, A2::Float64, B1::Float64, B2::Float64, h1::Float64, h2::Float64)::Matrix{Float64}
    (M, N) = size(matrix_apply_to, 1), size(matrix_apply_to, 2)
    result_matrix = zeros(Float64, M, N)
    # with padding, inside "picture"
    for i in 2:M-1
        for j in 2:N-1
            # here is (7) equation works
            result_matrix[i, j] =   (matrix_apply_to[i,j] * (2/(h1*h1) + 2/(h2*h2) + q(A1+ (i-1)*h1, B1+ (j-1)*h2)) + 
                                    matrix_apply_to[i-1,j] * (-1/(h1*h1)) +
                                    matrix_apply_to[i+1,j] * (-1/(h1*h1)) +
                                    matrix_apply_to[i,j-1] * (-1/(h2*h2)) +
                                    matrix_apply_to[i,j+1] * (-1/(h2*h2)))
        end
    end
    for i in 2:M-1
        # it's (10) equations
        # i=1,M-1
        # top applying
        result_matrix[i, N] = 2/(h2*h2) * (matrix_apply_to[i, N] - matrix_apply_to[i, N-1]) +
                            ( q(A1+ (i-1)*h1,B1+ (N-1)*h2) + 2/h2 ) * matrix_apply_to[i, N] -
                            1/(h1*h1)*(matrix_apply_to[i+1, N] - matrix_apply_to[i, N] - matrix_apply_to[i, N]+ matrix_apply_to[i-1, N])
        # bottom applying
        result_matrix[i, 1] = -2/(h2*h2) * (matrix_apply_to[i, 2] - matrix_apply_to[i, 1]) +
                            ( q(A1+ (i-1)*h1,B1+ 0*h2) + 2/h2 ) * matrix_apply_to[i, 1] -
                            1/(h1*h1)*(matrix_apply_to[i+1, 1] - matrix_apply_to[i, 1] - matrix_apply_to[i, 1]+ matrix_apply_to[i-1, 1])
    end

    for j in 2:N-1
        # it's (9) equations
        # j=1,N-1
        # right applying
        result_matrix[M, j] = 2/(h1*h1) * (matrix_apply_to[M, j] - matrix_apply_to[M-1, j]) + 
                            (q(A1+ (M-1)*h1,B1+ (j-1)*h2) + 2/h1) * matrix_apply_to[M, j] - 
                            1/(h2*h2)*(matrix_apply_to[M, j+1] - matrix_apply_to[M, j] - matrix_apply_to[M, j]+ matrix_apply_to[M, j-1])
        # left applying
        result_matrix[1, j] = -2/(h1*h1) * (matrix_apply_to[2, j] - matrix_apply_to[1, j]) + 
                            (q(A1+ 0*h1,B1+ (j-1)*h2) + 2/h1) * matrix_apply_to[1, j] - 
                            1/(h2*h2)*(matrix_apply_to[1, j+1] - matrix_apply_to[1, j] - matrix_apply_to[1, j]+ matrix_apply_to[1, j-1]);
    end

    # remaining corner points
    # bottom left
    # it's (11) equation
    result_matrix[1, 1] = -2/(h1*h1)*(matrix_apply_to[2, 1] - matrix_apply_to[1, 1]) - 
                         2/(h2*h2)*(matrix_apply_to[1, 2] - matrix_apply_to[1, 1]) +
                         (q(A1+ 0*h1,B1+ 0*h2) + 2/h1 + 2/h2) * matrix_apply_to[1, 1]

    # it's (12) equation
    # bottom right
    result_matrix[M, 1] = 2/(h1*h1)*(matrix_apply_to[M, 1] - matrix_apply_to[M-1, 1]) - 
                        2/(h2*h2)*(matrix_apply_to[M, 2] - matrix_apply_to[M, 1]) +
                        (q(A1+ (M-1)*h1,B1+ 0*h2) + 2/h1 + 2/h2) * matrix_apply_to[M, 1];
    # it's (13) equation
    # top right
    result_matrix[M, N] = 2/(h1*h1)*(matrix_apply_to[M, N] - matrix_apply_to[M-1, N]) +
                        2/(h2*h2)*(matrix_apply_to[M, N] - matrix_apply_to[M, N-1]) +
                        (q(A1+ (M-1)*h1,B1+ (N-1)*h2) + 2/h1 + 2/h2) * matrix_apply_to[M, N];
    # it's (14) equation
    # top left
    result_matrix[1, N] = -2/(h1*h1)*(matrix_apply_to[2, N]- matrix_apply_to[1, N])+
                         2/(h2*h2)*(matrix_apply_to[1, N]- matrix_apply_to[1, N-1])+
                         (q(A1+ 0*h1,B1+ (N-1)*h2) + 2/h1 + 2/h2) * matrix_apply_to[1, N];
    return result_matrix
end
"""
    get_B(M::Int, N::Int, h1::Float64, h2::Float64, A1::Float64, A2::Float64, B1::Float64, B2::Float64) -> Matrix{Float64}

Generates a matrix `B` based on the boundary conditions of a rectangular domain. The matrix `B` is computed by applying either the `psi` function (for boundary edges) or the `F` function (for interior points) based on the location of each matrix element relative to the defined boundaries.

The function initializes a matrix of size `M` by `N`, where each element is computed as follows:
- **Edges**: For elements on the boundary of the rectangular domain, the function uses the `psi` function to determine the value.
- **Interior**: For elements inside the boundary, the function uses the `F` function.

# Arguments
- `M::Int`: Number of rows in the resulting matrix.
- `N::Int`: Number of columns in the resulting matrix.
- `h1::Float64`: Weighting factor for the boundary conditions.
- `h2::Float64`: Another weighting factor for the boundary conditions.
- `A1::Float64`: The x-coordinate of the left boundary of the rectangular domain.
- `A2::Float64`: The x-coordinate of the right boundary of the rectangular domain.
- `B1::Float64`: The y-coordinate of the bottom boundary of the rectangular domain.
- `B2::Float64`: The y-coordinate of the top boundary of the rectangular domain.

# Returns
- `Matrix{Float64}`: A matrix of size `M` by `N` with values computed using the `psi` or `F` functions based on boundary conditions.
"""
function get_B(M::Int, N::Int, h1::Float64, h2::Float64, A1::Float64, A2::Float64, B1::Float64, B2::Float64)::Matrix{Float64}
    result_matrix = zeros(Float64, M, N)
    for i in 2:M-1
        for j in 2:N-1
            result_matrix[i,j] = F(A1+ (i-1)*h1, B1+ (j-1)*h2)
        end
    end
    for i in 2:M-1
        # it's (10) equations
        # i=1,M-1
        # top applying

        result_matrix[i, N] = psi(A1+ (i-1)*h1, B1+ (N-1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F(A1 + (i-1)*h1, B1 + (N-1)*h2)
        # bottom applying
        result_matrix[i, 1] = psi(A1+ (i-1)*h1, B1+ 0*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F(A1 + (i-1)*h1, B1 + 0*h2)
    end

    for j in 2:N-1
        # it's (9) equations
        # j=1,N-1
        # right applying
        result_matrix[M, j] = psi(A1+ (M-1)*h1, B1+ (j-1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F(A1 + (M-1)*h1, B1 + (j-1)*h2)
        # left applying
        result_matrix[1, j] = psi(A1+ 0*h1, B1+ (j-1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F(A1 + 0*h1, B1 + (j-1)*h2)
    end

    # remaining corner points
    # bottom left
    # it's (11) equation
    result_matrix[1,1] = psi(A1+ 0*h1, B1+ 0*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + 0*h2)
    # it's (12) equation
    # bottom right
    result_matrix[M,1] = psi(A1+ (M-1)*h1, B1+ 0*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + (M-1)*h1, B1 + 0*h2)
    # it's (13) equation
    # top right
    result_matrix[M, N] = psi(A1+ (M-1)*h1, B1+ (N-1)*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + (M-1)*h1, B1 + (N-1)*h2)
    # it's (14) equation
    # top left
    result_matrix[1, N] = psi(A1+ 0*h1, B1+ (N-1)*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + (N-1)*h2);
    return result_matrix
end

function ρ(index::Int, M_or_N::Int)::Float64
    if (index == 1 || index == M_or_N)
        return 0.5
    else
        return 1.0
    end
end
"""
   scalar_product(A::Matrix{Float64}, B::Matrix{Float64}, h1::Float64, h2::Float64) -> Float64

Computes the scalar product of two matrices `A` and `B` by performing element-wise multiplication and then summing the results. The function uses weighting factors `h1` and `h2` to adjust the product based on the index of each element.

The scalar product is computed as follows:
1. For each element `(i, j)` in the matrices `A` and `B`, compute the product `A[i, j] * B[i, j]`.
2. Apply a weight based on the indices `(i, j)` using the parameters `h1` and `h2`.
3. Sum all the weighted products to obtain the final scalar result.

# Arguments
- `A::Matrix{Float64}`: The first matrix for the scalar product calculation.
- `B::Matrix{Float64}`: The second matrix for the scalar product calculation.
- `h1::Float64`: A weighting factor used in the computation.
- `h2::Float64`: Another weighting factor used in the computation.

# Returns
- `Float64`: The scalar result of the element-wise product of `A` and `B`, weighted by `h1` and `h2`, and summed.
   scalar_product(A::Matrix{Float64}, B::Matrix{Float64}, h1::Float64, h2::Float64) -> Float64

Computes the scalar product of two matrices `A` and `B` by performing element-wise multiplication and then summing the results. The function uses weighting factors `h1` and `h2` to adjust the product based on the index of each element.

The scalar product is computed as follows:
1. For each element `(i, j)` in the matrices `A` and `B`, compute the product `A[i, j] * B[i, j]`.
2. Apply a weight based on the indices `(i, j)` using the parameters `h1` and `h2`.
3. Sum all the weighted products to obtain the final scalar result.

# Arguments
- `A::Matrix{Float64}`: The first matrix for the scalar product calculation.
- `B::Matrix{Float64}`: The second matrix for the scalar product calculation.
- `h1::Float64`: A weighting factor used in the computation.
- `h2::Float64`: Another weighting factor used in the computation.

# Returns
- `Float64`: The scalar result of the element-wise product of `A` and `B`, weighted by `h1` and `h2`, and summed.
"""
function scalar_product(A::Matrix{Float64}, B::Matrix{Float64}, h1::Float64, h2::Float64)::Float64
    sum = 0.0
    @assert size(A) == size(B)
    (M, N) = size(A, 1), size(A, 2)
    result_matrix = zeros(Float64, M, N)
    for i in 1:M
        for j in 1:N
            sum += h1*h2*ρ(i, M)*ρ(j, N)*A[i,j] * B[i,j];
        end
    end
    return sum

end

function multiply_by_num(A::Matrix{Float64}, num::Float64)::Matrix{Float64}
    return A*num
end

"""
    proximity_search(M::Int, N::Int)

Performs an iterative search to solve a given problem using an iterative method to update the solution based on residuals. This function iteratively refines the approximation `omega` to minimize the residual between `A * omega` and `B` until the change between successive iterations is below a specified tolerance.

The function performs the following steps:
1. Initializes the parameters, matrices, and iteration variables.
2. Computes the initial matrix `B` using boundary conditions.
3. Iteratively updates the solution `omega` using the following formula:
   \$\$ \\omega_{ij}^{(k+1)} = \\omega_{ij}^{(k)} - \\tau_{k+1} r_{ij}^{(k)} \$\$
   where the residual \$r^{(k)}\$ is calculated as \$r^{(k)} = A \\omega^{(k)} - B\$ and the iteration parameter \$\\tau_{k+1}\$ is computed as:
   \$\$ \\tau_{k+1} = \\frac{[A r^{(k)}, r^{(k)}]}{||A r^{(k)}||_{E}^{2}} \$\$
4. Stops the iteration when the change between successive iterations is less than the specified tolerance \$\\varepsilon\$.
5. Outputs the number of iterations and the maximum difference between the computed solution and the analytical solution.

# Arguments
- `M::Int`: Number of rows in the matrix for the problem.
- `N::Int`: Number of columns in the matrix for the problem.

# Returns
- None: The function prints the iteration results and maximum difference to the console.

"""
function proximity_search(M::Int, N::Int)
    epsilon::Float64 = 0.000001
    A1::Float64 = 0.0
    A2::Float64 = 4.0
    B1::Float64 = 0.0
    B2::Float64 = 3.0
    h1::Float64 = (A2 - A1)/M
    h2::Float64 = (B2 - B1)/N
    tau::Float64 = 0.0
    sq_eps::Float64 = epsilon * epsilon
    squared_difference::Float64 = sq_eps
    omega = zeros(Float64, M+1, N+1)
    omega_next = zeros(Float64, M+1, N+1)
    B = zeros(Float64, M+1, N+1)
    A_omega = zeros(Float64, M+1, N+1)
    r = zeros(Float64, M+1, N+1)
    A_r = zeros(Float64, M+1, N+1)
    tau_r = zeros(Float64, M+1, N+1)
    difference_omega = zeros(Float64, M+1, N+1)
    B = get_B(M+1, N+1, h1, h2, A1, A2, B1, B2);
    count = 0
    while (squared_difference >= sq_eps)
        omega = copy(omega_next)
        A_omega = apply_A(omega, A1, A2, B1, B2, h1, h2)
        r = A_omega - B
        A_r = apply_A(r, A1, A2, B1, B2, h1, h2)
        tau = scalar_product(A_r, r, h1, h2) / scalar_product(A_r, A_r, h1, h2)
        tau_r = multiply_by_num(r, tau)
        omega_next = omega - tau_r
        squared_difference = scalar_product(tau_r, tau_r, h1, h2)
        if (count % 1000 == 0)
            @printf("n:%d, diff:%.10f\n", count,sqrt(squared_difference))
        end
        count += 1
    end
    @printf("total count: %d\n", count)
    max_::Float64 = 0.0
    for i in 1:M+1
        for j in 1:N+1
            item = abs(omega_next[i,j]-u(h1*(i-1),h2*(j-1)))
            max_ = item > max_ ? item : max_
        end
    end
    @printf("diff:%.10f\n", max_)
end
end
