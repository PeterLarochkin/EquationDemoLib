module EquationDemoLibParallel
using Printf
using Base.Threads
export q, u, F, psi, apply_A, get_B, scalar_product, proximity_search


function q(x::Float64,y::Float64)::Float64
    return x + y
end 

function u(x::Float64,y::Float64)::Float64
    return sqrt(4.0+x*y)
end 

function F(x::Float64,y::Float64)::Float64
    u::Float64 = sqrt(4.0+x*y)
    return 1/(4*u*u*u)*(x*x + y*y) + (x + y)*u
end 

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

function apply_A(matrix_apply_to::Matrix{Float64}, A1::Float64, A2::Float64, B1::Float64, B2::Float64, h1::Float64, h2::Float64)::Matrix{Float64}
    (M, N) = size(matrix_apply_to, 1), size(matrix_apply_to, 2)
    result_matrix = zeros(Float64, M, N)
    # with padding, inside "picture"
    @threads for i in 2:M-1
        @threads for j in 2:N-1
            # here is (7) equation works
            result_matrix[i, j] =   (matrix_apply_to[i,j] * (2/(h1*h1) + 2/(h2*h2) + q(A1+ (i-1)*h1, B1+ (j-1)*h2)) + 
                                    matrix_apply_to[i-1,j] * (-1/(h1*h1)) +
                                    matrix_apply_to[i+1,j] * (-1/(h1*h1)) +
                                    matrix_apply_to[i,j-1] * (-1/(h2*h2)) +
                                    matrix_apply_to[i,j+1] * (-1/(h2*h2)))
        end
    end
    @threads for i in 2:M-1
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

    @threads for j in 2:N-1
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

function get_B(M::Int, N::Int, h1::Float64, h2::Float64, A1::Float64, A2::Float64, B1::Float64, B2::Float64)::Matrix{Float64}
    result_matrix = zeros(Float64, M, N)
    @threads for i in 2:M-1
        @threads for j in 2:N-1
            result_matrix[i,j] = F(A1+ (i-1)*h1, B1+ (j-1)*h2)
        end
    end
    @threads for i in 2:M-1
        # it's (10) equations
        # i=1,M-1
        # top applying

        result_matrix[i, N] = psi(A1+ (i-1)*h1, B1+ (N-1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F(A1 + (i-1)*h1, B1 + (N-1)*h2)
        # bottom applying
        result_matrix[i, 1] = psi(A1+ (i-1)*h1, B1+ 0*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F(A1 + (i-1)*h1, B1 + 0*h2)
    end

    @threads for j in 2:N-1
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


function scalar_product(A::Matrix{Float64}, B::Matrix{Float64}, h1::Float64, h2::Float64)::Float64
    # sum = 0.0
    @assert size(A) == size(B)
    (M, N) = size(A, 1), size(A, 2)
    result_matrix = zeros(Float64, M, N)
    @threads for i in 1:M
        @threads for j in 1:N
            result_matrix[i,j] = h1*h2*ρ(i, M)*ρ(j, N)* A[i,j] * B[i,j]
        end
    end
    return sum(result_matrix)
end

function multiply_by_num(A::Matrix{Float64}, num::Float64)::Matrix{Float64}
    return A*num
end

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
    @threads for i in 1:M+1
        @threads for j in 1:N+1
            item = abs(omega_next[i,j]-u(h1*(i-1),h2*(j-1)))
            max_ = item > max_ ? item : max_
        end
    end
    @printf("diff:%.10f\n", max_)
end
end
