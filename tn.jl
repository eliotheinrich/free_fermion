using ITensors

function create_1d_tensor_network(L::Int, d::Int, D::Int)
    # Create arrays to hold indices and tensors
    s = [Index(d, "s$i") for i in 1:L]          # External (site) indices
    l = [Index(D, "l$i") for i in 1:(L-1)]      # Internal (link) indices
    A = Vector{ITensor}(undef, L)               # The tensor network

    for i in 1:L
        if i == 1
            # First tensor has one site and one internal (right) index
            A[i] = randomITensor(s[i], l[i])
        elseif i == L
            # Last tensor has one site and one internal (left) index
            A[i] = randomITensor(s[i], l[i-1])
        else
            # Middle tensors have one site and two internal (left and right) indices
            A[i] = randomITensor(s[i], l[i-1], l[i])
        end
    end

    return A
end

using LinearAlgebra

function fermionic_creation_op(i::Int, L::Int)
    @assert 1 ≤ i ≤ L "Site index i must be in 1..L"

    Z  = [ 1.0  0.0; 0.0 -1.0]
    I  = [ 1.0  0.0; 0.0  1.0]
    sp = [ 0.0  1.0; 0.0  0.0]

    # Build the operator list
    ops = Vector{Matrix{Float64}}(undef, L)
    for j in 1:L
        if j < i
            ops[j] = Z
        elseif j == i
            ops[j] = sp
        else
            ops[j] = I
        end
    end

    op = ops[1]
    for j in 2:L
        op = kron(op, ops[j])
    end

    return op
end

L = 2
i = 2
c_dag = fermionic_creation_op(i, L)

println("Fermionic creation operator c†_$i for L = $L qubits:\n")
for row in 1:size(c_dag, 1)
    for col in 1:size(c_dag, 2)
        print(rpad(c_dag[row, col], 8), " ")
    end
    println()
end

# Example usage
L = 5      # Length of the tensor network
d = 2      # Physical index dimension
D = 3      # Internal bond dimension

A = create_1d_tensor_network(L, d, D)

function contract_tensor_network(A::Vector{ITensor})
    T = A[1]
    for i in 2:length(A)
        T *= A[i]
    end
    return T
end

T = contract_tensor_network(A)

print(T)


# Optionally, print the tensors and their indices
#for (i, tensor) in enumerate(A)
#    println("A[$i] has indices: ", inds(tensor))
#end

