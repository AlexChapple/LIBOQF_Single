function modulo_func(z)

    return sqrt.(real.(z).^2 + imag.(z).^2)


end 
# print(sum((modulo_func(g_1).^2) .+ (modulo_func(e_1).^2))) 

function main()

    g_1 = [1 2 3]
    e_1 = [1 2 3]

    total :: Float64 = 0 
    total += sum((modulo_func(g_1).^2) .+ (modulo_func(e_1).^2)) 

    print(total)

end 

function main2()

    a = spzeros(ComplexF64, 10)

    a[1] += (1 + 2im)
    a[2] += (2 + 3im)

    print(a)

end
using SparseArrays

function main3()

    a = rand(Float64, 10)
    print(a)

end 

main3()