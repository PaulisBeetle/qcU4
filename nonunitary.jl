using Yao
using Yao.ConstGate
using BitBasis

g = matblock(rand_unitary(1<<2))
reg = rand_state(2)

p0 = put(2, 1=>ConstGate.P0)
normalize_factor(reg) = probs(reg) |> sum

function run_()
    f = 1.0
    for i=1:10
        reg |> g |> p0
        f *= normalize_factor(reg)
        reg |> normalize!
    end
    return f
end
