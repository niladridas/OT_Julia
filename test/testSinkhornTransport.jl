# Read a,b,K,U,lambda from dat file
a = readdlm("sinkhorn_nil/matlabreplica/a.dat")
b = readdlm("sinkhorn_nil/matlabreplica/b.dat")
K = readdlm("sinkhorn_nil/matlabreplica/K.dat")
U = readdlm("sinkhorn_nil/matlabreplica/U.dat")
Lam = readdlm("sinkhorn_nil/matlabreplica/lambda.dat")

using('sinkhorn_nil/matlabreplica/sinkhornTransport.jl')
D, u, v = sinkhornTransport(a,b,K,U,Lam,"marginalDifference",Inf,.5e-2,5000,0)
