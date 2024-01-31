# OperatorApproximation.jl
A framework for approximating functions, operators and solving operator equations

This is largely motivated by ApproxFun.jl and is a reimplemenation of many of those ideas for use by my research group.

An example:

```
using OperatorApproximation, Plots, SpecialFunctions

D = Derivative();
E = Evaluation();
X = Multiplication(x -> x);
B = BoundaryFunctional(canonicalBC(1,1)...)
L = E*D*D - X*E
R = 40
b = [airyai(-R);airyai(R)]
f = x -> 0
setbasis(Ultraspherical(0,ChebyshevMappedInterval(-R,R)));
setN("adaptive")
@time u = [B;L]\[b,f]
u.c

plot(u;dx=0.001)
```
