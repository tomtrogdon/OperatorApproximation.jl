# OperatorApproximation.jl
[![codecov](https://codecov.io/gh/tomtrogdon/OperatorApproximation.jl/graph/badge.svg?token=QXMLQ083L3)](https://codecov.io/gh/tomtrogdon/OperatorApproximation.jl) ![CI](https://github.com/github/docs/actions/workflows/CI.yml/badge.svg)

A framework for approximating functions, operators and solving operator equations

This is largely motivated by ApproxFun.jl and is a reimplemenation of many of those ideas for use by my research group.

An example:

```
using OperatorApproximation, Plots, SpecialFunctions

R = 30;
sp = Ultraspherical(0.0,UltraMappedInterval(-R,R,2.0)); 
sp2 = Ultraspherical(2.0,UltraMappedInterval(-R,R,2.0));
M = Multiplication(x -> x);
D = Derivative();
Op = D^2 - Conversion(sp2)*M
lbdry = FixedGridValues([-R],ChebyshevMappedInterval(-R,R)) |> Conversion;
rbdry = FixedGridValues([R],ChebyshevMappedInterval(-R,R)) |> Conversion;
setbasis(sp)
u = (lbdry ⊘ rbdry ⊘ Op)\[[airyai(-R)];[airyai(R)]; x->0]
u = u[1]

plot(u;dx=0.001)
```


The structure of the code is as follows:

* The core of the code is the Domain and Bases.   The Bases directory contains one file for each basis that is implemented.  This has the specific functions for performing all the core tasks for a given basis.  
* The implementation of operators occurs in two steps.  The first is the definition of an AbstractOperator that is then turned to a ConcreteOperator of an appropriate type when it acts on a Basis.
* The details of this conversion should be implemented for each AbstractOperator in a file in the subdirectory of the Operators directory corresponding to the basis under consideration.  When implementing a new operator, it should be categorized by basis on which it acts.
