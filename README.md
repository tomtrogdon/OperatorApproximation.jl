# OperatorApproximation.jl
A framework for approximating functions, operators and solving operator equations

This is largely motivated by ApproxFun.jl and is a reimplemenation of many of those ideas for use by my research group.

An example:

```
using OperatorApproximation, Plots, SpecialFunctions

R = 30;
sp = Ultraspherical(0.0,ChebyshevMappedInterval(-R,R));
gv = GridValues(ChebyshevMappedInterval(-R,R));
E = Conversion(gv);
M = Multiplication(x -> x);
D = Derivative();
Op = E*D^2 - M*E
Matrix(Op*sp,10,10)

lbdry = FixedGridValues([-R],ChebyshevMappedInterval(-R,R)) |> Conversion;
rbdry = FixedGridValues([R],ChebyshevMappedInterval(-R,R)) |> Conversion;

setbasis(sp)
setgrid(gv)

u = [lbdry;rbdry;Op]\[[airyai(-R)];[airyai(R)]; x->0]

plot(u;dx=0.001)
```


The structure of the code is as follows:

* The core of the code is the Domain and Bases.   The Bases directory contains one file for each basis that is implemented.  This has the specific functions for performing all the core tasks for a given basis.  
* The implementation of operators occurs in two steps.  The first is the definition of an AbstractOperator that is then turned to a ConcreteOperator of an appropriate type when it acts on a Basis.
* The details of this conversion should be implemented for each AbstractOperator in a file in the subdirectory of the Operators directory corresponding to the basis under consideration.  When implementing a new operator, it should be categorized by basis on which it acts.