{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "fd63d8b7",
   "metadata": {},
   "outputs": [],
   "source": [
    "#I think this just says we can go from a Rational to a DiscreteBasis\n",
    "function isconvertible(b1::Rational,b2::DiscreteBasis)  # probably can be simplified\n",
    "    iscompatible(b1.GD,b2.GD)\n",
    "end\n",
    "\n",
    "#I think this just says we can go from a Rational basis to another Rational basis\n",
    "function isconvertible(b1::Rational,b2::Rational)\n",
    "    iscompatible(b1.GD,b2.GD)\n",
    "end\n",
    "\n",
    "#This converts the Rational basis to values on the grid; note that GridValues is a DiscreteBasis\n",
    "function conversion(b1::Rational,b2::GridValues)\n",
    "    basegrid =  n -> b2.GD.grid(n)\n",
    "    # In principle, we need to do this:\n",
    "    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))\n",
    "    # but we are checking that the two grid domains are compatible\n",
    "    # and currently this forces the composition of the maps to\n",
    "    # be the identity\n",
    "    \n",
    "    #!!!!!!!!!!!!!!Need to change to some kind of RationalEvaluationOperator???\n",
    "    Op = FourierEvaluationOperator(basegrid) #QUESTION: What does \"FourierEvalutationOperator\" do and where defined?\n",
    "    ConcreteOperator(b1,b2,Op)\n",
    "end\n",
    "\n",
    "#This does the same as above but for fixed values on the grid\n",
    "function conversion(b1::Rational,b2::FixedGridValues)\n",
    "    # See conversion remark above.\n",
    "    Op = FixedGridFourierEvaluationOperator(b2.pts) #QUESTION: Where does \"FixedGridFourierEvaluationOperator\" live?\n",
    "    ConcreteOperator(b1,b2,Op)\n",
    "end\n",
    "\n",
    "#This converts the Fourier basis to itself\n",
    "function conversion(b1::Rational,b2::Rational)\n",
    "    # TODO:  identity operator\n",
    "    ConcreteOperator(b1,b2,IdentityOperator())\n",
    "end"
   ]
  }
 ],
 "metadata": {
  "@webio": {
   "lastCommId": null,
   "lastKernelId": null
  },
  "kernelspec": {
   "display_name": "Julia 1.8.5",
   "language": "julia",
   "name": "julia-1.8"
  },
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.8.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
