{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "fd63d8b7",
   "metadata": {},
   "outputs": [],
   "source": [
    "#This builds the derivative operator and applies it appropriately\n",
    "function *(D::Derivative,domain::Rational)\n",
    "    if D.order == 1\n",
    "        range = domain\n",
    "        dom = domain.GD.D\n",
    "        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? 2im*pi*j/(dom.b-dom.a) : 0im ))\n",
    "    else\n",
    "        Derivative(D.order-1)*(Derivative(1)*domain)\n",
    "    end\n",
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
