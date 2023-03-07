# poisson-square-spectral

MATLAB demo of FFT fast solver for Poisson's equation with homogeneous Dirichlet BCs in the square.

Run `spectralfft2d` without arguments, which performs a set of demo tests
and produces the following figure outputs. The function is documented below.


### Demo outputs

1.

![fig 1: analytic solution at single frequency](figs/fig1.png)

2.

![fig 2: spectral convergence for numerically localized Gaussian forcing](figs/fig2.png)

3.

![fig 3: algebraic convergence for constant forcing, due to weak corner singularities](figs/fig3.png)





### Documentation

```
% SPECTRALFFT2D  Poisson solver for Dirichlet on square via spectral FFT
%
% u = spectralfft2d(rhsfun,n) returns u values on the (n+1)^2 grid defined
%  below, given handle to rhsfun for the right-hand side of form
%  f = rhsfun(x,y), on unit square domain D = [0,1]^2, with homogeneous
%  Dirichlet BCs. u solves the Poisson equation -Delta u = f with these BCs.
%
% The returned 2D grid definition is (i/n,j/n) for i,j=0,..,n. The solution
%  is zero at 1st or last point in either dimension, as per BCs.
%
% [u info] = spectralfft2d(f,n) also returns a debug struct info giving:
%  info.f    - rhsfun on 2n*2n expanded grid. Note: x fast, y slow.
%  info.u    - u on 2n*2n expanded grid.
%  info.fhat - 2D DFT of f
%  info.kg   - 1D k-grid
%
% Notes:
% 1) this answers the question in Fortunato-Townsend 2020 IMAJNA about
%    existence of a fast spectrally-accurate
%    solver for the Dirichlet square, in a simpler way (see 1st, 2nd tests).
% 2) Convergence should deteriorate to merely algebraic
%    when f does not reflect to [0,2]^2 as a smooth function (see 3rd test).
```

