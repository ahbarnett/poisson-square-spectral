function fft_vs_cheb()
% script to compare spectralfft2d vs chebfun2.poisson of Fortunato-Townsend.
%
% This loops over four RHS functions and records both L^infty and L^2 errors.
% It is packaged as function.
% Author: Dan Fortunato, March 2023.

rng(0)
ns = 2.^(2:6);  % convergence study

clf
tiledlayout(1, 4);
set(gcf, 'position', [500 700 1600 400]);

f = chebfun2(@(x,y) 1+0*x, [0 1 0 1]);
makeplot(f, ns, 'One');

f = randnfun2(1, [0 1 0 1]);     % a smooth random chebfun2 function
makeplot(f, ns, 'Random');

f = chebfun2(@(x,y) exp(-200*((x-0.6).^2 + (y-0.55).^2)), [0 1 0 1]);
makeplot(f, ns, 'Compactly supported');  % same RHS as fig2 of spectralfft2d

x = chebfun2(@(x,y) x, [0 1 0 1]);
y = chebfun2(@(x,y) y, [0 1 0 1]);
u = x.*(1-x).*y.*(1-y).*randnfun2(1, [0 1 0 1]);
f = -lap(u);
makeplot(f, ns, 'Smooth solution');

nexttile(1)
lw = 1;
hold on
loglog(ns, 0.5*ns.^-2, 'r--', 'linewidth', lw)
loglog(ns, ns.^-4,     'b--', 'linewidth', lw)
loglog(ns, 20*ns.^-6,  'k--', 'linewidth', lw)
hold off

leg = legend('FFT, 2-norm', 'FFT, \infty-norm', 'Chebyshev, 2-norm', 'Chebyshev, \infty-norm', '1/n^2', '1/n^4', '1/n^6');
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';

shg  % brings current window forward

end

function makeplot(f, ns, ttl)
% Make convergence plot for RHS f, in a sequence of tiled plot.
% Inputs:
%   f - chebfun2 object giving the RHS on [0,1]^2.
%   ns - list of numbers of nodes (per dimension) for convergence study
%   ttl - title string for the tiled plot

    % Adaptively compute a "true" solution:
    utrue = chebfun2.poisson(-f);

    nns = length(ns);
    % FFT soln (called u)
    u = cell(length(ns), 1);
    for k = 1:nns
      t=tic;
      [~, info] = spectralfft2d(@(x,y) f(x,y), ns(k));
      tu = toc(t);
      usym = chebfun2(info.u.', [0 2 0 2], 'trig');
      u{k} = restrict(usym, [0 1 0 1]);
    end

    % Chebyshev soln (called v)
    v = cell(length(ns), 1);
    for k = 1:length(ns)
      t=tic;
      v{k} = chebfun2.poisson(-f, 0, ns(k));     % use ns(k) nodes
      tv = toc(t);
    end

    fprintf("timing of last run: FFT(u) %.3g s, vs Cheb(v) %.3g s (FFT is %.3g times faster)\n", tu, tv, tv/tu)
    
    % Compute error
    uerr_2 = zeros(length(ns), 1); uerr_inf = zeros(length(ns), 1);
    verr_2 = zeros(length(ns), 1); verr_inf = zeros(length(ns), 1);
    normU_2 = norm(utrue, 2);
    normU_inf = norm(utrue, inf);
    for k = 1:length(ns)
        uerr_2(k) = norm(u{k} - utrue, 2) / normU_2;
        verr_2(k) = norm(v{k} - utrue, 2) / normU_2;
        uerr_inf(k) = norm(u{k} - utrue, inf) / normU_inf;
        verr_inf(k) = norm(v{k} - utrue, inf) / normU_inf;
    end

    nexttile
    uerr_2(uerr_2 == 0) = 1e-100;
    verr_2(verr_2 == 0) = 1e-100;
    uerr_inf(uerr_inf == 0) = 1e-100;
    verr_inf(verr_inf == 0) = 1e-100;
    lw = 1;
    loglog(ns, abs(uerr_2),   'ro-', 'markerfacecolor', 'r', 'linewidth', lw), hold on
    loglog(ns, abs(uerr_inf), 'r^-', 'markerfacecolor', 'r', 'linewidth', lw)
    loglog(ns, abs(verr_2),   'bo-', 'markerfacecolor', 'b', 'linewidth', lw)
    loglog(ns, abs(verr_inf), 'b^-', 'markerfacecolor', 'b', 'linewidth', lw)
    ylim([1e-14 1])
    hold off
    xlabel('N')
    ylabel('Relative error')
    title(ttl)
    axis square
    set(gca, 'fontsize', 14)
    set(gca, 'linewidth', lw)

end
