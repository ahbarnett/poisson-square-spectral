function [u info] = spectralfft2d(rhsfun,n)
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

% Barnett 3/6/23
if nargin==0, test_spectralfft2d; return; end
g = (1:n-1)/n;            % interior 1D grid (n-1 pts)
[xx yy] = ndgrid(g,g);
f = rhsfun(xx,yy);        % eval interior
fsym = zeros(2*n);
fsym(2:n,2:n) = f; fsym(n+2:end,2:n) = -flipud(f);  % unfold symm to L=2 box
fsym(2:n,n+2:end) = -fliplr(f); fsym(n+2:end,n+2:end) = flipud(fliplr(f));
kg = (-n:n-1)*pi;         % full k-grid assoc to FFT, spacing 2pi/L = pi
[kx ky] = ndgrid(kg,kg);
ik2 = 1./(kx.^2+ky.^2);   % spectral filter ||k||^{-2} on full k-grid
ik2(n+1,n+1) = 1.0;       % innocuous value for k = 0-vector
info.fhat = ifft2(fsym);  % Euler-Fourier formula quad w/ h=1/n per dim
usym = fft2(info.fhat .* fftshift(ik2));       % Poisson solve
usym = real(usym);
u = usym(1:n+1,1:n+1);    % extract 1st quarter, including BC vals
info.f = fsym; info.u = usym; info.kg = kg;   % dump other info

%%%%%%%
function test_spectralfft2d

disp('fig1: single-freq analytic test...')
n = 17;    % small and arb
rhsfun = @(x,y) sin(pi*x).*sin(3*pi*y);     % odd-reflection symm, single freq
[u info] = spectralfft2d(rhsfun,n);
g = (0:n)/n; [xx yy] = ndgrid(g,g);
uex = rhsfun(xx,yy) / (pi^2 + (3*pi)^2);    % exact soln for single freq
% u(2,2) / uex(2,2)  % debug
fprintf("u rel err vs analytic = %.3g\n", norm(u(:)-uex(:))/norm(u(:)))
gsym = (0:2*n-1)/n;            % 1d grid for unfolded driving func fsym
figure(1); clf; set(gcf,'position',[200 900 1400 300]);
subplot(1,4,1);
imagesc(gsym,gsym,info.f'); xlabel('x'); ylabel('y');
hold on; plot([0 1 1 0 0], [0 0 1 1 0], 'k--'); title('fsym and domain');
axis tight equal xy; colorbar
subplot(1,4,2);
kg = info.kg;
imagesc(kg,kg,abs(fftshift(info.fhat.'))); xlabel('k_x'); ylabel('k_y');
axis tight equal xy; colorbar; title('f Fourier coeffs')
subplot(1,4,3);
imagesc(gsym,gsym,info.u'); xlabel('x'); ylabel('y');
hold on; plot([0 1 1 0 0], [0 0 1 1 0], 'k--'); title('usym and domain');
axis tight equal xy; colorbar
subplot(1,4,4);
imagesc(g,g,u'); xlabel('x'); ylabel('y'); title('u as returned');
axis tight equal xy; colorbar

disp('fig2: conv rate on C-infty compact-eps-supp f...')
rhsfun = @(x,y) exp(-200*((x-0.6).^2 + (y-0.55).^2));  % 1e-14 symm trunc err
ns = 10:10:100;        % conv test, all n must be even
us = nan*ns;
for i=1:numel(ns), n=ns(i);
  u = spectralfft2d(rhsfun,n);
  us(i) = u(n/2+1,n/2+1);       % extract center val; always at same (0.5,0.5)
end
figure(2); clf; set(gcf,'position',[200 550 1000 300]);
subplot(1,3,1); semilogy(ns(1:end-1),abs(us(1:end-1)-us(end))/us(end),'+-');
xlabel('n'); ylabel('est u err'); title('u ptwise conv');
subplot(1,3,2); g = (0:n)/n; [xx yy] = ndgrid(g,g);
imagesc(g,g,rhsfun(xx,yy)'); xlabel('x'); ylabel('y'); title('f')
axis tight equal xy; colorbar
subplot(1,3,3); g = (0:n)/n; [xx yy] = ndgrid(g,g);
imagesc(g,g,u'); xlabel('x'); ylabel('y'); title('u');
axis tight equal xy; colorbar

disp('fig3: conv rate on C-infty f whose reflection to [0,2]^2 non-smooth...')
rhsfun = @(x,y) 1 + 0*x;     % identically 1; reflects to a "1,-1,1,-1 4-junc"
ns = 2.^(4:10);        % conv test, all n must be even, now geom spaced
us = nan*ns;
for i=1:numel(ns), n=ns(i);
  [u info] = spectralfft2d(rhsfun,n);
  us(i) = u(3*n/4+1,5*n/8+1);   % extract val at same (0.75,0.625), offcenter
end
% note we could also use analytic torsion soln (sum of sinh*cosh) here...
figure(3); clf; set(gcf,'position',[200 200 1400 300]); ns = ns(1:end-1);
subplot(1,4,1); loglog(ns,abs(us(1:end-1)-us(end))/us(end),'+-');
xlabel('n'); ylabel('est u err'); hold on; plot(ns,ns.^-2,'r--');
legend('u err', '1/n^2'); title('u ptwise conv');
subplot(1,4,2); g = (0:n)/n; [xx yy] = ndgrid(g,g);
imagesc(g,g,rhsfun(xx,yy)'); xlabel('x'); ylabel('y'); title('f')
axis tight equal xy; colorbar
subplot(1,4,3); kg = fftshift(info.kg); j = 2:2:n; % k-indices to slice in 1D
loglog(kg(j), abs(info.fhat(2,j)), '.-'); xlabel('k_y'); ylabel('|fhat|');
hold on; plot(kg(j),kg(j).^(-1),'r--'); legend('abs fhat','1/k');
title('f Fourier coeff decay')
subplot(1,4,4); g = (0:n)/n; [xx yy] = ndgrid(g,g);
imagesc(g,g,u'); xlabel('x'); ylabel('y'); title('u');
axis tight equal xy; colorbar
%
% We see fhat ~ 1/k, so that uhat ~ 1/k^3, which matches Boussinesq's
% analytic series for Poiseuille flow in rectangular pipe.
% the tail sum of such appears to give u error of 1/kmax^2, ie 1/n^2.
%
% Question is: what is case for 2D solver based on product Cheby nodes?
% (Fortunato-Townsend 2020 IMAJNA).
