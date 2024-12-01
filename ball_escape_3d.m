% -------------------------------------------------------------------------
% Escape From a Ball in 3D with Fully Absorbing Boundaries 
%
% This script is used to replicate the results of Section 4.2.1 of Bounding 
% Escape Rates and Approximating Quasi-Stationary Distributions of 
% Brownian Dynamics by Jason J. Bramburger.
%
% This script bounds the principle eigenvalue and approximates the leading
% eigenfunction for Brownian motion in the unit ball in 3D. The entire 
% boundary is taken to be absorbing and so the generator acts on functions
% with Dirichlet boundary conditions. The exact value of the leading
% eigenvalue is pi^2, the square of the first root of the 0th spherical
% Bessel function j_0(x) = sin(x)/x.
%
% Packages required: YALMIP and MOSEK
%
% Written by J. Bramburger.
%
% -------------------------------------------------------------------------

% Clean workspace
clear; close all; clc
yalmip clear
format long

%% Polynomial Optimization Setup

% Degree
d = 8;

% SDP variables
sdpvar x y z u ux uy uz lam % ux = du/dx, uy = du/dy, and  uz = du/dz 

% Create function f1(x,y) odd in x
fpowers1 = monpowers(3,d);
fpowers1 = fpowers1(find(mod(fpowers1(:,1),2)==1),:); % equivariant wrt x --> -x
fpowers1 = fpowers1(find(mod(fpowers1(:,2),2)==0),:); % invariant wrt y --> -y
fpowers1 = fpowers1(find(mod(fpowers1(:,3),2)==0),:); % invariant wrt z --> -z
% Build f1 polynomial
cf1 = sdpvar(size(fpowers1,1),1); % f1 poylnomial coefficients
f1 = 0;
for i = 1:size(fpowers1,1)
    f1 = f1 + cf1(i)*prod([x y z].^fpowers1(i,:));
end

% Create function f2(x,y) odd in y
fpowers2 = monpowers(3,d);
fpowers2 = fpowers2(find(mod(fpowers2(:,1),2)==0),:); % invariant wrt x --> -x
fpowers2 = fpowers2(find(mod(fpowers2(:,2),2)==1),:); % equivariant wrt y --> -y
fpowers2 = fpowers2(find(mod(fpowers2(:,3),2)==0),:); % invariant wrt z --> -z
% Build f2 polynomial
cf2 = sdpvar(size(fpowers2,1),1); % f2 poylnomial coefficients
f2 = 0;
for i = 1:size(fpowers2,1)
    f2 = f2 + cf2(i)*prod([x y z].^fpowers2(i,:));
end

% Create function f3(x,y) odd in z
fpowers3 = monpowers(3,d);
fpowers3 = fpowers3(find(mod(fpowers3(:,1),2)==0),:); % invariant wrt x --> -x
fpowers3 = fpowers3(find(mod(fpowers3(:,2),2)==0),:); % invariant wrt y --> -y
fpowers3 = fpowers3(find(mod(fpowers3(:,3),2)==1),:); % equivariant wrt z --> -z
% Build f3 polynomial
cf3 = sdpvar(size(fpowers3,1),1); % f3 poylnomial coefficients
f3 = 0;
for i = 1:size(fpowers2,1)
    f3 = f3 + cf3(i)*prod([x y z].^fpowers3(i,:));
end

% Derivative terms
dfx = jacobian(f1,x);
dfy = jacobian(f2,y);
dfz = jacobian(f3,z);
df = (1 - x^2 - y^2 - z^2)*dfx + 2*x*f1 + (1 - x^2 - y^2 - z^2)*dfy + 2*y*f2 + (1 - x^2 - y^2 - z^2)*dfz + 2*z*f3 ; % total derivative with rational ansatz 

% S-procedure
Spowers = monpowers(7,d+2);
Spowers = Spowers(find(mod(Spowers(:,1)+Spowers(:,5),2)==0),:); % invariant under (x,ux) --> -(x,ux)
Spowers = Spowers(find(mod(Spowers(:,2)+Spowers(:,6),2)==0),:); % invariant under (y,uy) --> -(y,uy)
Spowers = Spowers(find(mod(Spowers(:,3)+Spowers(:,7),2)==0),:); % invariant under (y,uz) --> -(y,uz)
Spowers = Spowers(find(Spowers(:,1)<=d+2),:); % max degree in x is d
Spowers = Spowers(find(Spowers(:,2)<=d+2),:); % max degree in y is d
Spowers = Spowers(find(Spowers(:,4)+Spowers(:,5)+Spowers(:,6)+Spowers(:,7)==2),:); % max total degree in u,ux,uy,uz is 2
% Build S polynomial
cs = sdpvar(size(Spowers,1),1); % s poylnomial coefficients
s = 0;
for i = 1:size(Spowers,1)
    s = s + cs(i)*prod([x y z u ux uy uz].^Spowers(i,:));
end

% Objective
p = ux^2 + uy^2 + uz^2 - lam*u^2;
obj = [sos( p*(1 - x^2 - y^2 - z^2)^2 + df*u^2 + 2*(ux*f1 + uy*f2 + uz*f3)*u*(1 - x^2 - y^2 - z^2) - (1 - x^2 - y^2 - z^2)*s ); sos(s)];
solvesos(obj,-lam,[],[cf1; cf2; cf3; cs; lam;]) 

%% Solve SOS problem

fprintf('Computed value of lambda: %f\n', value(lam))
fprintf('Compare with the exact value: %f\n', pi^2) % value should be 
 





