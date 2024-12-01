% -------------------------------------------------------------------------
% Escape From a Ball in 2D with Fully Absorbing Boundaries 
%
% This script is used to replicate the results of Section 4.2.1 of Bounding 
% Escape Rates and Approximating Quasi-Stationary Distributions of 
% Brownian Dynamics by Jason J. Bramburger.
%
% This script bounds the principle eigenvalue and approximates the leading
% eigenfunction for Brownian motion in the unit ball in 2D. The entire 
% boundary is taken to be absorbing and so the generator acts on functions
% with Dirichlet boundary conditions. The exact value of the leading
% eigenvalue is 5.7832, the square of the first root of the 0th Bessel
% function J_0.
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
d = 6;

% SDP variables
sdpvar x y u v w lam % v = du/dx and w = du/dy

% Create function f1(x,y) odd in x
fpowers1 = monpowers(2,d);
fpowers1 = fpowers1(find(mod(fpowers1(:,1),2)==1),:); % equivariant wrt x --> -x
fpowers1 = fpowers1(find(mod(fpowers1(:,2),2)==0),:); % invariant wrt y --> -y
% Build f1 polynomial
cf1 = sdpvar(size(fpowers1,1),1); % f1 poylnomial coefficients
f1 = 0;
for i = 1:size(fpowers1,1)
    f1 = f1 + cf1(i)*prod([x y].^fpowers1(i,:));
end

% Create function f2(x,y) odd in y
fpowers2 = monpowers(2,d);
fpowers2 = fpowers2(find(mod(fpowers2(:,2),2)==1),:); % equivariant wrt y --> -y
fpowers2 = fpowers2(find(mod(fpowers2(:,1),2)==0),:); % invariant wrt x --> -x
% Build f2 polynomial
cf2 = sdpvar(size(fpowers2,1),1); % f2 poylnomial coefficients
f2 = 0;
for i = 1:size(fpowers2,1)
    f2 = f2 + cf2(i)*prod([x y].^fpowers2(i,:));
end

% Derivative terms
dfx = jacobian(f1,x);
dfy = jacobian(f2,y);
df = (1 - x^2 - y^2)*dfx + 2*x*f1 + (1 - x^2 - y^2)*dfy + 2*y*f2; % total derivative with rational ansatz 

% S-procedure
Spowers = monpowers(5,d+2);
Spowers = Spowers(find(mod(Spowers(:,1)+Spowers(:,4),2)==0),:); % invariant under (x,v) --> -(x,v)
Spowers = Spowers(find(mod(Spowers(:,2)+Spowers(:,5),2)==0),:); % invariant under (y,w) --> -(y,w)
Spowers = Spowers(find(Spowers(:,1)<=d),:); % max degree in x is d
Spowers = Spowers(find(Spowers(:,2)<=d),:); % max degree in y is d
Spowers = Spowers(find(Spowers(:,3)+Spowers(:,4)+Spowers(:,5)<=2),:); % max total degree in u,v, w is 2
% Build S polynomial
cs = sdpvar(size(Spowers,1),1); % s poylnomial coefficients
s = 0;
for i = 1:size(Spowers,1)
    s = s + cs(i)*prod([x y u v w].^Spowers(i,:));
end

%% Solve SOS problem
p = v^2 + w^2 - lam*u^2;
obj = [sos( p*(1 - x^2 - y^2)^2 + df*u^2 + 2*(v*f1 + w*f2)*u*(1 - x^2 - y^2) - (1 - x^2 - y^2)*s ); sos(s)];
solvesos(obj,-lam,[],[cf1; cf2; cs; lam;]) 

%% Extract and plot approximate eigenfunction 

% Radial and phase variables
rend = 0.999999;
r = 0:rend/1000:rend;
th = 0:pi/500:2*pi;
xplt = r'*cos(th);
yplt = r'*sin(th);
tstep = 1e-3;
t = 0:tstep:1;

% Recast sdpvar functions to evaluate meshgrid values
f1eval = sdisplay(replace(f1,cf1,value(cf1)));
f1eval = f1eval{1};
f1eval = replace(f1eval,'^','.^');    % element-wise exponentiation
f1eval = replace(f1eval,'*','.*');    % element-wise multiplication
f1eval = str2func(['@(x,y)' f1eval]); % convert function string to function handle
% same for f2
f2eval = sdisplay(replace(f2,cf2,value(cf2)));
f2eval = f2eval{1};
f2eval = replace(f2eval,'^','.^');    % element-wise exponentiation
f2eval = replace(f2eval,'*','.*');    % element-wise multiplication
f2eval = str2func(['@(x,y)' f2eval]); % convert function string to function handle

loggradu = [];
logu = [];
for ind1 = 1:length(r)
    for ind2 = 1:length(th)
        loggradu = xplt(ind1,ind2)*f1eval(xplt(ind1,ind2)*t, yplt(ind1,ind2)*t)./(1 - (xplt(ind1,ind2)*t).^2 - (yplt(ind1,ind2)*t).^2) + yplt(ind1,ind2)*f2eval(xplt(ind1,ind2)*t, yplt(ind1,ind2)*t)./(1 - (xplt(ind1,ind2)*t).^2 - (yplt(ind1,ind2)*t).^2);
        logu(ind1,ind2) = -trapz(t,loggradu);
    end
end

% Exponentiate to get the original u^*
ustar = exp(logu);

% Plot the result
figure(1)
surf(xplt,yplt,ustar)
view(0,90)
shading interp
axis off
colorbar
set(gca,'FontSize',16)

%% Use Rayleigh quotient to bound eigenvalue from above

% Compute L2 norm of u^* using polar coordinates
temp = [];
for ind1 = 1:length(r)
    temp(ind1) = trapz(th,r(ind1)*ustar(ind1,:).^2);
end
normustar = trapz(r,temp);

% Compute L2 norm of du^*/dx using polar coordinates
temp = [];
for ind1 = 1:length(r)
    vint = ustar(ind1,:).*f1eval(xplt(ind1,:), yplt(ind1,:))./(1 - (xplt(ind1,:)).^2 - (yplt(ind1,:)).^2);
    temp(ind1) = trapz(th,r(ind1)*vint.^2);
end
normvstar = trapz(r,temp);

% Compute L2 norm of du^*/dy using polar coordinates
temp = [];
for ind1 = 1:length(r)
    wint = ustar(ind1,:).*f2eval(xplt(ind1,:), yplt(ind1,:))./(1 - (xplt(ind1,:)).^2 - (yplt(ind1,:)).^2);
    temp(ind1) = trapz(th,r(ind1)*wint.^2);
end
normwstar = trapz(r,temp);

upperBnd = (normvstar + normwstar)/normustar;

%% Print results

fprintf('Bounds using degree %i polynomials:\n',d)
fprintf('\t Computed lower bound: %f\n', value(lam))
fprintf('\t Computed upper bound: %f\n', upperBnd)
fprintf('Compare with the exact value: %f\n', 2.4048255577^2)% value should be 2.4048^2 = 5.783063039999999







