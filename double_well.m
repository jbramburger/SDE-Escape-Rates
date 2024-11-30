% -------------------------------------------------------------------------
% Escape From a Double Well Potential 
%
% This script is used to replicate the results of Section 4.1 of Bounding 
% Escape Rates and Approximating Quasi-Stationary Distributions of 
% Brownian Dynamics by Jason J. Bramburger.
%
% This script bounds the principle eigenvalue and approximates the leading
% eigenfunction for a stochastic double well problem. The well is
% parametrized by a value alpha (denoted al below) and the domain is fixed
% to be -L <= x <= L. This script works with the rescaled domain -1 <= x'
% <= 1 by defining x = Lx'. 
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

%% System parameters

% Noise strength
sig = 1;

% Alpha parameter 
%  - values in manuscript range between 0 and 1
al = 1.0; 

% Length of space
% - value in manuscript is fixed to 3
L = 3; 

%% SOS Setup

% Degree
d = 10;

% Variables and objectives
sdpvar lam x u v % v = du/dx

% Witten Laplacian potential function
w = -0.5*( -al*(2*L*x - 1) + (2 - 3*L*x)*L*x ) + (1/(sig^2))*( (L*x)*(1 - L*x)*(L*x + al) )^2; 

% Variational formulation to be optimized
p = (sig^2)*v^2/L^2 + w*u^2 - lam*u^2;

% PDR function
[f,cf] = polynomial(x,d);
jacf = jacobian(f,x);
df = jacf*(1 - x^2) + 2*f*x; % derivative of f/(1-x^2) wrt x variable

% S-procedure
powers = monpowers(3,d+2);
powers = powers(find(powers(:,2)+powers(:,3)<=2),:); % max total degree in u,v is 2
% Build q1 polynomial
cs = sdpvar(size(powers,1),1); % s poylnomial coefficients
s = 0;
for i = 1:size(powers,1)
    s = s + cs(i)*prod([x u v].^powers(i,:));
end

%% Solve SOS problem

obj = [sos( p*(1-x^2)^2 + df*u^2 + 2*f*v*u*(1 - x^2)  - (1 - x^2)*s); sos(s)];
solvesos(obj,-lam,[],[cf; cs; lam]) 

%% Extract approximation of leading eigenfunction 

% Plotting variables
xend = 1 - 0.0001; % smaller than 1 to avoid singularity with f1, f2 functions
xplt = -xend:xend/1000:xend;
tstep = 1e-3;
t = 0:tstep:1;

% Recast sdpvar functions to evaluate meshgrid values
feval = sdisplay(replace(f,cf,value(cf)));
feval = feval{1};
feval = replace(feval,'^','.^');    % element-wise exponentiation
feval = replace(feval,'*','.*');    % element-wise multiplication
feval = str2func(['@(x)' feval]); % convert function string to function handle

loggradu = [];
logu = [];
for ind1 = 1:length(xplt)
   
    loggradu = xplt(ind1)*feval(xplt(ind1)*t)./( 1 - (xplt(ind1)*t).^2 );
    logu(ind1) = -trapz(t,loggradu); 
    
end

% Exponentiate to get the original u^*
ustar = exp(logu*L^2/sig^2);

% Undo symmetrization using weight function 
potential = (1/12)*(L*L*xplt.^2).*(-al*(6 - 4*L*xplt) + L*xplt.*(3*L*xplt - 4));

% Plot the potential landscape
figure(1)
plot(L*xplt,potential,'k','LineWidth',2)
xlabel('$x$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$V(x)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',16)

vstar = ustar.*exp(-0.5*potential/(sig^2/L^2));

% Plot the approximation of the leading eigenfunction
figure(2)
plot(L*xplt,vstar,'k','LineWidth',2)
xlabel('$x$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$u^*(x)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',16)

%% Use the Rayleigh Quotient to bound lambda^* from above

% L2 norm of u^* 
normustar = trapz(xplt,ustar.^2);

% L2 norm of du^*/dx
normvstar = (L^2/sig^2)*trapz(xplt, (ustar.*feval(xplt)./( 1 - xplt.^2 )).^2 );

% Recast sdpvar functions to evaluate meshgrid values
weval = sdisplay(w);
weval = weval{1};
weval = replace(weval,'^','.^');    % element-wise exponentiation
weval = replace(weval,'*','.*');    % element-wise multiplication
weval = str2func(['@(x)' weval]); % convert function string to function handle

% L2 norm of potential
normpotential = trapz(xplt, weval(xplt).*ustar.^2 );

% Approximate eigenvalue from above
upperBnd = (normvstar + normpotential)/normustar;

%% Print results
fprintf('Bounds using degree %i polynomials:\n',d)
fprintf('\t Upper bound = %f\n', upperBnd)
fprintf('\t Lower bound = %f\n', value(lam))

