% -------------------------------------------------------------------------
% Escape From a Ball in 2D with 2 Holes at (+/-1,0) 
%
% This script is used to replicate the results of Section 4.2.2 of Bounding 
% Escape Rates and Approximating Quasi-Stationary Distributions of 
% Brownian Dynamics by Jason J. Bramburger.
%
% We bound the principle eigenvalue and approximates the leading
% eigenfunction for Brownian motion in the unit ball in 2D that can only
% escape through two small balls of radius r centred at (1,0) and (0,1) on
% the boundary.
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

%% Variables and global polynomial optimization info

% Degree
d = 20;

% SDP variables
sdpvar x y u v w lam % v = du/dx and w = du/dy

% Poincare inequality portion
p = v^2 + w^2 - lam*u^2;

% Need to break up into two parts: -b <= y <= b inside circle and complement
% Compute b based on r
r = 0.5;
xb = (2 - r^2)/2;
b = sqrt(1 - xb^2);

%% Inside -b <= y <= b (i.e. b^2 - y^2 >= 0)

% Create function f1(x,y) odd in x and even in y
fpowers1 = monpowers(2,d-2);
fpowers1 = fpowers1(find(mod(fpowers1(:,1),2)==1),:); % equivariant wrt x --> -x
fpowers1 = fpowers1(find(mod(fpowers1(:,2),2)==0),:); % invariant wrt y --> -y
% Build f1 polynomial
cf1 = sdpvar(size(fpowers1,1),1); % f1 poylnomial coefficients
f1 = 0;
for i = 1:size(fpowers1,1)
    f1 = f1 + cf1(i)*prod([x y].^fpowers1(i,:));
end

% Create function f2(x,y) odd in y and even in x
fpowers2 = monpowers(2,d-2);
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
Spowers = monpowers(5,d);
Spowers = Spowers(find(mod(Spowers(:,1)+Spowers(:,4),2)==0),:); % invariant under (x,v) --> -(x,v)
Spowers = Spowers(find(mod(Spowers(:,2)+Spowers(:,5),2)==0),:); % invariant under (y,w) --> -(y,w)
Spowers = Spowers(find(Spowers(:,1)<=d),:); % max degree in x is d
Spowers = Spowers(find(Spowers(:,2)<=d),:); % max degree in y is d
Spowers = Spowers(find(Spowers(:,3)+Spowers(:,4)+Spowers(:,5)==2),:); % max total degree in u,v, w is 2

% Build S1 polynomial
cs1 = sdpvar(size(Spowers,1),1); % s poylnomial coefficients
s1 = 0;
for i = 1:size(Spowers,1)
    s1 = s1 + cs1(i)*prod([x y u v w].^Spowers(i,:));
end

% Build S2 polynomial
cs2 = sdpvar(size(Spowers,1),1); % s poylnomial coefficients
s2 = 0;
for i = 1:size(Spowers,1)
    s2 = s2 + cs2(i)*prod([x y u v w].^Spowers(i,:));
end

% Objective
obj = [sos( p*(1 - x^2 - y^2)^2 + df*u^2 + 2*(v*f1 + w*f2)*u*(1 - x^2 - y^2) - (1 - x^2 - y^2)*s1 - (b^2 - y^2)*s2 ); sos(s1); sos(s2)];

%% Outside -b <= y <= b (i.e. y^2 - b^2 >= 0)

% Create function g1(x,y) odd in x and even in y
gpowers1 = monpowers(2,d);
gpowers1 = gpowers1(find(mod(gpowers1(:,1),2)==1),:); % equivariant wrt x --> -x
gpowers1 = gpowers1(find(mod(gpowers1(:,2),2)==0),:); % invariant wrt y --> -y
% Build g1 polynomial
cg1 = sdpvar(size(gpowers1,1),1); % g1 poylnomial coefficients
g1 = 0;
for i = 1:size(gpowers1,1)
    g1 = g1 + cg1(i)*prod([x y].^gpowers1(i,:));
end

% Create function g2(x,y) odd in y and even in x
gpowers2 = monpowers(2,d);
gpowers2 = gpowers2(find(mod(gpowers2(:,2),2)==1),:); % equivariant wrt y --> -y
gpowers2 = gpowers2(find(mod(gpowers2(:,1),2)==0),:); % invariant wrt x --> -x
% Build g2 polynomial
cg2 = sdpvar(size(gpowers2,1),1); % g2 poylnomial coefficients
g2 = 0;
for i = 1:size(gpowers2,1)
    g2 = g2 + cg2(i)*prod([x y].^gpowers2(i,:));
end

% Derivative terms
dgx = jacobian(g1,x);
dgy = jacobian(g2,y);
dg = dgx + dgy; % total derivative

% Build S3 polynomial
cs3 = sdpvar(size(Spowers,1),1); % s poylnomial coefficients
s3 = 0;
for i = 1:size(Spowers,1)
    s3 = s3 + cs3(i)*prod([x y u v w].^Spowers(i,:));
end

% Build S4 polynomial
cs4 = sdpvar(size(Spowers,1),1); % s poylnomial coefficients
s4 = 0;
for i = 1:size(Spowers,1)
    s4 = s4 + cs4(i)*prod([x y u v w].^Spowers(i,:));
end

% Objective part 2
obj = [obj; sos( p + dg*u^2 + 2*(v*g1 + w*g2)*u - (1 - x^2 - y^2)*s3 - (y^2 - b^2)*s4 ); sos(s3); sos(s4)];

%% Additional condition on g1 and g2 to make them satisfy Neumann conditions

% Neuman boundary condition
z = 2*x*g1 + 2*y*g2;

% S-procedure 
Spowers2 = monpowers(2,d);
Spowers2 = Spowers2(find(mod(Spowers2(:,1),2)==0),:); % invariant under x --> -x
Spowers2 = Spowers2(find(mod(Spowers2(:,2),2)==0),:); % invariant under y --> -y

% Build S5 polynomial
cs5 = sdpvar(size(Spowers2,1),1); % s poylnomial coefficients
s5 = 0;
for i = 1:size(Spowers2,1)
    s5 = s5 + cs5(i)*prod([x y].^Spowers2(i,:));
end

% Build S6 polynomial
cs6 = sdpvar(size(Spowers2,1),1); % s poylnomial coefficients
s6 = 0;
for i = 1:size(Spowers2,1)
    s6 = s6 + cs6(i)*prod([x y].^Spowers2(i,:));
end
% Build S7 polynomial
cs7 = sdpvar(size(Spowers2,1),1); % s poylnomial coefficients
s7 = 0;
for i = 1:size(Spowers2,1)
    s7 = s7 + cs7(i)*prod([x y].^Spowers2(i,:));
end

% Build S8 polynomial
cs8 = sdpvar(size(Spowers2,1),1); % s poylnomial coefficients
s8 = 0;
for i = 1:size(Spowers2,1)
    s8 = s8 + cs8(i)*prod([x y].^Spowers2(i,:));
end

obj = [obj; sos( z - (1 - x^2 - y^2)*s5 - (y^2 - b^2)*s6 ); sos( -z - (1 - x^2 - y^2)*s7 - (y^2 - b^2)*s8 ); sos(s6); sos(s8)];

%% Continuity of f_i and g_i at y = b

f1b = replace(f1,y,b);
f2b = replace(f2,y,b);
g1b = replace(g1,y,b);
g2b = replace(g2,y,b);

[s9, cs9] = polynomial(x,d,0,1);
s9 = x*s9;
[s10, cs10] = polynomial(x,d,0,1);
s10 = x*s10;
[s11, cs11] = polynomial(x,d,0,1);
s11 = x*s11;
[s12, cs12] = polynomial(x,d,0,1);
s12 = x*s12;

obj = [obj; sos( f1b - g1b*(1 - x^2 - b^2) - (xb^2 - x^2)*s9 ); sos( -f1b + g1b*(1 - x^2 - b^2) - (xb^2 - x^2)*s10 ); sos( f2b - g2b*(1 - x^2 - b^2) - (xb^2 - x^2)*s11 ); sos( -f2b + g2b*(1 - x^2 - b^2) - (xb^2 - x^2)*s12 ); sos(s9); sos(s10); sos(s11); sos(s12)];

%% Solve the problem

sol = solvesos(obj,-lam,[],[cf1; cf2; cg1; cg2; cs1; cs2; cs3; cs4; cs5; cs6; cs7; cs8; cs9; cs10; cs11; cs12; lam;]) 

%% Extract approximate eigenfunction 

% Radial and phase variables
rend = 0.999999; % smaller than 1 to avoid singularity with f1, f2 functions
r = 0:rend/1200:rend;
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
% same for g1
g1eval = sdisplay(replace(g1,cg1,value(cg1)));
g1eval = g1eval{1};
g1eval = replace(g1eval,'^','.^');    % element-wise exponentiation
g1eval = replace(g1eval,'*','.*');    % element-wise multiplication
g1eval = str2func(['@(x,y)' g1eval]); % convert function string to function handle
% same for g2
g2eval = sdisplay(replace(g2,cg2,value(cg2)));
g2eval = g2eval{1};
g2eval = replace(g2eval,'^','.^');    % element-wise exponentiation
g2eval = replace(g2eval,'*','.*');    % element-wise multiplication
g2eval = str2func(['@(x,y)' g2eval]); % convert function string to function handle

% Continuity condition
contint = b*f2eval(0*t, b*t)./(1 - (0*t).^2 - (b*t).^2);
contcond = -trapz(t,contint);

loggradu = [];
logu = [];
for ind1 = 1:length(r)
    for ind2 = 1:length(th)
        if yplt(ind1,ind2) <= b && yplt(ind1,ind2) >= -b
            
            loggradu = xplt(ind1,ind2)*f1eval(xplt(ind1,ind2)*t,yplt(ind1,ind2)*t)./(1 - (xplt(ind1,ind2)*t).^2 - (yplt(ind1,ind2)*t).^2) + yplt(ind1,ind2)*f2eval(xplt(ind1,ind2)*t, yplt(ind1,ind2)*t)./(1 - (xplt(ind1,ind2)*t).^2 - (yplt(ind1,ind2)*t).^2);
            logu(ind1,ind2) = -trapz(t,loggradu);
            
        elseif yplt(ind1,ind2) > b 
        
            loggradu = xplt(ind1,ind2)*g1eval(xplt(ind1,ind2)*t,b + (yplt(ind1,ind2) - b)*t) + (yplt(ind1,ind2) - b)*g2eval(xplt(ind1,ind2)*t, b + (yplt(ind1,ind2) - b)*t);
            logu(ind1,ind2) = contcond - trapz(t,loggradu);
            
        elseif yplt(ind1,ind2) < -b 
        
            loggradu = xplt(ind1,ind2)*g1eval(xplt(ind1,ind2)*t,-b + (yplt(ind1,ind2) + b)*t) + (yplt(ind1,ind2) + b)*g2eval(xplt(ind1,ind2)*t, -b + (yplt(ind1,ind2) + b)*t);
            logu(ind1,ind2) = contcond - trapz(t,loggradu);
            
        end
        
    end
end

% Exponentiate to get the original u^*
ustar = exp(logu);

% Surface plot of u^*
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
vint = [];
for ind1 = 1:length(r)
    for ind2 = 1:length(th)
        if yplt(ind1,ind2) <= b && yplt(ind1,ind2) >= -b
            
            vint(ind1,ind2) = ustar(ind1,ind2)*f1eval(xplt(ind1,ind2), yplt(ind1,ind2))/(1 - (xplt(ind1,ind2))^2 - (yplt(ind1,ind2))^2);
           
        else    
        
            vint(ind1,ind2) = ustar(ind1,ind2)*g1eval(xplt(ind1,ind2),yplt(ind1,ind2));
            
        end
    end
    temp(ind1) = trapz(th,r(ind1)*vint(ind1,:).^2);
end
normvstar = trapz(r,temp);

% Compute L2 norm of du^*/dy using polar coordinates
temp = [];
wint = [];
for ind1 = 1:length(r)
    for ind2 = 1:length(th)
        if yplt(ind1,ind2) <= b && yplt(ind1,ind2) >= -b
            
            wint(ind1,ind2) = ustar(ind1,ind2)*f2eval(xplt(ind1,ind2), yplt(ind1,ind2))/(1 - (xplt(ind1,ind2))^2 - (yplt(ind1,ind2))^2);
           
        else    
        
            wint(ind1,ind2) = ustar(ind1,ind2)*g2eval(xplt(ind1,ind2),yplt(ind1,ind2));
            
        end
    end
    temp(ind1) = trapz(th,r(ind1)*wint(ind1,:).^2);
end
normwstar = trapz(r,temp);

upperBnd = (normvstar + normwstar)/normustar;
