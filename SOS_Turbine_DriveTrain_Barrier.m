%% Barrier Tester on Turbine System! (using Sondergaard et al- 2012)
clear all
close all

%Move on to 5 state system! 
%% Drive Train Sub system
%% States - 
% omega_r - Rotor angular velocity [rad/s]
% omega_g - Generator angular velocity [rad/s]
% feta_delta - Drive train torsional angle [rad]
% feta_beta - Blade-pitch [deg]
% omega_beta - Blade-pitch rate [deg/s]

%% Constants
J_r = 19.38e6; % Rotor inertia [kg m2]
B_r = 150e3; % Rotor friction [Nm/(rad/s)]
K_a = 15e6; % Drive train total spring constant [Nm/rad]
B_a = 12.15e6; % Drive train total torsional friction [Nm/(rad/s)]
N = 1/97; % Gear ratio
J_g = 534.12; % Generator inertia [kg m2]
B_g = 0; % Generator friction [Nm/(rad/s)] (Frictionless)
tau_aero = 0; %Ignore Wind Disturbance for now! (otherwise look up in table)
% Also ignoring tower top motion and lead lag bending for simplicity
% Values from paper are designed to resemble NREL 5-MW turbine from FAST

lambda_r = 2.025; % Max rotor velocity [rad/s]
lambda_delta = 441.42e-3; % Ultimate load limit of drive train torsion [rad]

%% MAKE SURE SOSTOOLs & SeDuMi (AND YALMIP?) ARE IN PATH

options.solver='SeDuMi';

syms x1 x2 x3 x4 x5;

x = [x1;x2;x3;x4;x5];
%Create SOS program
prog = sosprogram(x);

% Define set of monomials for input to sos poly variables - INCLUDE 0th
% Order!!!
VEC_4 = monomials(x,0:4); % B up to order 4!
VEC_2 = monomials(x,0:2); %Sigmas up to order 2

%Define initial and unsafe spaces
% Initial state
%g_x0 = x'*Q*x;
[prog, g_x0] = sospolyvar(prog,monomials(x,2),'wscoeff');

%Unsafe Space
g_xu_1 = -(x1 - lambda_r); % <0
g_xu_2 = -(x3^2 - lambda_delta^2);


%Construct the vector field (x' = f(x))
f = [(1/J_r)*(tau_aero - B_r*x1 - K_a*x3 - B_a*(x1 - N*x2));
    (1/J_g)*(K_a*N*x3 + B_a*N*(x1 - N*x2) - B_g*x2); 
    x1 - N*x2; 
    x5;
    -0.6*x5 - 0.0894*x4];


%Create Barrier Function B, and sigma_1 and sigma_2
[prog, B] = sospolyvar(prog,VEC_4);

[prog,sig_1] = sossosvar(prog,VEC_2);
[prog,sig_2] = sossosvar(prog,VEC_2);

%Constrain so B is positive/negative in respective regions and gradB is negative 
% Inequalities are >= 0 
prog = sosineq(prog,-B -0.1 + (g_x0-1));
prog = sosineq(prog, B -0.1 + sig_1*g_xu_1);
prog = sosineq(prog, B -0.1 + sig_2*g_xu_2);

prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2) + diff(B,x3)*f(3) + diff(B,x4)*f(4) + diff(B,x5)*f(5)));

% %Define Objective Function
prog = sossetobj(prog,coeff_1+coeff_3+coeff_6+coeff_10+coeff_15); %Minimise Tr(Q) = q1 + q3

%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,B)
g_x = sosgetsol(prog,g_x0)

%Plot Solution
xx1 = -2:0.05:3;
xx2 = -500:10:500;
xx3 = -30:0.6:30;
[x1,x2] = meshgrid(xx1,xx2);

% Set x4 and x5 to set values to plot B clearly!
% x4 = 10;
% x5 = 0;
% x3 = 0;

B = - 3.99e-6*x1^4 - 2.669e-5*x1^3*x2 - 2.323e-5*x1^3*x3 - 1.327e-9*x1^3*x4 - 2.678e-11*x1^3*x5 - 0.0001425*x1^3 - 2.446e-5*x1^2*x2^2 - 6.053e-5*x1^2*x2*x3 - 8.194e-11*x1^2*x2*x4 + 6.918e-11*x1^2*x2*x5 + 3.677e-5*x1^2*x2 - 7.198e-5*x1^2*x3^2 + 1.129e-9*x1^2*x3*x4 - 1.298e-9*x1^2*x3*x5 + 0.0001955*x1^2*x3 - 0.006393*x1^2*x4^2 - 0.002274*x1^2*x4*x5 + 7.339e-9*x1^2*x4 - 0.0068*x1^2*x5^2 + 2.12e-8*x1^2*x5 + 8.591*x1^2 - 1.455e-5*x1*x2^3 - 4.656e-5*x1*x2^2*x3 - 2.062e-10*x1*x2^2*x4 + 6.44e-11*x1*x2^2*x5 + 3.305e-5*x1*x2^2 - 6.181e-5*x1*x2*x3^2 - 3.956e-10*x1*x2*x3*x4 + 1.704e-10*x1*x2*x3*x5 + 2.87e-6*x1*x2*x3 - 0.0001249*x1*x2*x4^2 - 4.374e-5*x1*x2*x4*x5 + 3.595e-10*x1*x2*x4 - 0.0001711*x1*x2*x5^2 + 4.418e-13*x1*x2*x5 + 0.01497*x1*x2 - 9.223e-5*x1*x3^3 + 5.123e-9*x1*x3^2*x4 - 3.362e-10*x1*x3^2*x5 + 0.0003412*x1*x3^2 - 0.006026*x1*x3*x4^2 + 0.002485*x1*x3*x4*x5 - 8.501e-9*x1*x3*x4 - 0.006914*x1*x3*x5^2 - 1.146e-8*x1*x3*x5 + 1.595*x1*x3 + 3.476e-9*x1*x4^3 + 2.163e-10*x1*x4^2*x5 + 0.04681*x1*x4^2 - 3.945e-9*x1*x4*x5^2 + 0.007792*x1*x4*x5 + 2.048e-7*x1*x4 - 7.323e-10*x1*x5^3 + 0.0871*x1*x5^2 + 3.395e-7*x1*x5 + 0.001009*x1 + 1.003e-7*x2^4 - 1.799e-5*x2^3*x3 + 3.462e-12*x2^3*x4 + 1.343e-12*x2^3*x5 - 1.863e-7*x2^3 - 2.033e-5*x2^2*x3^2 - 2.542e-10*x2^2*x3*x4 + 7.177e-11*x2^2*x3*x5 - 1.775e-5*x2^2*x3 - 2.682e-9*x2^2*x4^2 + 4.003e-7*x2^2*x4*x5 - 3.181e-11*x2^2*x4 - 3.623e-7*x2^2*x5^2 + 6.969e-11*x2^2*x5 + 0.0001789*x2^2 - 3.357e-5*x2*x3^3 - 4.114e-10*x2*x3^2*x4 + 1.163e-10*x2*x3^2*x5 - 1.022e-5*x2*x3^2 - 0.0001422*x2*x3*x4^2 - 4.876e-5*x2*x3*x4*x5 - 3.27e-10*x2*x3*x4 - 0.0001855*x2*x3*x5^2 - 6.039e-11*x2*x3*x5 - 0.01244*x2*x3 + 8.071e-12*x2*x4^3 - 2.153e-11*x2*x4^2*x5 + 0.0001193*x2*x4^2 - 2.044e-11*x2*x4*x5^2 - 0.000101*x2*x4*x5 - 1.22e-10*x2*x4 + 2.523e-12*x2*x5^3 + 0.0004753*x2*x5^2 + 6.378e-10*x2*x5 + 3.376e-6*x2 - 3.154e-5*x3^4 + 1.369e-9*x3^3*x4 - 5.38e-10*x3^3*x5 + 0.0002092*x3^3 - 0.002285*x3^2*x4^2 + 0.0008615*x3^2*x4*x5 - 5.881e-9*x3^2*x4 - 0.002881*x3^2*x5^2 - 8.051e-9*x3^2*x5 + 9.085*x3^2 + 2.783e-9*x3*x4^3 + 2.114e-9*x3*x4^2*x5 + 0.009463*x3*x4^2 - 7.792e-10*x3*x4*x5^2 - 0.03435*x3*x4*x5 + 3.55e-7*x3*x4 - 2.226e-10*x3*x5^3 + 0.03882*x3*x5^2 + 1.419e-7*x3*x5 - 0.0003923*x3 - 4.093e-6*x4^4 - 2.988e-5*x4^3*x5 + 6.185e-9*x4^3 - 3.718e-5*x4^2*x5^2 + 1.871e-8*x4^2*x5 + 2.84*x4^2 - 2.224e-5*x4*x5^3 + 4.32e-8*x4*x5^2 + 4.243*x4*x5 + 9.796e-10*x4 + 1.642e-6*x5^4 + 9.759e-9*x5^3 + 7.686*x5^2 + 1.439e-9*x5 - 0.4334;
figure
mesh(x1,x2,B)
hold on

lengt = 101;
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,x2,x1*500,CO)
hold on
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt),CO)

xlabel('omega_r')
ylabel('omega_g')
x4 = 10
x5 = 0
g_x0 = @(x1,x2,x3) -1 + 0.3159*x1.^2 - 0.02004*x1.*x2 + 0.07848*x1.*x3 + 1.657e-8*x1.*x4 + 1.497e-8*x1.*x5 + 0.00122*x2.^2 - 0.003586*x2.*x3 - 1.125e-9*x2.*x4 - 1.178e-9*x2.*x5 + 6.17*x3.^2 - 1.673e-9*x3.*x4 - 1.374e-9*x3.*x5 + 0.0008966*x4.^2 + 0.00079*x4.*x5 + 0.002201*x5.^2;
fimplicit3(g_x0,[-3 3 -1000 1000 -0.5 0.5])
hold on
lengt = 101;
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,x2,x1*8,CO)
ylim([-50 50])
zlim([-0.4 0.4])


xlabel('Rotor Angular Velocity (rad/s)')
ylabel('Generator Angular Velocity (rad/s)')
zlabel('Drive Trani Torsion angle (rad)')
title('Drive Train Subsystem Operating Region (pitch = 10 deg and pitch rate =0)')
