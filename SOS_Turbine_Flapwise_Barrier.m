%% Barrier Tester on Turbine System! (using Sondergaard et al- 2012)
clear all
close all

%Start on 2 state system! 
% Flapewise Blade Bending
% States - Flapwise - Blade Tip Displacement (x1) & Blade tip Velocity(x2) (makes
% sense)
% K_flap = Flapwise blade bending spring constant [N/m]
% B_flap = Flapwise blade bending damping constant [N/(m/s)]
% M_flap = Fictitious mass of the blade beyond the break point [kg]
% F_aero = Aerodynamic Thrust Force (N)
% Values from paper are designed to resemble NREL 5-MW turbine from FAST

K_flap = 55.25e3;
B_flap = 250.00e3;
M_flap = 66.00e3;
F_aero = 0; %Ignore Wind Disturbance for now! (otherwise look up in table)
% Also ignoring tower top motion for simplicity

% Max flapwise blade tip displacement = 11.57m
lambda_sqrd = 11.57^2;

%% MAKE SURE SOSTOOLs & SeDuMi (AND YALMIP?) ARE IN PATH

options.solver='SeDuMi';

syms x1 x2;

x = [x1;x2];

%Create SOS program
prog = sosprogram([x1;x2]);
[prog,Q] = sospolymatrixvar(prog,monomials(x,0),[2 2],'symmetric');

%Define initial ans unsafe spaces
% Initial state
%g_x0 = x'*Q*x;
[prog, g_x0] = sospolyvar(prog,monomials(x,2),'wscoeff');
g_xu = -(x1^2-lambda_sqrd); % >0

%Construct the vector field (x' = f(x))
f = [x2; (1/M_flap)*(F_aero - K_flap*x1 - B_flap*x2)];

% Define set of monomials for input to sos poly variables - INCLUDE 0th
% Order!!!
VEC_4 = monomials([x1; x2],0:4); % B up to order 4!
VEC_2 = monomials([x1; x2],0:2); %Sigmas up to order 2
VEC_1 = monomials([x1; x2],0:1); % Up to order 1

%Create Barrier Function B, and sigma_1 and sigma_2
[prog, B] = sospolyvar(prog,VEC_2);

[prog,sig_1] = sossosvar(prog,VEC_1);
[prog,sig_2] = sossosvar(prog,VEC_1);

%Constrain so B is positive/negative in respective regions and gradB is negative 
% Inequalities are >= 0 
prog = sosineq(prog,-B -0.1 + (g_x0-1));
prog = sosineq(prog, B -0.1 + sig_2*g_xu);

prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2)));

% %Define Objective Function
 prog = sossetobj(prog,coeff_4+coeff_6); %Minimise Tr(Q) = q1 + q3

%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,B)
g_x0 = sosgetsol(prog,g_x0)

%Plot Solution
xx1 = -30:0.5:30;
xx2 = -30:0.5:30;
[x1,x2] = meshgrid(xx1,xx2);

B = 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1.1;  
figure
mesh(x1,x2,B)
hold on

%Plot plane of B = 0
lengt = length(x1);
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+11.57,x2,x1/8,CO)
mesh(zeros(lengt)-11.57,x2,x1/8,CO)
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt),CO)

xlim([-15 15])
ylim([-30 30])
zlim([-2 5])
xlabel('Flapwise Tip Displacement')
ylabel('Flapwise Tip Velocity')
title('Flapwise Blade Bending Barrier Function, Green = operating region, Red = unsafe region')

% Contour Plot
%plot unsafe and operation regions
x_u_1 = nsidedpoly(4, 'Center', [11.57+30 0], 'Sidelength', 60);
x_u_2 = nsidedpoly(4, 'Center', [-11.57-30 0], 'Sidelength', 60);
g_x0 = @(x1,x2) 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1;


figure
%Plot Solution
xx1 = -30:2:30;
xx2 = -30:2:30;
[x1,x2] = meshgrid(xx1,xx2);
B = 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1.1;

[U,V2] = gradient(-B); 
%contour(x1,x2,B)
hold on
quiver(x1,x2,U,V2)
hold on
plot(x_u_1, 'FaceColor', 'r')
hold on
plot(x_u_2,'FaceColor', 'r')
hold on
fimplicit(g_x0,[-15 15 -25 25],'MeshDensity',300,'Color','g')

xlim([-15 15])
ylim([-30 30])
xlabel('Flapwise Tip Displacement')
ylabel('Flapwise Tip Velocity')
title('Flapwise Blade Bending State Space, Green = operating region, Red = unsafe region')
