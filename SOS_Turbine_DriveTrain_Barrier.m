%% Barrier Tester on Turbine System! (using Sondergaard et al- 2012)
%Model based off of OpenFAST 5MW turbine
clear all
close all
echo on
%Move on to 5 state system! 
%% Drive Train Sub system
%% States - 
% omega_r - Rotor angular velocity [rad/s]
% omega_g - Generator angular velocity [rad/s]
% feta_delta - Drive train torsional angle [rad]
% feta_beta - Blade-pitch [deg] (actually 90 minus pitch)
% omega_beta - Blade-pitch rate [deg/s] (actually minus pitch rate)

%% Constants
J_r = 19.38e6; % Rotor inertia [kg m2]
B_r = 150e3; % Rotor friction [Nm/(rad/s)]
K_a = 15e6; % Drive train total spring constant [Nm/rad]
B_a = 12.15e6; % Drive train total torsional friction [Nm/(rad/s)]
N = 1/97; % Gear ratio
J_g = 534.12; % Generator inertia [kg m2]
B_g = 0; % Generator friction [Nm/(rad/s)] (Frictionless)

%Ignore Wind Disturbance for now! (otherwise look up in table)
rho = 1.22521;
A = 11.88e3;
R = 61.5;
v_wind = 15;
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
VEC_4 = monomials(x,0:2); % B up to order 4!
VEC_2 = monomials(x,0:2); %Sigmas up to order 2

%Define initial and unsafe spaces
% Initial state
%g_x0 = x'*Q*x;
[prog, g_x0] = sospolyvar(prog,monomials(x,2),'wscoeff');
%g_x0 = x1^2/(1^2) + x2^2/(97^2) + x3^2/(0.1^2) + x4^2/(10^2) + x5^2/(0.01^2)  ; % = 1

%Unsafe Space
g_xu_1 = -x1 + lambda_r; % <0
g_xu_2 = -x3^2 +lambda_delta^2;


%Approximated C_q coefficients
C_q_Fit = [-0.0240889231402065,0.00108816365174505,0.0234530529270319,-9.14549613647725e-06,-0.000613719856304388,-0.00158790611519949];

% Simple one
C_q_Fit = [0 0 0 0 0 0];

%Linear
C_q_Fit = [0.0887009882823798, -0.000994478306871089, -0.0111240840652831, 0 0 0];


%Construct the vector field (x' = f(x))
f = [(1/J_r)*(0.5*rho*A*R*v_wind^2*(C_q_Fit(1) + C_q_Fit(2).*(90-x4) + C_q_Fit(3).*(R*x1/v_wind) + C_q_Fit(4).*(90-x4).^2 +C_q_Fit(5).*(90-x4).*(R*x1/v_wind) + C_q_Fit(6).*(R*x1/v_wind).^2)- B_r*x1 - K_a*x3 - B_a*(x1 - N*x2)) ;
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

prog = sosineq(prog, B -0.1 + sig_2*g_xu_2) ;

prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2) + diff(B,x3)*f(3) + diff(B,x4)*f(4) + diff(B,x5)*f(5)));

% %Define Objective Function
prog = sossetobj(prog,coeff_1*10 + coeff_3/100 + coeff_6*100 + coeff_10 + coeff_15/1000); %Minimise Tr(Q) = q1 + q3

%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,B)
g_x = sosgetsol(prog,g_x0)
sig_1 = sosgetsol(prog,sig_1)
sig_2 = sosgetsol(prog,sig_2)

%% Plot Barrier Function
xx1 = -2:0.05:3;
xx2 = -500:10:500;
xx3 = -0.5:0.01:0.5;
xx4 = -150:3:150;
[x1,x4] = meshgrid(xx1,xx4);

%Convert Sym variable to function handle
B = matlabFunction(SOLV);

figure
mesh(x1,x4,B(x1,97*x1,0,90-x4,0)) %plot B for pitch (90 - x2) and rotor speed, pitch rate = 0, torsion angle = 0
hold on

lengt = 101;
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,x4,x1*10,CO)
hold on
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x4,zeros(lengt)+0.1,CO)

xlim([-2 3])
ylim([-10 120])
xlabel('Rotor Speed (rad/s)')
ylabel('Initial Pitch (deg)')


%% Plot Operation Region


%Convert Sym variable to function handle
g_x0_plot = matlabFunction(g_x-1);

%Change variables so it's in pitch and pitch rate :))
% (And set some variables so we can plot in just 3D)
% Set pitch rate and torsion angle?
g_x0_plot1 = @(x1,x4,x3) g_x0_plot(x1,97*x1,x3,90-x4,0);


figure
fimplicit3(g_x0_plot1,[-2 3 -20 110 -0.5 0.5 ])

hold on
[x1,x4] = meshgrid(xx1,xx4);
lengt = 101;

CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,x4,x1*8,CO)
hold on
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x4,zeros(lengt)+0.44142,CO)
zlim([-0.5 0.5])

xlabel('Rotor Angular Velocity (rad/s)')
ylabel('Initial Pitch (deg)')
zlabel('Torsion Angle (rad)')
title('Drive Train Subsystem Operating Region')



% Q = [0.23, -4e-8, -3.1e-8, 7.38e-9 -3.80e-8;
% -4e-8 5.5e-7 -1.16e-10 -1.30e-11 1.07e-12;
% -3.1e-8 -1.16e-10 861.7 -3.41e-12 5.97e-13;
% 7.38e-9 -1.30e-11 -3.41e-12 0.56e-4 -3.23e-11;
% -3.80e-8 1.07e-12 5.97e-13 -3.23e-11 0.13e-3];
% 
% figure
% x3 = 0; % Pitch = 10 so, feta = (90 - pitch) = 80
% x5 = 0;
% % REMEMBER THE '-1'IN g_x0
% g_x0 = @(x1,x2,x4) -1 + [x1, x2, x3, x4,x5]*Q*[x1, x2, x3,x4,x5]';
% fimplicit3(g_x0,[-2.5 2.5 -2000 2000 -120 120])



