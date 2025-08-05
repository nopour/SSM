% =========================================================================
% SSM Reduction for a 2-DOF Nonlinear Shaw-Pierre System
% =========================================================================
%
% Description:
% This script applies the equation-driven Spectral Submanifold (SSM) model 
% reduction methodology to a classic 2-DOF spring-mass-damper system with 
% a cubic nonlinearity. It performs the reduction from first principles
% using MATLAB's Symbolic Math Toolbox, without relying on the external
% SSMTool package.
%
% --- Script Workflow ---
% 1. System Definition: Construct the 4D Full Order Model (FOM).
% 2. Modal Analysis: Find the slow and fast subspaces of the linear system.
% 3. Symbolic SSM Calculation: Derive the 3rd-order manifold z = h(y).
% 4. Symbolic ROM Formulation: Derive the 2D nonlinear reduced dynamics.
% 5. Numerical Simulation: Solve both FOM and ROM ODEs.
% 6.  Visualization & Comparison: Generate journal-quality plots to

%
%% --- 0. Initial Setup and Configuration ---
clear;
clc;
close all;
format long g;
tic_script_total = tic;

% --- Plotting Style Configuration for Professional, Journal-Quality Plots ---
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', 'w');

%% --- 1. System Definition (Full Order Model - FOM) ---
disp('Part 1: Defining the 2-DOF nonlinear system (FOM)...');

% System Parameters
m1 = 1.0; m2 = 1.0;         % Masses
k1 = 2.0; k2 = 1.0; k3 = 2.0; % Linear stiffnesses
c1 = 0.02; c2 = 0.01; c3 = 0.02;% Linear damping
gamma_nl = 0.5;             % Cubic nonlinearity coefficient

% Physical System Matrices
M = [m1, 0; 0, m2];
C = [c1+c2, -c2; -c2, c2+c3];
K = [k1+k2, -k2; -k2, k2+k3];

% State-space formulation: x_dot = A*x + f_nl(x)
% State vector x = [q1; q2; q1_dot; q2_dot]
A_fom = [zeros(2), eye(2);
         -M\K,    -M\C];

% Define symbolic variables for nonlinearity
syms q1_s q2_s real
F_nl_sym = [-gamma_nl * q1_s^3; 0];
f_nl_sym = [zeros(2,1); M\F_nl_sym];

fprintf('FOM defined. System has 4 dimensions.\n\n');

%% --- 2. Modal Analysis and Subspace Selection ---
disp('Part 2: Performing modal analysis to find slow/fast subspaces...');

% Eigen-analysis of the linear part
[V_eig, D_eig] = eig(A_fom);
lambda = diag(D_eig);

% Sort eigenvalues by real part to identify slow vs. fast modes
[~, sort_idx] = sort(real(lambda), 'descend');
lambda_sorted = lambda(sort_idx);
V_sorted = V_eig(:, sort_idx);

% Construct the real modal transformation matrix V_modal
V_modal = [real(V_sorted(:,1)), imag(V_sorted(:,1)), ...
           real(V_sorted(:,3)), imag(V_sorted(:,3))];
V_inv_modal = inv(V_modal);

% Transform the linear system into modal coordinates (eta = [y; z])
A_modal = V_inv_modal * A_fom * V_modal;
Ay = A_modal(1:2, 1:2); % 2x2 slow dynamics matrix
Az = A_modal(3:4, 3:4); % 2x2 fast dynamics matrix

fprintf('<strong>Linear Slow Subsystem (ROM) Matrix Ay:</strong>\n'); disp(Ay);
fprintf('<strong>Linear Fast Subsystem Matrix Az:</strong>\n'); disp(Az);

%% --- 3. Symbolic SSM Manifold Calculation (z = h(y)) ---

% Define symbolic modal coordinates
syms y1 y2 real;
y_sym = [y1; y2];
z_sym = sym('z', [2, 1]); assume(z_sym, 'real');
eta_sym = [y_sym; z_sym];

% Express nonlinearity in modal coordinates
x_phys_from_modal = V_modal * eta_sym;
q1_from_modal = x_phys_from_modal(1);
f_nl_in_modal = V_inv_modal * subs(f_nl_sym, q1_s, q1_from_modal);

% Extract nonlinear terms for slow (fy) and fast (fz) dynamics
fy_sym = f_nl_in_modal(1:2);
fz_sym = f_nl_in_modal(3:4);

% --- Solve the Invariance (Homological) Equation ---
% We seek z = h(y) = h3(y), where h3 contains cubic terms in y1, y2.
% The invariance equation is: Dh(y)*[Ay*y + fy(y,h)] = Az*h(y) + fz(y,h)
% For 3rd order, this simplifies to: Dh3(y)*Ay*y - Az*h3(y) = fz_3(y,0)

% Define h3 with unknown coefficients
C_h3 = sym('c', [2, 4]); % 2 slave coords, 4 cubic monomials
monomials_y3 = [y1^3; y1^2*y2; y1*y2^2; y2^3];
h3_sym = C_h3 * monomials_y3;

% Jacobian of h3
Dh3_sym = jacobian(h3_sym, y_sym);

% Extract cubic part of the fast forcing fz
fz_at_z0 = subs(fz_sym, z_sym, [0;0]);
fz3 = taylor(fz_at_z0, y_sym, 'Order', 4) - taylor(fz_at_z0, y_sym, 'Order', 3);

% Form the homological equation
homological_eq = Dh3_sym * Ay * y_sym - Az * h3_sym - fz3;

% Solve for coefficients by equating monomials to zero
% This is done by taking derivatives w.r.t. y1 and y2 and evaluating at (0,0)
eqs_to_solve = sym([]);
for i = 0:3
    j = 3 - i; % Powers of y1 and y2 sum to 3
    % The coefficient of y1^i * y2^j in a polynomial P is given by
    % (1/(i!*j!)) * d^(i+j)P / (dy1^i * dy2^j) evaluated at (0,0)
    
    % Extract coefficient from the first component of the homological equation
    coeff1 = diff(homological_eq(1), y1, i);
    coeff1 = diff(coeff1, y2, j);
    eq1 = subs(coeff1, y_sym, [0;0]) == 0;
    
    % Extract coefficient from the second component
    coeff2 = diff(homological_eq(2), y1, i);
    coeff2 = diff(coeff2, y2, j);
    eq2 = subs(coeff2, y_sym, [0;0]) == 0;
    
    eqs_to_solve = [eqs_to_solve; eq1; eq2];
end

% Solve the linear system for the manifold coefficients
sol = solve(eqs_to_solve, C_h3(:));
C_h3_solved = reshape(struct2array(sol), 2, 4);
h3_final = C_h3_solved * monomials_y3;

fprintf('<strong>Symbolic 3rd-Order SSM Manifold z = h(y):</strong>\n');
fprintf('h1(y1, y2) = \n'); disp(vpa(h3_final(1), 4));
fprintf('h2(y1, y2) = \n'); disp(vpa(h3_final(2), 4));
fprintf('SSM calculation complete.\n\n');

%% --- 4. Symbolic ROM Formulation ---
% ROM equation: y_dot = Ay*y + fy(y, h(y))
fy_on_manifold = subs(fy_sym, z_sym, h3_final);
ROM_eqs_sym = Ay * y_sym + fy_on_manifold;

fprintf('<strong>Symbolic 2D ROM Equations y_dot = f_rom(y):</strong>\n');
fprintf('dy1/dt = \n'); disp(vpa(ROM_eqs_sym(1), 4));
fprintf('dy2/dt = \n'); disp(vpa(ROM_eqs_sym(2), 4));
fprintf('ROM formulation complete.\n\n');

%% --- 5. Numerical Simulation Setup ---
% Function handle for the manifold z = h(y)
% Input: a 2xN vector y_in. Output: a 2xN vector z.
h_func = matlabFunction(h3_final, 'Vars', {[y1; y2]});

% Function handle for the ROM ODEs dy/dt = f_rom(y)
% Input: t (scalar), y_in (2x1 vector). Output: 2x1 vector.
rom_ode_func = matlabFunction(ROM_eqs_sym, 'Vars', {'t', [y1; y2]});

% Function handle for the FOM ODEs dx/dt = f_fom(x)
% Input: t (scalar), x_in (4x1 vector). Output: 4x1 vector.
syms q1d_s q2d_s real; 
x_phys_sym_vec = [q1_s; q2_s; q1d_s; q2d_s];
fom_eqs_sym = A_fom * x_phys_sym_vec + subs(f_nl_sym, 'q1_s', q1_s);
fom_ode_func = matlabFunction(fom_eqs_sym, 'Vars', {'t', x_phys_sym_vec});

% Simulation parameters
t_span = [0, 150];
ode_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);

% Initial conditions ON the manifold
y0 = [1.0; 0.2]; % Initial slow coordinates
z0 = h_func(y0); % Call with a single vector
eta0 = [y0; z0]; % Initial modal coordinates
x0_fom = V_modal * eta0; % Initial physical coordinates for FOM

% --- Run Simulations ---
% FOM Simulation
tic_fom = tic;
[t_fom, x_fom] = ode45(fom_ode_func, t_span, x0_fom, ode_options);
time_fom = toc(tic_fom);
fprintf('FOM simulation complete. Elapsed time: %.4f seconds.\n', time_fom);

% ROM Simulation
tic_rom = tic;
[t_rom, y_rom] = ode45(rom_ode_func, t_span, y0, ode_options);
time_rom = toc(tic_rom);
fprintf('ROM simulation complete. Elapsed time: %.4f seconds.\n\n', time_rom);

% Reconstruct full state from ROM for comparison
z_rom_recon = h_func(y_rom'); % y_rom is Nx2, so y_rom' is 2xN
eta_rom_recon = [y_rom'; z_rom_recon];
x_rom_recon = V_modal * eta_rom_recon;
x_rom_recon = x_rom_recon';

fprintf('Simulations complete.\n\n');

%% --- 6. Visualization and Comparison ---

% --- Figure 1: 3D Manifold Visualization ---
figure('Name', '3D SSM Manifold', 'Position', [100, 100, 900, 700]);
[y1_grid, y2_grid] = meshgrid(linspace(-1.2, 1.2, 50), linspace(-1.2, 1.2, 50));
grid_points = [y1_grid(:)'; y2_grid(:)']; % Create a 2x(50*50) matrix of grid points
h_grid_vals = h_func(grid_points); % Evaluate all points at once
% Reshape the output back to the grid size for plotting
h_manifold_z1 = reshape(h_grid_vals(1,:), size(y1_grid));

s = surf(y1_grid, y2_grid, h_manifold_z1, 'FaceColor', [0.7, 0.85, 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
s.FaceLighting = 'gouraud';
light('Position',[-1 -1 1],'Style','infinite');
hold on;

% Project FOM trajectory onto this space
eta_fom = V_inv_modal * x_fom';
y_fom_traj = eta_fom(1:2, :)';
z_fom_traj = eta_fom(3:4, :)';
plot3(y_fom_traj(:,1), y_fom_traj(:,2), z_fom_traj(:,1), 'k-', 'LineWidth', 2.5, 'DisplayName', 'FOM Trajectory');

% Plot reconstructed ROM trajectory
plot3(y_rom(:,1), y_rom(:,2), z_rom_recon(1,:)', 'r--', 'LineWidth', 2.5, 'DisplayName', 'ROM Trajectory');

scatter3(y0(1), y0(2), z0(1), 150, 'g', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Start');
hold off;
xlabel('$y_1$ (Slow Modal Coord)');
ylabel('$y_2$ (Slow Modal Coord)');
zlabel('$z_1$ (Fast Modal Coord)');
title('\textbf{Nonlinear Spectral Submanifold (SSM) and System Dynamics}');
legend('show', 'Location', 'northeast');
view(35, 25);
grid on;
axis tight;
box on;
fprintf('3D Manifold plot generated.\n');

% --- Figure 2: Trajectory Comparison in Physical Coordinates ---
figure('Name', 'Physical Trajectory Comparison', 'Position', [1050, 100, 900, 700]);
sgtitle('\textbf{Comparison of FOM and ROM in Physical Coordinates}');

% q1 vs time
subplot(2, 2, 1);
plot(t_fom, x_fom(:,1), 'k-', 'DisplayName', 'FOM');
hold on;
plot(t_rom, x_rom_recon(:,1), 'r--', 'DisplayName', 'ROM');
hold off;
grid on;
xlabel('Time (s)');
ylabel('$q_1(t)$');
title('Displacement of Mass 1');
legend('show', 'Location', 'southeast');

% q2 vs time
subplot(2, 2, 2);
plot(t_fom, x_fom(:,2), 'k-', 'DisplayName', 'FOM');
hold on;
plot(t_rom, x_rom_recon(:,2), 'r--', 'DisplayName', 'ROM');
hold off;
grid on;
xlabel('Time (s)');
ylabel('$q_2(t)$');
title('Displacement of Mass 2');
legend('show', 'Location', 'southeast');

% Phase Portrait
subplot(2, 2, 3);
plot(x_fom(:,1), x_fom(:,3), 'k-', 'DisplayName', 'FOM');
hold on;
plot(x_rom_recon(:,1), x_rom_recon(:,3), 'r--', 'DisplayName', 'ROM');
scatter(x0_fom(1), x0_fom(3), 100, 'g', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Start');
hold off;
grid on;
xlabel('$q_1$');
ylabel('$\dot{q}_1$');
title('Phase Portrait for Mass 1');
legend('show', 'Location', 'northwest');

% Error Plot
subplot(2, 2, 4);
% Interpolate ROM to match FOM time steps for error calculation
x_rom_interp = interp1(t_rom, x_rom_recon, t_fom);
error_q1 = abs(x_fom(:,1) - x_rom_interp(:,1));
error_q2 = abs(x_fom(:,2) - x_rom_interp(:,2));
semilogy(t_fom, error_q1, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Error in $q_1$');
hold on;
semilogy(t_fom, error_q2, 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'Error in $q_2$');
hold off;
grid on;
xlabel('Time (s)');
ylabel('Absolute Error');
title('Error between FOM and ROM');
legend('show', 'Location', 'northeast');
fprintf('Physical coordinate comparison plots generated.\n');

% --- Figure 3: Performance Comparison ---
figure('Name', 'Computational Time Comparison', 'Position', [1050, 600, 600, 500]);
bar_data = [time_fom; time_rom];
b = bar(bar_data, 'FaceColor', 'flat');
b.CData(1,:) = [0.2 0.2 0.2]; % Dark gray for FOM
b.CData(2,:) = [1 0.4 0.4]; % Red for ROM
set(gca, 'xticklabel', {'FOM (4-DOF)', 'ROM (2-DOF)'});
ylabel('Elapsed Time for Simulation (s)');
title('\textbf{Computational Performance}');
grid on;
text(b.XEndPoints, b.YEndPoints, string(round(b.YData,4))+' s', 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 12);
fprintf('Computational time comparison plot generated.\n\n');

%% --- 7. Final ---
fprintf('Total script execution time: %.4f seconds\n', toc(tic_script_total));


