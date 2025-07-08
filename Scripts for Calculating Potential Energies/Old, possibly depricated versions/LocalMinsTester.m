%% MATLAB Script: Energy Minimization and Stable Configurations
% This script finds stable configurations by minimizing the total energy 
% using fminunc, based on specified starting points and angles.

clear; clc;

%% Define Simulation Parameters
% Number of bending segments
num_segments = 2; 
EI = [5e-5, 5e-6]; % Bending stiffness (Nm^2)
L = [0.05, 0.05]; % Lengths of each bending segment (m)

% Magnetic parameters
mu0 = 4 * pi * 1e-7; % Permeability of free space (T*m/A)
M = 5.44e-3; % Dipole moment magnitude (Am^2)
m_ext = [-8.24899; 0; 1.1652e2]; % External magnet dipole moment (mx, my, mz)
start_x = 0.05; start_z = -0.11; % Start position
start_angle = 2*pi/3; % Starting angle

%% Define Initial Guess for Angles (Input from User)
initial_theta = [pi/6, -pi/6]; % Example initial guess

%% Optimize Energy Using fminunc
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
[optimal_theta, optimal_energy] = fminunc(@(theta) total_energy(theta, num_segments, EI, L, start_x, start_z, start_angle, mu0, M, m_ext), initial_theta, options);

%% Compute Energy Components at Optimized Angles
[E_elastic, E_magnetic, E_m1, E_m2] = total_energy_components(optimal_theta, num_segments, EI, L, start_x, start_z, start_angle, mu0, M, m_ext);

%% Display Results
fprintf('Start Position: (%.3f, %.3f)\n', start_x, start_z);
fprintf('Start Angle: %.3f radians\n', start_angle);
fprintf('Optimized Angles for Stable Configuration (radians):\n');
disp(optimal_theta);
fprintf('Elastic Potential Energy (J): %.6e\n', E_elastic);
fprintf('Magnetic Potential Energy of First Magnet (J): %.6f\n', E_m1);
fprintf('Magnetic Potential Energy of Second Magnet (J) (Flipped Sign): %.6f\n', E_m2);
fprintf('Total Magnetic Potential Energy (J): %.6f\n', E_magnetic);
fprintf('Total Optimized Energy (J): %.6f\n', optimal_energy);
%fprintf('Optimized Angles: %f, %f\n', optimal_theta(1), optimal_theta(2));



%% Visualization of Final Configuration
figure; hold on;
colors = ['r', 'b'];
x0 = start_x; z0 = start_z;
angle_accumulated = start_angle;
t = linspace(0, 1, 50);

for i = 1:num_segments
    R = L(i) / optimal_theta(i);
    x = x0 + R * (sin(angle_accumulated + optimal_theta(i) * t) - sin(angle_accumulated));
    z = z0 + R * (cos(angle_accumulated) - cos(angle_accumulated + optimal_theta(i) * t));
    plot(x, z, colors(i), 'LineWidth', 2);
    x0 = x(end); z0 = z(end);
    angle_accumulated = angle_accumulated + optimal_theta(i);
end

%% Add Magnetic Field Vectors
[X, Z] = meshgrid(linspace(start_x - 0.1, start_x + 0.1, 10), linspace(start_z - 0.1, start_z + 0.1, 10));
U = zeros(size(X));
V = zeros(size(Z));

for i = 1:numel(X)
    B = compute_dipole_field(m_ext, [X(i); 0; Z(i)], mu0);
    magB = norm([B(1), B(3)]);
    if magB > 0
        U(i) = B(1) / magB * 0.02;
        V(i) = B(3) / magB * 0.02;
    end
end

quiver(X, Z, U, V, 0, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);

%% ADD RED RECTANGULAR OUTLINE WITH CENTER (0,0), WIDTH=0.05, HEIGHT=0.06
wRect = 0.05;
hRect = 0.06;
rectangle('Position', [-wRect/2, -hRect/2, wRect, hRect], ...
          'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'none');

grid on;
axis equal;
title('Optimized Robot Shape in 2D');
xlabel('X Position (m)');
ylabel('Z Position (m)');
legend('Segment 1', 'Segment 2','Field Vectors');
hold off;

%% Function: Compute Total Energy
function E_total = total_energy(theta, num_segments, EI, L, start_x, start_z, start_angle, mu0, M, m_ext)
    [E_elastic, E_magnetic, ~, ~] = total_energy_components(theta, num_segments, EI, L, start_x, start_z, start_angle, mu0, M, m_ext);
    E_total = E_elastic + E_magnetic;
end

%% Function: Compute Total Energy Components
function [E_elastic, E_magnetic, E_m1, E_m2] = total_energy_components(theta, num_segments, EI, L, start_x, start_z, start_angle, mu0, M, m_ext)
    E_elastic = compute_elastic_energy(num_segments, EI, L, theta);
    [p1, p2, end_angle1, end_angle2] = compute_magnet_positions(L, start_angle, theta, start_x, start_z);
    B1 = compute_dipole_field(m_ext, p1, mu0);
    B2 = compute_dipole_field(m_ext, p2, mu0);
    F1 = [cos(end_angle1); 0; sin(end_angle1)];
    F2 = -[cos(end_angle2); 0; sin(end_angle2)];
    m1 = F1 * M;
    m2 = F2 * M;
    E_m1 = -dot(m1, B1);
    E_m2 = -dot(m2, B2);
    E_magnetic = E_m1 + E_m2;
     fprintf('Elastic Energy: %f\n', E_elastic);
     fprintf('Optimized Angles: %f, %f\n', theta);

end

%% Function: Compute Elastic Energy
function total_energy = compute_elastic_energy(num_segments, EI, L, theta)
    elastic_energy = zeros(1, num_segments);
    for i = 1:num_segments
        if theta(i) == 0  % Handle zero bending case
            elastic_energy(i) = 0;
        else
            kappa = theta(i) / L(i);
            elastic_energy(i) = 0.5 * EI(i) * (kappa^2) * L(i);
        end
    end
    total_energy = sum(elastic_energy);
   
end

%% Function: Compute Dipole Field
function B = compute_dipole_field(m, p, mu0)
    r = norm(p);
    p_hat = p / r;
    B = (mu0 / (4 * pi)) * (3 * dot(m, p_hat) * p_hat - m) / r^3;
end

%% Function: Compute Magnet Positions
function [p1, p2, end_angle1, end_angle2] = compute_magnet_positions(L, start_angle, theta, start_x, start_z)
    x0 = start_x; 
    z0 = start_z; 
    angle_accumulated = start_angle;

    % Compute first segment
    if theta(1) == 0
        x1 = x0 + L(1) * cos(angle_accumulated);
        z1 = z0 + L(1) * sin(angle_accumulated);
    else
        R1 = L(1) / theta(1);
        x1 = x0 + R1 * (sin(angle_accumulated + theta(1)) - sin(angle_accumulated));
        z1 = z0 + R1 * (cos(angle_accumulated) - cos(angle_accumulated + theta(1)));
    end
    p1 = [x1; 0; z1];
    end_angle1 = angle_accumulated + theta(1);
    angle_accumulated = end_angle1;

    % Compute second segment
    if theta(2) == 0
        x2 = x1 + L(2) * cos(angle_accumulated);
        z2 = z1 + L(2) * sin(angle_accumulated);
    else
        R2 = L(2) / theta(2);
        x2 = x1 + R2 * (sin(angle_accumulated + theta(2)) - sin(angle_accumulated));
        z2 = z1 + R2 * (cos(angle_accumulated) - cos(angle_accumulated + theta(2)));
    end
    p2 = [x2; 0; z2];
    end_angle2 = angle_accumulated + theta(2);
end

