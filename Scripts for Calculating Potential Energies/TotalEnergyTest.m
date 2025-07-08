%% MATLAB Script: Total Energy Calculation and Visualization
% This script calculates both the elastic and magnetic potential energy 
% of a continuum robot in 2D and visualizes the setup including magnetic field vectors.

clear; clc;

%% Define Simulation Parameters
% Number of bending segments
num_segments = 2; 
EI = [5e-3, 5e-3]; % Bending stiffness (Nm^2)
L = [0.05, 0.05]; % Lengths of each bending segment (m)
theta = [pi/4, -pi]; % Bending angles for each segment (radians)

% Magnetic parameters
mu0 = 4 * pi * 1e-7; % Permeability of free space (T*m/A)
M = 5.44e-3; % Dipole moment magnitude (Am^2)
m_ext = [-8.24899; 0; 1.1652e2]; % External magnet dipole moment (mx, my, mz)
start_x = 0.001; start_z = -0.11; % Start position
start_angle = 0; % Starting angle

%% Compute Elastic Potential Energy
total_elastic_energy = compute_elastic_energy(num_segments, EI, L, theta);

%% Compute Magnetic Potential Energy
[p1, p2, end_angle1, end_angle2] = compute_magnet_positions(L, start_angle, theta, start_x, start_z);
B1 = compute_dipole_field(m_ext, p1, mu0);
B2 = compute_dipole_field(m_ext, p2, mu0);
F1 = [cos(end_angle1); 0; sin(end_angle1)];
F2 = -[cos(end_angle2); 0; sin(end_angle2)];
m1 = F1 * M;
m2 = F2 * M;
E_m1 = -dot(m1, B1);
E_m2 = -dot(m2, B2); % Flip sign for the opposite magnet
total_magnetic_energy = E_m1 + E_m2;

%% Display Results
fprintf('Total Elastic Potential Energy (J): %.6f\n', total_elastic_energy);
fprintf('Magnetic Potential Energy of First Magnet (J): %.6f\n', E_m1);
fprintf('Magnetic Potential Energy of Second Magnet (J) (Flipped Sign): %.6f\n', E_m2);
fprintf('Total Magnetic Potential Energy (J): %.6f\n', total_magnetic_energy);

%% Visualization
visualize_setup(p1, p2, start_angle, theta, m_ext, mu0, L, start_x, start_z);

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

%% Function: Visualization
function visualize_setup(p1, p2, start_angle, theta, m_ext, mu0, L, start_x, start_z)
    figure; hold on;
    colors = ['r', 'b']; 
    x0 = start_x; 
    z0 = start_z; 
    angle_accumulated = start_angle;
    t = linspace(0, 1, 50);

    for i = 1:2
        if theta(i) == 0
            x = [x0, x0 + L(i) * cos(angle_accumulated)];
            z = [z0, z0 + L(i) * sin(angle_accumulated)];
        else
            R = L(i) / theta(i);
            x = x0 + R * (sin(angle_accumulated + theta(i) * t) - sin(angle_accumulated));
            z = z0 + R * (cos(angle_accumulated) - cos(angle_accumulated + theta(i) * t));
        end
        plot(x, z, colors(i), 'LineWidth', 2);
        x0 = x(end); 
        z0 = z(end); 
        angle_accumulated = angle_accumulated + theta(i);
    end

    [X, Z] = meshgrid(linspace(start_x - 0.1, start_x + 0.1, 10), linspace(start_z - 0.1, start_z + 0.1, 10));
    U = zeros(size(X)); V = zeros(size(Z)); C = zeros(size(Z));

    for i = 1:numel(X)
        B = compute_dipole_field(m_ext, [X(i); 0; Z(i)], mu0);
        magB = norm([B(1), B(3)]);
        if magB > 0
            U(i) = B(1) / magB * 0.02;
            V(i) = B(3) / magB * 0.02;
        end
        C(i) = magB;
    end

    C = C / max(C(:));
    quiver(X, Z, U, V, 0, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
    colormap(jet); colorbar;
    plot(p1(1), p1(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(p2(1), p2(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    xlabel('X Position (m)'); ylabel('Z Position (m)');
    title('Magnet Positions and Robot Shape with Field Vectors');
    legend('First Tube', 'Second Tube', 'Field Vectors');
    grid on; axis equal; hold off;
end
