%% MATLAB Script: Elastic Potential Energy for Piecewise Curvature (2D)
% This script calculates the total elastic potential energy of a
% continuum robot in 2D with multiple bending segments, assuming each
% segment follows a constant curvature (C-shape, S-shape, J-shape, etc.).

clear; clc;

%% Define Segment Properties
% Number of bending segments
num_segments = 2; % Change to 1 for single C-shape, 2 for S-shape, etc.

% Define bending stiffness (EI) for each segment
EI = [5e-3, 5e-3]; % Adjust based on material properties (Nm^2)

% Define lengths of each bending segment (m)
L = [0.05, 0.05]; % Adjust for each segment

% Define total bending angles for each segment (radians)
theta = [pi/4, -pi/2]; % Example: S-shape with opposite curvatures

%% Compute Elastic Potential Energy for Each Segment
elastic_energy = zeros(1, num_segments);
for i = 1:num_segments
    kappa = theta(i) / L(i); % Compute curvature
    elastic_energy(i) = 0.5 * EI(i) * (kappa^2) * L(i); % Compute elastic energy
end

%% Compute Total Elastic Potential Energy
total_energy = sum(elastic_energy);

%% Display Results
fprintf('Elastic Potential Energy per Segment (J):\n');
disp(elastic_energy);
fprintf('Total Elastic Potential Energy (J): %.6f\n', total_energy);

%% Visualization (Ensuring Continuity Between Segments)
figure; hold on;
colors = ['r', 'b']; % Colors for segments
x0 = 0; y0 = 0; % Starting position
angle_accumulated = 0; % Track accumulated angle to ensure continuity

t = linspace(0, 1, 50); % Parameter for arc plotting
for i = 1:num_segments
    % Compute arc shape using circular arc formula
    R = L(i) / theta(i); % Radius of curvature
    x = x0 + R * (sin(angle_accumulated + theta(i) * t) - sin(angle_accumulated));
    y = y0 + R * (cos(angle_accumulated) - cos(angle_accumulated + theta(i) * t));
    plot(x, y, colors(i), 'LineWidth', 2);
    x0 = x(end);
    y0 = y(end);
    angle_accumulated = angle_accumulated + theta(i); % Ensure continuity
end

grid on;
axis equal;
title('Continuum Robot Shape in 2D');
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Segment 1', 'Segment 2');