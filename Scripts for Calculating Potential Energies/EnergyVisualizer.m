%% MATLAB Script: 3D Energy Visualization
% This script generates a 3D plot of the total potential energy based on
% the angles of the two magnets.

clear; clc;

%% Input Parameters

%for Astar part
theta1_0 = 0;
theta2_0 = 0;
start_theta = [0, 0]; % Starting from the initial state

start_x = -0.001; start_z = -0.12; % Start position
start_angle = pi/2; % Starting angle
EI = [8.05e-4, 8.05e-4]; % Bending stiffness (Nm^2)
L = [0.05, 0.05]; % Lengths of each bending segment (m)

mu0 = 4 * pi * 1e-7; % Permeability of free space (T*m/A)
M = 5.44e-3; % Dipole moment magnitude (Am^2)
m_ext = [-8.24899; 0; 1.1652e2]; % External magnet dipole moment (mx, my, mz)

%% Define Angle Ranges for Visualization
theta1_range = linspace(-pi/2, pi/2, 100);
theta2_range = linspace(-pi/2, pi/2, 100);

[Theta1, Theta2] = meshgrid(theta1_range, theta2_range);
Energy = zeros(size(Theta1));

%% Calculate Energy for Each Angle Combination
for i = 1:length(theta1_range)
    for j = 1:length(theta2_range)
        theta = [Theta1(i, j), Theta2(i, j)];
        
        % Compute Elastic and Magnetic Potential Energy
        total_elastic_energy = compute_elastic_energy(2, EI, L, theta);
        [p1, p2, end_angle1, end_angle2] = compute_magnet_positions(L, start_angle, theta, start_x, start_z);
        B1 = compute_dipole_field(m_ext, p1, mu0);
        B2 = compute_dipole_field(m_ext, p2, mu0);
        F1 = [cos(end_angle1); 0; sin(end_angle1)];
        F2 = -[cos(end_angle2); 0; sin(end_angle2)];
        m1 = F1 * M;
        m2 = F2 * M;
        E_m1 = -dot(m1, B1);
        E_m2 = -dot(m2, B2);
        total_magnetic_energy = E_m1 + E_m2;
        
        % Total Potential Energy
        Energy(i, j) = total_elastic_energy + total_magnetic_energy;
    end
end

minEnergyPerturbation([start_x, start_z], start_angle, Theta1, Theta2, start_theta, Energy);

%% Plotting 3D Energy Surface
figure;
surf(rad2deg(Theta1), rad2deg(Theta2), Energy, 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('Magnet 1 Angle (°)');
ylabel('Magnet 2 Angle (°)');
zlabel('Total Potential Energy (J)');
title('3D Plot of Total Potential Energy vs. Magnet Angles');
grid on; view(45, 30);

%% Functions (Same as in the Provided Script)
function total_energy = compute_elastic_energy(num_segments, EI, L, theta)
    elastic_energy = zeros(1, num_segments);
    for i = 1:num_segments
        if theta(i) == 0
            elastic_energy(i) = 0;
        else
            kappa = theta(i) / L(i);
            elastic_energy(i) = 0.5 * EI(i) * (kappa^2) * L(i);
        end
    end
    total_energy = sum(elastic_energy);
end

function B = compute_dipole_field(m, p, mu0)
    r = norm(p);
    p_hat = p / r;
    B = (mu0 / (4 * pi)) * (3 * dot(m, p_hat) * p_hat - m) / r^3;
end

function [p1, p2, end_angle1, end_angle2] = compute_magnet_positions(L, start_angle, theta, start_x, start_z)
    x0 = start_x; 
    z0 = start_z; 
    angle_accumulated = start_angle;

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





% Function to calculate the minimum energy perturbation required to shift to a different stable configuration
function minEnergyPerturbation(start_Pos, start_angle, Theta1, Theta2, original_Thetas, energies)
    
    % Preallocate Stable Configurations
    stableConfigs = zeros([20,20, 2]);
    
    % Evaluate Stability for Each Configuration
for i = 1:size(Theta1, 1)
    for j = 1:size(Theta1, 2)
        config = [Theta1(i, j), Theta2(i, j)];
        temp = CallableMinsTester(start_angle, start_Pos, config(1), config(2), false, 0, 0);
        stableConfigs(i, j, 1) = temp(1);
        stableConfigs(i, j, 2) = temp(2);
    end
end
    
    % Identify Starting Configuration Stability
    startStableConfig = CallableMinsTester(start_angle, start_Pos, original_Thetas(1), original_Thetas(2), false, 0, 0);
    
    % Tolerance for comparing stable configurations
    sameTolerance = 0.4;
    differentTolerance = 0.4;
    
    % Find Border Points
    isBorder = false(size(Theta1));
    for i = 2:size(Theta1, 1)-1
        for j = 2:size(Theta2, 2)-1
            % Check if the current point has the same stable config as the start
            if all(abs(squeeze(stableConfigs(i, j, :)) - startStableConfig(:)) < sameTolerance)
                % Check neighboring points
                neighbors = [squeeze(stableConfigs(i-1, j, :)); 
                             squeeze(stableConfigs(i+1, j, :));
                             squeeze(stableConfigs(i, j-1, :));
                             squeeze(stableConfigs(i, j+1, :))];
                      
                
                % If any neighbor has a sufficiently different stable config, mark as border
                for k = 1:size(neighbors, 1)
                   % fprintf("hi");
                   % disp(abs(neighbors(k, :) - startStableConfig) >= tolerance);
                    if any(all(abs(neighbors(k, :) - startStableConfig(:)) > differentTolerance))
                       isBorder(i, j) = true;
                        break;
                    end
                end
            end
        end
    end

 
    
    % Extract Border Energies and Interpolate
    borderEnergies = energies(isBorder);
    borderTheta1 = Theta1(isBorder);
    borderTheta2 = Theta2(isBorder);

    if isempty(borderTheta1) || isempty(borderTheta2)
    error('No border points found. Please check the stability logic or input configurations.');
    end
    
    % Perform Interpolation for Refinement
    resolution = max(size(Theta1));
    refinedTheta1 = linspace(min(borderTheta1), max(borderTheta1), 2*resolution);
    refinedTheta2 = linspace(min(borderTheta2), max(borderTheta2), 2*resolution);
    [RefinedTheta1, RefinedTheta2] = meshgrid(refinedTheta1, refinedTheta2);
    refinedEnergies = griddata(borderTheta1, borderTheta2, borderEnergies, RefinedTheta1, RefinedTheta2, 'cubic');
    
    % Find Minimum Energy Perturbation
    [minEnergy, minIndex] = min(refinedEnergies(:));
    minTheta1 = RefinedTheta1(minIndex);
    minTheta2 = RefinedTheta2(minIndex);

    borderTheta1 = Theta1(isBorder);
    borderTheta2 = Theta2(isBorder);
    borderEnergies = energies(isBorder);

    figure;
    surf(Theta1, Theta2, energies);
    hold on;
    scatter3(borderTheta1, borderTheta2, borderEnergies, 50, 'red', 'filled');
       legend('Bifurcation Point');
    xlabel('Theta1'); ylabel('Theta2'); zlabel('Energy');
    title('Energy Landscape with Minimum Perturbation');
    hold off;
    
    % Display Results
    fprintf('Minimum Energy Perturbation Found: %f\n', minEnergy);
    fprintf('Configuration: Theta1 = %f, Theta2 = %f\n', minTheta1, minTheta2);
    
    % Visualize Results
    figure;
    surf(Theta1, Theta2, energies);
    hold on;
    
    plot3(minTheta1, minTheta2, minEnergy, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Theta1'); ylabel('Theta2'); zlabel('Energy');
    title('Energy Landscape with Minimum Perturbation');
    hold off;
     
  
     %plotBorderPoints(Theta1, Theta2, isBorder);
    

end


function plotBorderPoints(Theta1, Theta2, isBorder)
    % Extract border point coordinates
    borderTheta1 = Theta1(isBorder);
    borderTheta2 = Theta2(isBorder);
    
    % Plot the border points
    figure;
    scatter(borderTheta1, borderTheta2, 50, 'red', 'filled');
    xlabel('Theta1 (rad)');
    ylabel('Theta2 (rad)');
    
 
    title('Border Points Visualization');
    grid on;
end

