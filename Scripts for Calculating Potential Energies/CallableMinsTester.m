function optimal_theta = CallableMinsTester(start_angle, start_pos, theta1_0, theta2_0, do_plot, iteration, max_iterations)
    % OPTIMIZE_MAGNET_ANGLES - Optimize and optionally visualize robot configuration.
    % Inputs:
    %   start_angle - initial angle (radians)
    %   start_pos   - [x, z] start position
    %   theta1_0, theta2_0 - initial guesses for magnet angles
    %   do_plot     - boolean (true/false) to enable plotting
    % Output:
    %   optimal_theta - optimized relative angles of the two magnets

    %% Parameters
    num_segments = 2;
    EI = [5.05e-5, 5.05e-5];
    L = [0.04, 0.04];
    mu0 = 4 * pi * 1e-7;
    M = 5.44e-3;
    m_ext = [-8.24899; 0; 1.1652e2];

    %% Optimization
    initial_theta = [theta1_0, theta2_0];
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'none');
    [optimal_theta, ~] = fminunc(@(theta) total_energy(theta), initial_theta, options);

    %% Visualization (if requested)
    if do_plot
       if do_plot
    drawMovement(start_angle, start_pos, optimal_theta, iteration, max_iterations);
       end
    end

    %% Nested Energy + Kinematics Functions
    function E_total = total_energy(theta)
        [E_elastic, E_magnetic] = energy_components(theta);
        E_total = E_elastic + E_magnetic;
    end

    function [E_elastic, E_magnetic] = energy_components(theta)
        E_elastic = sum(0.5 .* EI .* (theta.^2 ./ L));
        [p1, p2, a1, a2] = magnet_positions(theta);
        B1 = compute_dipole_field(m_ext, p1, mu0);
        B2 = compute_dipole_field(m_ext, p2, mu0);
        m1 = [cos(a1); 0; sin(a1)] * M;
        m2 = -[cos(a2); 0; sin(a2)] * M;
        E_magnetic = -dot(m1, B1) - dot(m2, B2);
    end

    function [p1, p2, a1, a2] = magnet_positions(theta)
        x0 = start_pos(1); z0 = start_pos(2); angle = start_angle;
        if theta(1) == 0
            x1 = x0 + L(1)*cos(angle); z1 = z0 + L(1)*sin(angle);
        else
            R1 = L(1)/theta(1);
            x1 = x0 + R1 * (sin(angle + theta(1)) - sin(angle));
            z1 = z0 + R1 * (cos(angle) - cos(angle + theta(1)));
        end
        p1 = [x1; 0; z1];
        a1 = angle + theta(1);

        if theta(2) == 0
            x2 = x1 + L(2)*cos(a1); z2 = z1 + L(2)*sin(a1);
        else
            R2 = L(2)/theta(2);
            x2 = x1 + R2 * (sin(a1 + theta(2)) - sin(a1));
            z2 = z1 + R2 * (cos(a1) - cos(a1 + theta(2)));
        end
        p2 = [x2; 0; z2];
        a2 = a1 + theta(2);
    end

    function B = compute_dipole_field(m, p, mu0)
        r = norm(p); p_hat = p / r;
        B = (mu0 / (4 * pi)) * (3 * dot(m, p_hat) * p_hat - m) / r^3;
    end

    function drawMovement(start_angle, start_pos, optimal_theta, iteration, max_iterations)
    % VISUALIZE_ROBOT_SHAPE - Plots the robot's shape and magnetic field.
    % Inputs:
    %   start_angle - initial angle (radians)
    %   start_pos   - [x, z] start position
    %   optimal_theta - optimized joint angles
    %   iteration   - current iteration number
    %   max_iterations - total number of iterations

    %% Parameters (should match those in the main script)
    num_segments = 2;
    L = [0.04, 0.04];
    mu0 = 4 * pi * 1e-7;
    m_ext = [-8.24899; 0; 1.1652e2];

    %% Plotting
    figure(1); hold on;

    % Set colors
   
    
        colors = ['r', 'b'];  % original
    

    % Draw robot
    x0 = start_pos(1); z0 = start_pos(2);
    angle_accum = start_angle;
    t = linspace(0, 1, 50);

    for i = 1:num_segments
        R = L(i) / optimal_theta(i);
        x = x0 + R * (sin(angle_accum + optimal_theta(i) * t) - sin(angle_accum));
        z = z0 + R * (cos(angle_accum) - cos(angle_accum + optimal_theta(i) * t));
        if ischar(colors)
            plot(x, z, colors(i), 'LineWidth', 2);
        else
            plot(x, z, 'Color', colors, 'LineWidth', 2);
        end
        x0 = x(end); z0 = z(end);
        angle_accum = angle_accum + optimal_theta(i);
    end

    % Plot magnetic field vectors once
    persistent field_plotted;
    if isempty(field_plotted)
        [X, Z] = meshgrid(linspace(-0.3, 0.3, 30), linspace(-.3, 0.3, 30));
        U = zeros(size(X)); V = zeros(size(Z));
        for i = 1:numel(X)
                if abs(X(i)) < 0.025 && abs(Z(i)) < 0.025
                   continue;  % Skip points within 0.025 on both axes
                 end
            B = compute_dipole_field(m_ext, [X(i); 0; Z(i)], mu0);
            magB = norm([B(1), B(3)]);
            if magB > 0
                U(i) = B(1) / magB * 0.02;
                V(i) = B(3) / magB * 0.02;
            end
        end
        quiver(X, Z, U, V, 0, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
        rectangle('Position', [-0.025, -0.03, 0.05, 0.06], ...
                  'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'none');
        grid on; axis equal;
        title('Multi-Run Robot Shape Visualization');
        xlabel('X Position (m)');
        ylabel('Z Position (m)');
        field_plotted = true;
    end
end

end


