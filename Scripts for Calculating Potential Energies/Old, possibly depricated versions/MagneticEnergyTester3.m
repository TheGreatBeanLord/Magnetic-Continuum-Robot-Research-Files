% Magnetic Potential Energy Calculation for Magnetic Continuum Robot

% Define constants
mu0 = 4 * pi * 1e-7; % Permeability of free space (T*m/A)
M = 5.44e-3; % Dipole moment magnitude (Am^2)
L_tube = 0.05; % Length of each bending tube (m)

% Define dipole moment of the external magnet
m_ext = [-8.24899; 0; 1.1652e2]; % (mx, my, mz)

% Define starting position and angle of the first tube
start_x = 0.1; % User-defined starting x position
start_z = 0.05; % User-defined starting z position
start_angle = 0; % Starting angle

% Define bending angles for the tubes
theta1 = pi/6; % Angle of first tube
theta2 = pi/4; % Angle of second tube relative to first

% Compute positions of the magnets based on bending configuration
[p1, p2, end_angle1, end_angle2] = compute_magnet_positions(L_tube, start_angle, theta1, theta2, start_x, start_z);

% Compute the magnetic field from the external dipole at the small magnet positions
B1 = compute_dipole_field(m_ext, p1, mu0);
B2 = compute_dipole_field(m_ext, p2, mu0);

% Define the unit vector orientations of the small magnets based on tube end angles
F1 = [cos(end_angle1); 0; sin(end_angle1)];
F2 = -[cos(end_angle2); 0; sin(end_angle2)]; % Reverse orientation for opposite direction

% Compute the dipole moments
m1 = F1 * M;
m2 = F2 * M;

% Compute the magnetic potential energy for each dipole, switching polarity for one
E_m1 = -dot(m1, B1);
E_m2 = dot(m2, B2); % Switched polarity
E_m = E_m1 + E_m2;

% Display the results
disp(['Magnetic potential energy of first dipole: ', num2str(E_m1), ' J']);
disp(['Magnetic potential energy of second dipole (switched polarity): ', num2str(E_m2), ' J']);
disp(['Total magnetic potential energy: ', num2str(E_m), ' J']);

% Visualize the magnet positions and field vectors
visualize_magnets(p1, p2, start_angle, theta1, theta2, m_ext, mu0, L_tube, start_x, start_z);

% Function to compute dipole field
function B = compute_dipole_field(m, p, mu0)
    r = norm(p);
    p_hat = p / r;
    B = (mu0 / (4 * pi)) * (3 * dot(m, p_hat) * p_hat - m) / r^3;
end

% Function to compute magnet positions with curvature displacement
function [p1, p2, end_angle1, end_angle2] = compute_magnet_positions(L, start_angle, theta1, theta2, start_x, start_z)
    % Initialize position and angle tracking
    x0 = start_x;
    z0 = start_z;
    angle_accumulated = start_angle;
    
    % Compute first segment
    R1 = L / theta1;
    x1 = x0 + R1 * (sin(angle_accumulated + theta1) - sin(angle_accumulated));
    z1 = z0 + R1 * (cos(angle_accumulated) - cos(angle_accumulated + theta1));
    p1 = [x1; 0; z1];
    end_angle1 = angle_accumulated + theta1;
    angle_accumulated = end_angle1;
    
    % Compute second segment
    R2 = L / theta2;
    x2 = x1 + R2 * (sin(angle_accumulated + theta2) - sin(angle_accumulated));
    z2 = z1 + R2 * (cos(angle_accumulated) - cos(angle_accumulated + theta2));
    p2 = [x2; 0; z2];
    end_angle2 = angle_accumulated + theta2;
end

% Function to visualize the magnet positions with angles and field vectors
function visualize_magnets(p1, p2, start_angle, theta1, theta2, m_ext, mu0, L_tube, start_x, start_z)
    figure;
    hold on;
    
    % Plot arc paths
    colors = ['r', 'b'];
    x0 = start_x; z0 = start_z;
    angle_accumulated = start_angle;
    t = linspace(0, 1, 50);
    
    for i = 1:2
        theta = [theta1, theta2];
        R = L_tube / theta(i);
        x = x0 + R * (sin(angle_accumulated + theta(i) * t) - sin(angle_accumulated));
        z = z0 + R * (cos(angle_accumulated) - cos(angle_accumulated + theta(i) * t));
        plot(x, z, colors(i), 'LineWidth', 2);
        x0 = x(end);
        z0 = z(end);
        angle_accumulated = angle_accumulated + theta(i);
    end
    
    % Generate and normalize field vectors
    [X, Z] = meshgrid(linspace(start_x - 0.1, start_x + 0.1, 10), linspace(start_z - 0.1, start_z + 0.1, 10));
    U = zeros(size(X));
    V = zeros(size(Z));
    C = zeros(size(Z));
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
    colormap(jet);
    colorbar;
    
    % Plot magnet positions
    plot(p1(1), p1(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(p2(1), p2(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    
    xlabel('X Position (m)');
    ylabel('Z Position (m)');
    title('Magnet Positions and Robot Shape with Field Vectors');
    legend('First Tube', 'Second Tube', 'Field Vectors');
    grid on;
    axis equal;
    hold off;
end