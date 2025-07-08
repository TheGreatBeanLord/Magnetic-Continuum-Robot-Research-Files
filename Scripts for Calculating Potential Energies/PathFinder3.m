
% Script: find_transition_path.m
% Description: Finds a path from a start configuration to a target configuration
% by varying the base position and orientation, calling a user-defined
% function to compute the stable magnet angles at each step.

clear CallableMinsTester;


% Parameters (set these manually for now)
start_angle = -pi/2; % radians
start_pos = [0.11; 0.13]; % [x; z]
theta_current = [0; -pi]; % initial guess for magnet angles

target_theta = [pi/4; +pi/4]; % desired magnet angles

%record original start stuff for later drawing
original_Angle = start_angle;
original_Pos = start_pos;
original_Theta = theta_current;


% Step sizes
step_x = 0.005;  % mm
step_z = 0.005;  % mm
step_angle = 0.1;  % radians

tolerance = deg2rad(30);  % radians
max_iterations = 100;
plot_interval = 2;  % plot after every N movements

% Constraints
max_offset = 0.15;  % mm


% Logging
path_log = [];
loss_log = [];
 fprintf('  Start Pos: (%.3f, %.3f), Start Angle: %.3f\n', ...
        start_pos(1), start_pos(2), start_angle);
 

for iter = 1:max_iterations
  
    % Evaluate current configuration
    gonnaPrint = false;
    if mod(iter, plot_interval) == 0
        gonnaPrint = true;
    end

    optimal_theta = CallableMinsTester(start_angle, start_pos, theta_current(1), theta_current(2), gonnaPrint, iter,  max_iterations);
    theta_current = optimal_theta;


    % Compute loss (distance to target)
    %loss(1) = sum((optimal_theta(1) - target_theta(1)).^2);
    %loss(2) = sum((optimal_theta(2) - target_theta(2)).^2);
      optimal_theta = optimal_theta(:);
       target_theta = target_theta(:);
    loss = norm((optimal_theta - target_theta));
  %  path_log = [path_log; start_pos', start_angle];
   
  %loss_log = [loss_log; loss];
    
    fprintf('Iteration %d\n', iter);
    fprintf('  Start Pos: (%.3f, %.3f), Start Angle: %.3f\n', ...
        start_pos(1), start_pos(2), (start_angle));
    fprintf('  Magnet Angles: [%.3f, %.3f]\n', ...
        (optimal_theta(1)), (optimal_theta(2)));
    fprintf('  Desired Angles: [%.3f, %.3f]\n', ...
        (target_theta(1)), (target_theta(2)));
    
    fprintf('  Loss: %.6f\n', loss);

    if loss < tolerance^2
        disp('Target reached!');
        break;
    end



    % Generate candidate perturbations
    dx = [-1, 0, 1] * step_x;
    dz = [-1, 0, 1] * step_z;
    dtheta = [-1, 0, 1] * step_angle;
    
    [DX, DZ, DTH] = ndgrid(dx, dz, dtheta);
    candidates = [DX(:), DZ(:), DTH(:)];
    candidates(all(candidates == 0, 2), :) = [];  % remove zero step

    best_loss = Inf;
    
    for i = 1:size(candidates, 1)
        new_pos = start_pos + [candidates(i,1); candidates(i,2)];
        new_angle = start_angle + candidates(i,3);

        % Check position constraints
        if any(abs(new_pos) > max_offset)
            continue;
        end
        
        theta_try = CallableMinsTester(new_angle, new_pos, theta_current(1), theta_current(2), false, iter, max_iterations);
        
          % Perform collision check
        if checkCollision(new_pos, new_angle, theta_try, [0, 0], 0.05) % Adjust the magnet position and radius
            continue;
        end
        
        theta_try = theta_try(:);
        target_theta = target_theta(:);
        candidate_loss = norm((theta_try - target_theta));
       

        if candidate_loss < best_loss
            best_loss = candidate_loss;
            best_candidate = [new_pos; new_angle];
            best_theta = theta_try;
        end

    end

    % If no valid candidate found, exit
    if best_loss == Inf
        disp('No valid candidates within bounds. Exiting.');
        break;
    end

       fprintf("Best Loss: %d, current loss: %d", best_loss, loss);

    if best_loss > loss  % no improvement
        fprintf('Stuck, trying random perturbation.\n');
        start_pos = start_pos + 0.01 * (rand(2,1)-0.5);
        start_angle = start_angle + 0.01 * (rand()-0.5);
        continue;
    
    else
    % Update to best candidate
      start_pos = best_candidate(1:2);
     start_angle = best_candidate(3);
     theta_current = best_theta;
  
    end

    % Plot every N steps
    %{
    if mod(iter, plot_interval) == 0
        figure(1); clf;
        subplot(2,1,1);
        plot(loss_log, '-o');
        xlabel('Iteration');
        ylabel('Loss (squared rad)');
        title('Loss vs. Iteration');
        grid on;

        subplot(2,1,2);
        plot(path_log(:,1), path_log(:,2), 'o-');
        xlabel('Start X (mm)');
        ylabel('Start Z (mm)');
        title('Start Position Path');
        grid on;

        drawnow;
    end
    %}
end

% Final plot
%{

figure;
plot(loss_log, '-o');
xlabel('Iteration');
ylabel('Loss (squared rad)');
title('Final Loss vs. Iteration');
grid on;
%}


quickdraw(optimal_theta, start_pos, start_angle, "Yellow");
quickdraw(original_Theta+0.0001, original_Pos, original_Angle, "k");

qw{1} = plot(nan, 'r');
qw{2} = plot(nan, 'b');
qw{3} = plot(nan, 'yellow');
qw{4} = plot(nan, 'k'); % You can add an extra element too
legend([qw{:}], {'Segment 1','Segment 2','End Configuration', 'Original Configuration'}, 'location', 'best')

function quickdraw(drawtheta, start_pos, start_angle, color)

  %% Parameters (should match those in the main script)
    num_segments = 2;
    L = [0.04, 0.04];
    mu0 = 4 * pi * 1e-7;
    m_ext = [-8.24899; 0; 1.1652e2];

    %% Plotting
    figure(1); hold on;

    % Set colors
   
        colors = color; 




    % Draw robot
    x0 = start_pos(1); z0 = start_pos(2);
    angle_accum = start_angle;
    t = linspace(0, 1, 50);

    for i = 1:num_segments
        R = L(i) / drawtheta(i);
        x = x0 + R * (sin(angle_accum + drawtheta(i) * t) - sin(angle_accum));
        z = z0 + R * (cos(angle_accum) - cos(angle_accum + drawtheta(i) * t));
        if ischar(colors)
            plot(x, z, colors(i), 'LineWidth', 2);
        else
            plot(x, z, 'Color', colors, 'LineWidth', 2);
        end
        x0 = x(end); z0 = z(end);
        angle_accum = angle_accum + drawtheta(i);
    end
end


clear CallableMinsTester;



function isCollision = checkCollision(pos, angle, drawtheta, magnetPos, magnetRadius)
    % Parameters
    num_segments = 2;
    L = [0.04, 0.04];
    
    % Initialize positions
    x0 = pos(1);
    z0 = pos(2);
    angle_accum = angle;
    isCollision = false;

    % Check each segment for collision
    for i = 1:num_segments
        R = L(i) / drawtheta(i);
        t = linspace(0, 1, 50);
        x = x0 + R * (sin(angle_accum + drawtheta(i) * t) - sin(angle_accum));
        z = z0 + R * (cos(angle_accum) - cos(angle_accum + drawtheta(i) * t));

        % Calculate distance from magnet center
        distances = sqrt((x - magnetPos(1)).^2 + (z - magnetPos(2)).^2);
        if any(distances < magnetRadius)
            isCollision = true;
            return;
        end

        x0 = x(end);
        z0 = z(end);
        angle_accum = angle_accum + drawtheta(i);
    end
end
