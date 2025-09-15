clear
clc
close all

R = 6.7;
Rm = 2.4;
theta = deg2rad(14.5);

[RTi,tau1] = inner_tol(R,Rm,theta);
[RTo,tau2] = outer_tol(R,Rm,theta);








function [RT_O,tau_sol] = inner_tol(R, Rm, theta)
    % Define constants
    r1 = 2* R * sin(theta/2);
    r2 = @(tau) 2 * Rm * sin(tau);

    b1 = 2*R*sin(theta/2).^2;
    b2 = @(tau) 2*Rm*sin(tau).^2;

    % Define f and g
    f = @(tau) ( b2(tau) - b1 ) ./ cos(theta);
    g = @(tau) sqrt( r2(tau).^2 + r1^2 - 2 * r2(tau).*r1 .* cos(tau - (theta/2))  );

    % Root function
    h = @(tau) f(tau) - g(tau);

    % Search interval for tau
    tau_range = linspace(0, pi/2, 500);  % can adjust interval if needed
    h_vals = h(tau_range);

    % Find sign changes -> possible roots
    idx = find(h_vals(1:end-1).*h_vals(2:end) <= 0);

    if isempty(idx)
        % No solution found, return R
        RT_O = R;
        tau_sol = pi/2;
    else
        % Solve numerically from first sign change (smallest root)
        tau_guess = tau_range(idx(1));
        tau_sol = fzero(h, tau_guess);
        RT_O = f(tau_sol);
    end
end


function [RT_O,tau_sol] = outer_tol(R, Rm, theta)
    % Evaluating Constants
    x1 = 2*R*sin(theta/2)^2; 
    x2 = @(tau) 2*Rm*cos(tau).^2;
    
    % Root function
    h = @(tau) 2*Rm*( ( sin(tau).^2 + (R/(2*Rm)) ).*tan(theta) - sin(tau).*cos(tau));

    % Search interval for tau
    tau_range = linspace(0, pi/2, 500);  % can adjust interval if needed
    h_vals = h(tau_range);

    
    % Find sign changes -> possible roots
    idx = find(h_vals(1:end-1).*h_vals(2:end) <= 0);

    if isempty(idx)
        % No solution found, return R
        RT_O = R;
        tau_sol = pi/2;
    else
        % Solve numerically from first sign change (smallest root)
        tau_guess = tau_range(idx(1));
        tau_sol = fzero(h, tau_guess);
        RT_O =  (x1 + x2(tau_sol)) / cos(theta);
       
    end
end