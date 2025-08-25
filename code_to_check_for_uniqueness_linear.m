%% Checking for multiplicity in the linear case
clc; clear; close all;

% Put in any parameter values
alpha = 1;     % precision of public prior
beta = 0.055;  %  precision of private signal
z = 0.6;       % mean of public prior

% Determinig gamma (equation (34))
gamma = (alpha^2*(beta+alpha))/(beta*(alpha+2*beta));

% Taking the "inverse" of the positive part of the standard normal pdf
x = linspace(0, 10, 1000); 
pdf_vals = normpdf(x); % Compute standard normal PDF

% Creating the inverse function via interpolation
inv_phi_pos = @(y) interp1(pdf_vals, x, y, 'pchip', NaN);
% y_vals = linspace(max(pdf_vals), min(pdf_vals), 1000); 
% x_inverse = arrayfun(inv_phi_pos, y_vals);

% Determining k** and k*** acoording to Proposition 4
% Simplifying some terms
sqrt_term = sqrt((alpha + 2*beta) * (alpha + beta) / (beta * alpha^2));
phi_inv_val = inv_phi_pos(1/sqrt(gamma));

% Computing k** and k***
k_star_star  = z - phi_inv_val * sqrt_term;
k_star_star_star  = z + phi_inv_val * sqrt_term;

fktvalue1 = normcdf( sqrt((beta*alpha^2)/((alpha+2*beta)*(alpha+beta))) * (k_star_star - z) )- (beta/(beta+alpha))*k_star_star - (alpha/(beta+alpha))*z;
fktvalue2 = normcdf( sqrt((beta*alpha^2)/((alpha+2*beta)*(alpha+beta))) * (k_star_star_star - z) )- (beta/(beta+alpha))*k_star_star_star - (alpha/(beta+alpha))*z;

% Applying Proposition 4 to determine the number of equilibria
if 2*pi >= gamma
    disp('A unique equilibrium')
elseif 2*pi < gamma && sign(fktvalue1) == sign(fktvalue2)
    disp('A unique equilibrium')
else
    disp('Multiple equilibria')
end



% Plottinh the b(k) function to check.

% Solving the implicit function of xi* in the Appendix for different values
% for k. 
% Simplifying some terms
A = sqrt(1 / (1/(alpha+beta) + 1/beta));

% Defining the implicit function
f = @(xi_star, k) normcdf( A * (k - ( (beta/(beta+alpha))*xi_star + (alpha/(beta+alpha))*z )) ) ...
                  - (beta/(beta+alpha))*xi_star ...
                  - (alpha/(beta+alpha))*z;

% Range of k values
k_vals = linspace(-10, 10, 200);

% Solve for xi_star for each k
xi_star_vals = NaN(size(k_vals));
for i = 1:length(k_vals)
    k = k_vals(i);

    % Solve f(xi_star, k) = 0 numerically
    try
        xi_star_vals(i) = fzero(@(xi) f(xi, k), 0);  % initial guess = 0
    catch
        xi_star_vals(i) = NaN;  % no solution found
    end
end

% Plot
figure;
plot(k_vals, xi_star_vals, 'LineWidth', 1.5);
xlabel('k');
ylabel('x_i^*');
title('The b(k) function');
grid on;
hold on;

% Had to increase visibility of the plot (thanks chatGPT)

% Compute padded axes limits
padding_x = 0;
padding_y = 0.05 * (max(xi_star_vals) - min(xi_star_vals));
xLimits = [min(k_vals)-padding_x, max(k_vals)+padding_x];
yLimits = [min(xi_star_vals)-padding_y, max(xi_star_vals)+padding_y];

% Fix axes limits
xlim(xLimits);
ylim(yLimits);

% Drawng some lines
line(xLimits, [0 0], 'Color', 'k', 'LineWidth', 0.8);  % x-axis
line([0 0], yLimits, 'Color', 'k', 'LineWidth', 0.8);  % y-axis
line(xLimits, xLimits, 'Color', 'r', 'LineWidth', 0.8);
hold off;