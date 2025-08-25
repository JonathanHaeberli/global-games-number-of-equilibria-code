%% Checking for multiplicity
clc; clear; close all;

% Put in any parameter values
alpha = 1;     % precision of public prior
beta = 0.1;    %  precision of private signal
z = 0.5;       % median of common prior
c = 0.8;       % ∈(0,1) cost of attacking
%c* = normcdf((sqrt(beta) / sqrt(beta + alpha)) * ((alpha / sqrt(beta)) * (0.5 - z) - norminv(0.5)));
%c=c* ensures multiple equilibria for any beta < alpha^2 / 2*pi




% Taking the "inverse" of the positive part of the standard normal pdf
x = linspace(0, 10, 1000); % 0 to 5 captures most of the density
pdf_vals = normpdf(x); % Compute standard normal PDF

% Creating the inverse function via interpolation
inv_phi_pos = @(y) interp1(pdf_vals, x, y, 'pchip', NaN);
y_vals = linspace(max(pdf_vals), min(pdf_vals), 1000); 
x_inverse = arrayfun(inv_phi_pos, y_vals);

% Getting theta1 and theta2 according to the formulas described in the
% paper Proposition 1
t2 = normcdf(inv_phi_pos(sqrt(beta)/alpha));
t1 = 1-t2;

% Using the formula in Proposition 2 do determine the c interval for
% multiplicity
term1 = sqrt(beta) / sqrt(beta + alpha);
term2_1 = (alpha / sqrt(beta)) * (t1 - z);
term2_2 = (alpha / sqrt(beta)) * (t2 - z);
term3_1 = norminv(t1); 
term3_2 = norminv(t2);

fktvalue1 = normcdf(term1 * (term2_1 - term3_1)) - c; 
fktvalue2 = normcdf(term1 * (term2_2 - term3_2)) - c;

c_star = fktvalue2 +c;
c_star_star = fktvalue1 +c;

% Determining the number of equilibria applying Proposition 1.
if beta >= ((alpha^2) /(2*pi))
    disp('A unique equilibrium')
elseif beta < ((alpha^2) /(2*pi)) && sign(fktvalue1) == sign(fktvalue2)
    disp('A unique equilibrium')
else
    disp('Multiple equilibria')
end


% Plottinh the b(k) function to check.
% Applying method in Section 5.2.1
x_vals = linspace(-10, 10, 1000);  % Range of x values for plotting
t_vals = zeros(size(x_vals));  % Preallocating space 

v = beta^(1/2);

% First solving for theta*(k) numerically (equation (16))
for i = 1:length(x_vals)
    % Using fminbnd to numerically solve the equation
    % minimizing the difference between normcdf(v*(x - t)) and t
    x = x_vals(i);
    t_vals(i) = fminbnd(@(t) abs(normcdf(v*(x - t)) - t), -10, 10);
end

% Now with values for theta*(k), calculating equation (9) over different values of
% k
y_vals = zeros(size(x_vals)); 
for i = 1:length(x_vals)
    x = x_vals(i);
    t = t_vals(i); 
    %inv_phi_c = norminv(1 - c);  % Inverse CDF of the normal distribution for (1 - c)
    %y_vals(i) = (inv_phi_c / sqrt(beta + alpha) + t - (alpha / (beta + alpha)) * z) / (beta / (beta + alpha));
    y_vals(i) = (t - norminv(c) / sqrt(beta + alpha) - (alpha / (beta + alpha)) * z) * ((beta + alpha) / beta);
    % above I encoded notation of Angeletos, below its mine for Angeletos
    % (2). results are obviously the same
end

% Plotting
figure(1);
plot(x_vals, y_vals, 'b-', 'LineWidth', 2);
xlabel('xj*');
ylabel('xi*');
title('Plot of b(k)');
grid on; 
hold on; % Keep existing plot
plot(x_vals, zeros(size(x_vals)), 'k-', 'LineWidth', 1); % x-axis (y = 0)
plot(zeros(size(y_vals)), x_vals, 'k-', 'LineWidth', 1); % y-axis (x = 0)
plot(x_vals, x_vals, 'r', 'LineWidth', 1);
hold off;

