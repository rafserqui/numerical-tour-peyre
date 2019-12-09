close all
clear all
clc
set(0,'defaulttextinterpreter','factory')
set(0,'defaultAxesFontName','Palatino LinoType')
figs_folder = 'figures/';
% Newton's method in unconstrained problems
g = @(x1,x2)(1-x1.^2) + 100.*(x2 - x1.^2).^2;

% Construct grid
npoints = 150;
x1 = linspace(-2,2,npoints);
x2 = linspace(-.5,3,npoints);
[X1,X2] = meshgrid(x1,x2);

G = g(X1,X2);
norcols = perform_hist_eq(G,'linear'); % Normalize colors

figure
s = surf(x2,x1,G,norcols);
s.EdgeColor = 'none';
set(gca, 'Color', 'none'); % Sets axes background
export_fig(sprintf('%srosenbrock',figs_folder),'-pdf','-transparent','-q101')