close all
clear all
clc

n = 300;
name = 'road2';
f = rescale(load_image(name,n));

figure
imageplot(f)

x0 = [14;161];
x1 = [293;148];
epsilon = 1e-2;
W = epsilon + abs(f-f(x0(1),x0(2)));

figure
imageplot(W)

options.nb_iter_max = Inf;
options.end_points = x1;

[D,S] = perform_fast_marching(1./W, x0, options);

figure
imageplot( convert_distance_color(D,f) );
hold on;
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);