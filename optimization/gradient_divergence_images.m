close all
clear all
clc
set(0,'defaulttextinterpreter','factory')
set(0,'defaultAxesFontName','Palatino LinoType')
figs_folder = 'figures/';

%==================================================================================================
% Gradient and divergence of Images
%==================================================================================================

% Simple example
n = 255;
nm = median(1:n);
I = zeros(n,n);
s = 100;
k = 30;
I(-s+nm:s+nm,-s+nm:s+nm) = 1;
I(-k+nm:k+nm,-k+nm:k+nm) = 0.5;

% Compute gradient using finite differences
grad = @(x)cat(3,x - x([end 1:end-1],:), x-x(:,[end 1:end-1]));

vI = grad(I);

% Example with real image
n = 256;
x0 = rescale(load_image('lena',n));

% Compute gradient using finite differences
grad = @(x)cat(3,x - x([end 1:end-1],:), x-x(:,[end 1:end-1]));

vX = grad(x0);

figure('units','normalized','outerposition',[0 0 0.95 0.95])
subplot(2,3,1)
imagesc(I)
colormap(gray)
axis square;
title('Original')
subplot(2,3,2)
imagesc(vI(:,:,2))
colormap(gray)
axis square;
title('dx')
subplot(2,3,3)
imagesc(vI(:,:,1))
colormap(gray)
axis square;
title('dy')
subplot(2,3,4)
imagesc(x0)
colormap(gray)
axis square;
title('Original')
subplot(2,3,5)
imagesc(vX(:,:,2))
colormap(gray)
axis square;
title('dx')
subplot(2,3,6)
imagesc(vX(:,:,1))
colormap(gray)
axis square;
title('dy')
export_fig(sprintf('%sgradient_descent_image',figs_folder),'-pdf','-transparent');

% Norm
magnitude = sqrt(sum(vX.^2,3));

figure
subplot(1,2,1)
imagesc(x0)
colormap(gray)
title('Original')
axis square
subplot(1,2,2)
imagesc(magnitude)
colormap(gray)
title('Magnitude')
axis square
export_fig(sprintf('%smagnitude_image',figs_folder),'-pdf','-transparent')

% Divergence Operator
div = @(v)v([2:end 1],:,1)-v(:,:,1) + v(:,[2:end 1],2)-v(:,:,2);

% Laplacian Operator
delta = @(x)div(grad(x));

figure
subplot(1,2,1)
imageplot(grad(x0))
colormap(gray)
axis square;
title('Gradient')
subplot(1,2,2)
imagesc(delta(x0))
colormap(gray)
axis square;
title('Laplacian')
export_fig(sprintf('%sgrad_lap',figs_folder),'-pdf','-transparent')