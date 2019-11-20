close all
clear all
clc

figs_folder = 'figures/';

%==================================================================================================
% Gradient and divergence of Images
%==================================================================================================
n = 256;
x0 = rescale(load_image('lena',n));

% Compute gradient using finite differences
grad = @(x)cat(3,x - x([end 1:end-1],:), x-x(:,[end 1:end-1]));

v = grad(x0);

[FX,FY] = gradient(x0);

% figure('units','normalized','outerposition',[0 0 0.95 0.95])
% subplot(2,2,1)
% imagesc(v(:,:,1))
% axis equal;
% title('d/dx')
% 
% subplot(2,2,3)
% imagesc(v(:,:,2))
% axis equal;
% title('d/dy')
% 
% subplot(2,2,2)
% imagesc(FX)
% axis equal;
% title('d/dx')
% 
% subplot(2,2,4)
% imagesc(FY)
% axis equal;
% title('d/dy')
% export_fig(sprintf('%sgradient_descent_image',figs_folder),'-pdf','-transparent');

dx = abs(FX - v(:,:,1));
dx = abs(FY - v(:,:,2));

[Gx,Gy] = imgradientxy(x0);

figure
subplot(2,2,1)
imagesc(Gx)
title('dx')

subplot(2,2,2)
imagesc(v(:,:,1))
title('dx')

subplot(2,2,3)
imagesc(Gy)
title('dy')

subplot(2,2,4)
imagesc(v(:,:,2))
title('dy')