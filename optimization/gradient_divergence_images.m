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

N = 16;
F = 1:N;
F = reshape(F,[sqrt(N) sqrt(N)])';
ddx = F([end 1:end-1],:);
ddy = F(:,[end 1:end-1]);

% figure('units','normalized','outerposition',[0 0 0.95 0.95])
% subplot(1,2,1)
% imagesc(v(:,:,1))
% axis square;
% axis equal;
% axis image;
% colormap(gray)
% title('d/dx')
% 
% subplot(1,2,2)
% imagesc(v(:,:,2))
% axis square;
% axis equal;
% axis image;
% colormap(gray)
% title('d/dy')
% export_fig(sprintf('%sgradient_descent_image',figs_folder),'-pdf','-transparent');