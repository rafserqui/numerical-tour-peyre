close all
clear all
clc

%==================================================================================================
% Gradient and divergence of Images
%==================================================================================================
n = 256;
x0 = rescale(load_image('lena',n));

% Compute gradient using finite differences
grad = @(x)cat(3,x - x([end 1:end-1],:), x-x(:,[end 1:end-1]));

v = grad(x0);

figure
subplot(2,2,1)
imagesc(v(:,:,1))
axis equal;
axis off;
title('d/dx')

subplot(2,2,2)
imagesc(v(:,:,2))
axis equal;
axis off;
title('d/dy')

subplot(2,2,3)
imageplot(v(:,:,1))
title('d/dx')

subplot(2,2,4)
imageplot(v(:,:,2))
title('d/dy')

