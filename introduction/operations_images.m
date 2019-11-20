close all
clear all
clc

figs_folder = 'figures/';

% Start with an image as follows
n = 10;
f = zeros(n,n);
f(3:7,4:8) = 90;
f(4,5) = 0;
f(8,2) = 90;

% Smooth it out
ftilde = zeros(n,n);
for row = 2:9 % Move over rows
	for col = 2:9 % Move over columns
		xij = f(row,col);
		K = f(row-1:row+1,col-1:col+1);
		K = sum(K(:))./9;
		ftilde(row,col) = K;
	end
end

% Plot
figure
subplot(1,2,1)
imagesc(f)
colormap(gray)
axis square;
axis equal;
axis image;

subplot(1,2,2)
imagesc(ftilde)
colormap(gray)
axis square;
axis equal;
axis image;
export_fig(sprintf('%saveraging_image',figs_folder),'-pdf','-transparent')


close all
clearvars -except figs_folder

% Box-filter
n = 512; % Size of image
F = imresize(rgb2gray(imread('peppers.png')),[n n]);

k = 21; % size of the kernel (the larger, the more blurred)
g = ones(k,k);
g = g/sum(g(:)); % normalize
fg = conv2(F,g,'same');

% Apply filter
figure
subplot(1,2,1)
imagesc(F)
colormap(gray)
title('Original')
axis square;
axis equal;
axis image;
subplot(1,2,2)
imagesc(fg)
colormap(gray)
title('Box-Filtered')
axis square;
axis equal;
axis image;

% Apply median filter
close all
nf = 11;		% Size of filter
J = medfilt2(F,[nf nf]);

figure
subplot(1,2,1)
imagesc(F)
colormap(gray)
title('Original')
axis square;
axis equal;
axis image;
subplot(1,2,2)
imagesc(J)
colormap(gray)
title('Median-Filtered')
axis square;
axis equal;
axis image;
export_fig(sprintf('%smedian_filtered',figs_folder),'-pdf','-transparent');

close all

% Gaussian Filter
sigma = 2;
J = imgaussfilt(F,sigma);

figure
subplot(1,2,1)
imagesc(F)
colormap(gray)
title('Original')
axis square;
axis equal;
axis image;
subplot(1,2,2)
imagesc(J)
colormap(gray)
title('Gaussian-Filtered')
axis square;
axis equal;
axis image;

close all

% Kernel filtering
n = 7;
% Define Impulse Image
F = zeros(n,n);
F(4,4) = 1;

K = linspace(0,1,9);
K = reshape(K,[3 3])';

G = xcorr2(F,K);
Gcon = conv2(F,K);

figure
subplot(2,2,1)
imagesc(F)
colormap(gray)
title('Impulse Image')
axis square;
axis equal;
axis image;
subplot(2,2,2)
imagesc(K)
colormap(gray)
title('Kernel')
axis square;
axis equal;
axis image;
subplot(2,2,3)
imagesc(G)
colormap(gray)
title('Cross-correlation')
axis square;
axis equal;
axis image;
subplot(2,2,4)
imagesc(Gcon)
colormap(gray)
title('Convolution')
axis square;
axis equal;
axis image;
export_fig(sprintf('%scross_corr',figs_folder),'-pdf','-transparent')

% Properties of convolution
%%
close all
clear all
clc

I = zeros(11,11);
I(6,6) = 1;

F = ones(11,11);
for row = 1:size(F,1)
	for col = 1:size(F,2)
		F(row,col) = row*col;
	end
end

figure
subplot(2,2,1)
imagesc(I)
colormap(gray)
title('Identity')
axis square;
axis equal;
axis image;
subplot(2,2,2)
imagesc(F)
title('Matrix')
colormap(gray)
axis square;
axis equal;
axis image;
subplot(2,2,3)
imagesc(conv2(F,I,'same'))
title('Convolution')
colormap(gray)
axis square;
axis equal;
axis image;
subplot(2,2,4)
imagesc(xcorr2(F,I))
title('Cross-correlation')
colormap(gray)
axis square;
axis equal;
axis image;