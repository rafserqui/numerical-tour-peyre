close all
clear all
clc

set(0,'defaulttextinterpreter','factory')
set(0,'defaultAxesFontName','Palatino LinoType')
figs_folder = 'figures/';

%==================================================================================================
% Reading and plotting images
%==================================================================================================

% Load and visualize signals and images
n = 256; % Size of image
M = load_image('lena',n);

figure
imageplot(M)


% The command imageplot implements subplots as well. 
figure
imageplot(M(1:50,1:50),'Zoom',1,2,1)
imageplot(-M,'Reversed contrast',1,2,2)

% Alternatively
figure
imageplot({M(1:50,50:100), -M},{'Zoom','Contrast'})

close all;

%==================================================================================================
% Image approximation with Fourier and Wavelets
%==================================================================================================
% Load image
% An image is a matrix f\in R^N of N = N0 x N0 pixels

name = 'lena';
n0   = 512;
f    = rescale(load_image(name,n0));

figure
imageplot(f,'Original Image')


% An image is a 2D array that can be modified just like matrices
figure
imageplot(-f,'Inverse of f',1,2,1)
imageplot(f(n0:-1:1,:),'Turned upside-down',1,2,2)

close all;

% NOTE (from Wikipedia): convolution is a mathematical operation on two functions (f and g) that 
% produces a third function expressing how the shape of one is modified by the other. The term 
% convolution refers to both the result function and to the process of computing it. It is defined 
% as the integral of the product of the two functions after one is reversed and shifted.
% The convolution of f and g is written f∗g, denoting the operator with the symbol ∗. It is 
% defined as the integral of the product of the two functions after one is reversed and shifted. 
% As such, it is a particular kind of integral transform:
% (f*g) \equiv \int^{\inf}_{-\inf} f(tau)g(t-tau)dtau


% Blurring is achieved by computing a convolution (f*g) with a kernel g.
k = 9; % size of the kernel (the larger, the more blurred)
g = ones(k,k);
g = g/sum(g(:)); % normalize

fg = perform_convolution(f,g);
fg2 = conv2(f,g,'same');

% Blurred image
figure
imageplot(fg,'Blurred (Peyre)',1,3,1)
imageplot(fg2,'Blurred (conv2)',1,3,2)
imageplot(f,'Original',1,3,3)
export_fig(sprintf('%sblurr_example',figs_folder),'-pdf','-transparent');

close all

% Sobel edge detector
n = 3;
g = zeros(n,n);
g(:,1) = [1;2;1];
g(:,3) = -g(:,1);

fg = conv2(f,g);

figure
imageplot(fg,'Sobel Edge Detector',1,2,1)
imageplot(f,'Original',1,2,2)


%==================================================================================================
% Fourier transform
%==================================================================================================
% Normalized Fast Fourier Transform
F = fft2(f)/n0;

% Check conservation of the image
fprintf('Energy of image: %5.5f \n',norm(f(:)))
fprintf('Energy of Fourier: %5.5f \n',norm(F(:)))

% Compute the log¡ of the Fourier magnitude for some small epsilon
epsi = 1e-2;

% Shift the zero frequency component to the center of the array
L = fftshift(log(abs(F)+1e-1));

% Display
clf;
imageplot(L,'Log(Fourier Transform)')
export_fig(sprintf('%slog_fourier_transform',figs_folder),'-pdf','-transparent');
%==================================================================================================
% Linear Fourier Approximation
%==================================================================================================
% Perform linear Fourier approximation for several number of coefficients

% Number of kept coefs
nsq = 8;
squares = 1;

% Plot
figure
imageplot(f,'Original',3,3,1)
for nit = 1:nsq
	% Number of coefficients
	M = (squares*2)^2;

	% Bound of interval
	q = sqrt(M);
	
	% Compute Centered Fourier transform
	F = fftshift(fft2(f));
	
	% Linear approximation pre-allocation
	F1 = zeros(n0,n0);
	
	% Choose a square in the middle of the image
	sel = (n0/2-q/2:n0/2+q/2) + 1;
	
	% Take the points of the square to the linear approx
	F1(sel,sel) = F(sel,sel);
	
	% Invert the Fourier and keep real terms
	fM = real(ifft2(fftshift(F1)));
	
	imageplot(fM,{['SNR = ',num2str(snr(f,fM),4)],['N. Coefs = ',num2str(M)]},3,3,nit+1)

	squares = 5+squares;
end
export_fig(sprintf('%sfourier_linear_approx',figs_folder),'-pdf','-transparent');

%==================================================================================================
% Non-linear Fourier
%==================================================================================================
T = linspace(0.25,0.75,8);
F = fft2(f)/n0;
figure
imageplot(f,'Original',3,3,1)
for thr = 1:length(T)
	thresh = T(end-thr+1);
	FT = F.*(abs(F)>thresh);
	fM = real(ifft2(FT)*n0);
	imageplot(clamp(fM), {['SNR=' num2str(snr(f,fM), 4)],['Threshold = ',num2str(thresh,3)]},3,3,thr+1);
end
export_fig(sprintf('%sfourier_nonlinear_approx',figs_folder),'-pdf','-transparent')