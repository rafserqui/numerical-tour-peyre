close all
clear all
clc

%==========================================================================
% Gradient descent method in the simplest fashion
%==========================================================================
tau = 2e-1;
f = @(x)(x.^2);
fgrad = @(x)(2.*x);
x = linspace(-10,10,100);
x0 = -500;
tol = 1e-6;
err = 10000;
it = 1;

x1(1) = x0;

while err > tol
	fprintf('Iteration number %2d \n',it)

	fgradx0 = fgrad(x0);

	fprintf('Gradient at x0 = %3.3f \n',fgradx0)
    
	x1(it+1) = x0-tau.*fgradx0;

	err = abs(fgradx0);

	if err > tol
		x0 = x1(it+1);
	end

	fprintf('New x = %3.3f \n',x0)
	it = it+1;
end

figure
plot(1:it,x1)
ylim([-0.5 0.5])

%==========================================================================
% Gradient descent method in 2-D
%==========================================================================
close all
clear all
clc

set(0,'defaulttextinterpreter','factory')
set(0,'defaultAxesFontName','Palatino LinoType')
figs_folder = 'figures/';

eta = 10;
f = @(x)(0.5.*(x(1).^2 + eta.*x(2).^2));

% Background image of the function
t     = linspace(-0.7,0.7,101);
[u,v] = meshgrid(t,t);
F     = (u.^2 + eta*v.^2)/2;

figure
imagesc(t,t,F)
hold on
contour(t,t,F, 20, 'k')

Gradf = @(x)[x(1); eta*x(2)];

% Step-size tau < 2/eta
tau = 1.8/eta;

% Perform the gradient descent using a fixed step size tau. 
% Display the decay of the energy f(x(k)) through the iteration. 
% Save the iterates so that |X(:,k)| corresponds to x(k).

tol = 1e-5;
x0 = [1;1];
x1 = [];
x1(:,1) = x0;
it = 1;
err = [1000;1000];

while err(1) > tol | err(2) > tol
	fx1(:,it) = f(x0);

	grad = Gradf(x0);

	err = abs(grad);
	x0 = x0 - tau*grad;
	it = it+1;
	x1(:,it) = x0;
end

figure
plot(log10(fx1))
xlim([1 length(fx1)])
title('Log_{10}(x^{(k)})')

figure
imagesc(t,t,F)
hold on
contour(t,t,F, 20,'k-');
plot(x1(1,:), x1(2,:),'k.-');

close all

% For several step-size parameters
symb = {'-','--','-.'};
cmap = lines(4);

taulist = [0.3 1 1.7]/eta;
xinit = [[0.7;0.7],[-0.7;0.5],[-0.7;-0.6]];

figure
imagesc(t,t,F)
colormap(summer(256))
hold on
contour(t,t,F, 20,'k-');
for tt = 1:length(taulist)
	x1 = [];
	x0 = xinit(:,tt);
	tau = taulist(tt);
	x1(:,1) = x0;

	err = [1000;1000];
	it = 1;

	while err > tol
		fx1(:,it) = f(x0);

		grad = Gradf(x0);

		err = max(abs(grad));
		x0 = x0 - tau*grad;
		it = it+1;
		x1(:,it) = x0;
	end
	plot(x1(1,:), x1(2,:),symb{tt},'Color',cmap(tt+1,:),'LineWidth',2);
end
set(gca,'YDir','normal')
axis equal;
export_fig(sprintf('%sgradient_descent_2D',figs_folder),'-pdf','-transparent');