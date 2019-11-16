%==================================================================================================
% Dijkstra and Fast Marching Algorithms
%==================================================================================================
close all
clear all
clc

% For the first time only
% addpath('../toolboxes/toolbox_signal')
% addpath('../toolboxes/toolbox_general')
% addpath('../toolboxes/toolbox_graph')

% Navigating on the grid
% We use a cartesian grid of size n*n, and defines operators to navigate in the grid.
% Use a single index i in {1,2,...,n^2} to index a position on the 2D grid.

n = 40;

% Four displacement vector to go to the four neighbors
neigh = [[1;0] [-1;0] [0;1] [0;-1]];

% Use periodic boundary conditions
bounds = @(x)mod(x-1,n)+1;

% For a given index |k| and a given neighboring index k in 1,2,3,4, |Neigh(k,i)| gives the
% corresponding grid neighboring index.

ind2sub1 = @(k)[rem(k-1, n)+1; (k - rem(k-1,n) - 1)/n + 1];
sub2ind1 = @(u)(u(2) - 1)*n + u(1);
Neigh = @(k,i)sub2ind1(bounds(ind2sub1(k) + neigh(:,i)));

%==================================================================================================
% Dijkstra Algorithm
%==================================================================================================

% The Dijkstra algorithm compute the geodesic distance on a graph. We use here a graph whose nodes
% are the pixels, and whose edges defines the usual 4-connectity relationship.
% 
% In the following, we use the notation i∼j to indicate that an index j is a neighbor of i on the
% graph defined by the discrete grid.
% 
% The metric W(x). We use here a constant metric.

W = ones(n);

% Set S = {x0} of initial points
x0 = [n/2;n/2];

% Initialize the stack of available indexes
I = sub2ind1(x0);

% Initialize the distance to +Inf excepted for the boundary conditions
D = zeros(n)+Inf;
D(I) = 0;

% Initialize the state to 0 (unexplored), excepted for the boundary point S (indexed by |I|) to 1 
% (front)
S = zeros(n);
S(I) = 1;

% 1.- Pop the from stack the element i with smallest current distance Di
[tmp,j] = sort(D(I));
j = j(1);
i = I(j);
I(j) = [];

% 2.- Update its state S to be dead (-1)
S(i) = -1;

% 3.- Retrieve list of the four neighbors
J = [Neigh(i,1); Neigh(i,2); Neigh(i,3); Neigh(i,4)];

% 4.- Remove those that are dead (no need to consider again)
J(S(J) == -1) = [];

% 5.- Those that are not yet considered (S(J) = 0) to the front stack I (state 1)
J1 = J(S(J) == 0);
I  = [I; J1];
S(J1) = 1;

% 6.- Update neighbor values. For each neighbor j of i, perform the update assuming the length of 
% the edge between j and k is Wj
% Dj <- min_{k~j} Dk + Wj

for j = J'
	dx = min(D([Neigh(j,1) Neigh(j,2)]));
	dy = min(D([Neigh(j,3) Neigh(j,4)]));
	D(j) = min(dx+W(j), dy+W(j));
end

% Implement the Dijkstra algorithm by iterating these step while the stack |I| is non empty.
% Display from time to time the front that propagates.

options.method = 'dijkstra';
options.svg_rate = n*6;
[D,Dsvg,Ssvg] = perform_dijkstra_fm(W, x0, options);

for p=1:4
	figure(1)
	subplot(2,2,p)
	d = Dsvg(:,:,2+p);
	d(d==Inf) = 0;
	imagesc(d)

	figure(2)
	subplot(2,2,p)
	imageplot(d)
	colormap jet(256)
end


%==================================================================================================
% Fast Marching Method Algorithm
%==================================================================================================
% The implementation of the FMM is basically the same as Dijkstra's Algorithm with a continuous 
% surface. The key is that the FMM updates graph in a different manner. 
% Dijkstra's updating => D(j) = min(dx+W(j), dy+W(j))
% FMM updating => Eikonal Update
% Delta = 2*W(j) - (dx - dy)^2;
% if Delta >= 0
% 	D(j) = (dx + dy + sqrt(Delta))/2;
% else
% 	D(j) = min(dx+W(j), dy+W(j));
% end

% Implement FMM as before

options.method = 'fm';
options.svg_rate = n*6;
[D,Dsvg,Ssvg] = perform_dijkstra_fm(W, x0, options);

for p = 1:4
	figure(3)
	subplot(2,2,p)
	d = Dsvg(:,:,2+p);
	d(d==Inf) = 0;
	imageplot(d)
	colormap jet(256)

	figure(4)
	subplot(2,2,p)
	imagesc(d)
end

%==================================================================================================
% Geodesic paths 
%==================================================================================================
% Suppose now there is a bump in the middle (a mountain).
n = 100;
x = linspace(-1,1,n);
[Y, X] = meshgrid(x,x);
sigma = .2;
W = 1 + 8.*exp(-(X.^2 + Y.^2)/(2*sigma^2));

% Plot it
close all
clf
figure(1)
imageplot(W);

figure(2)
imagesc(W);

% Starting points
x0 = round([.1;.1]*n);

% Compute distance map to these starting point using FMM
options.method = 'fm';
D = perform_dijkstra_fm(W, x0, options);
k = 8;

displ = @(D)cos(2*pi*k*D/max(D(:)));

figure(3)
imageplot(displ(D))
colormap jet(256)

figure(4)
imagesc(displ(D))

% Once the geodesic distance map to S has been computed, the geodesic curve between any point
% x1 and S is extracted through gradient descent,
% Compute the gradient G_0(x) = Nabla D(x) in R^2 of the distance map. Use centered differences.
options.order = 2;
G0 = grad(D,options);

% Normalize gradient to obtained G(x) = G0(x)/norm(G0(x)) in order to have unit speed geodesic 
% curve (parameterized by arc length)
G = G0./repmat(sqrt(sum(G0.^2, 3)),[1 1 2]);

% The geodesic curve is then numerically computed using a discretized gradient descent.
% Step size tau for the gradient descent
tau = .8;

% Initialize the path with an ending point
x1 = round([.9;.88]*n);
gamma = x1;

% Define a shortcut to interpolate G at 2D points. interp2 switches the role of the axes
Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1));
			interp2(1:n,1:n,G(:,:,2),x(2),x(1))];

% Compute the gradient at the last point in the path using interpolation
g = Geval(G,gamma(:,end));

% Perform the descent and add the new point to the path
gamma(:,end+1) = gamma(:,end) - tau*g;

% Perform the full geodesic path extraction by iterating the gradient descent. You must be very
% careful when the path become close to x0, because the distance function is not differentiable at
% this point. You must stop the iteration when the path is close to x0.

gamma = x1;

for p = 1:1.5*n/tau
	gamma(:,end+1) = gamma(:,end) - tau*Geval(G,gamma(:,end));
	if norm(gamma(:,end) - x0) < 1
		break;
	end
end

% Display geodesic curve
figure
subplot(1,2,1)
imageplot(W)
hold on 
colormap gray(256)
plot(gamma(2,:),gamma(1,:),'--','LineWidth',2)
plot(x0(2),x0(1),'.r','MarkerSize',25)
plot(x1(2),x1(1),'.b','MarkerSize',25)
title('Using imageplot Command')

subplot(1,2,2)
imagesc(W)
hold on 
colormap gray(256)
plot(gamma(2,:),gamma(1,:),'--','LineWidth',2)
plot(x0(2),x0(1),'.r','MarkerSize',25)
plot(x1(2),x1(1),'.b','MarkerSize',25)
title('Using imagesc Command')