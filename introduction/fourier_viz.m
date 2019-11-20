close all
clear all
clc

% Suppose time is in milisecs
npoints = 1500;
Fs      = 50;			% Sampling period

% 5 secs in intervals of Fs milisecs
t       = (0:npoints-1).*(1/Fs);	% Time vector

% Parameters of the different signals
signals = 9;
freq0   = 1;
freq1   = 13;
freq    = freq0 + (freq1-freq0).*rand(signals,1);
freq 	= round(freq);

% Circle of radius 1
r=1;
x0=0;
y0=0;
circ = [x0 + r*cos(2.*pi.*t);
        y0 + r*sin(2.*pi.*t)];
    

g = ones(signals,npoints);

figure
hold on 
for ff = 1:length(freq)
	w = freq(ff);
	g(ff,:) = cos(2.*pi.*w.*t);

	% Signal in circle
	om = 1./w;
	gc = g(ff,:).*exp(-2.*pi.*1i.*om.*t);

	subplot(3,3,ff)
	plot(real(gc),imag(gc))
	hold on
	plot(circ(1,:),circ(2,:),'--')
	axis equal;
	xline(0,'--k');
	yline(0,'--k');
	title(['Frequency = ',num2str(w)])
end

% Fourier transform of the signal
g = g';
G = sum(g,2);

% Plot the sum of signals
figure
plot(Fs.*t,G)

% Compute Fourier transform and show the frequencies
Ghat = fft(G);
freq_domain = Fs*(0:length(Ghat)-1)/length(Ghat);

figure
plot(freq_domain,abs(Ghat))
hold on
for ff = 1:length(freq)
	xline(freq(ff),'--k');
end

% Center Fourier to remove signal's negative frequencies
fshift = (-npoints/2:npoints/2-1)*(Fs/npoints);
Ghatshift = fftshift(Ghat);

figure
plot(fshift,abs(Ghatshift))