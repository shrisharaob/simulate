dt = 1e-3;
t = 0:dt:1e5;
tau = 3;
spks = 10;
y = zeros(1, length(t));
z = zeros(1, length(t));
tic
for n = 2 : length(t)
	  y(n) = (t(n) > spks) * (exp(-1 * (t(n) - spks) / tau));
end
toc
disp('done')
tic
  cnst =   exp(-1 * dt / tau);
for n = 2 : length(t)
	  z(n) = z(n - 1) * cnst; 
if(t(n) == spks)
  z(n) = z(n) + 1;
end
end
toc
%%figure, plot(t, y, 'k');
%%hold on, plot(t, z, 'or');
