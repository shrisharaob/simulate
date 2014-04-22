dt = 1e-1;
tStart = 0;
tStop = 100;
t = tStart: dt : tStop;
x = zeros(length(t), 1);
y = zeros(length(t), 1);
z = zeros(length(t), 1);
x(1) = 0;
y(1) = 0;
z(1) = 0;
for n = 2 : length(t)
    x(n) = x(n-1) - dt *  x(n-1);
    %    y(n) = y(n-1) + dt * (-1 * y(n-1) + randn);
    z(n) = z(n-1) + dt * randn;
    y(n) = y(n-1) + sqrt(dt) * randn;
end
plot(t, z, 'b', t, y, 'k');
%plot(t, x,'k', t, y, 'r', t, z, 'g');

