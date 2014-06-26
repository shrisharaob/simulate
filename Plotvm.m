clear vm cm v i;
dt = 0.05;
fb = '~/Documents/cnrs/simResults/';
vm = importdata([fb, 'vm.csv']);
v = vm(:, 2:end);
%i = vm(:, [2:2:end]);

cm = importdata([fb, 'conMat.csv']);
st = importdata([fb, 'spkTimes.csv']);
t = dt: dt: dt * size(v, 1);
for r = 1 : size(v, 2) 
    for c = 1 : size(v, 2)
        if((cm(r, c) == 1) && (r ~= c))
            plot(t, v(:, r), 'r', t, v(:, c), 'k');
            drawnow;
            drawnow;
            title([num2str(r), '-->', num2str(c)]);
            waitforbuttonpress;
            clf;
        end
    end
end

