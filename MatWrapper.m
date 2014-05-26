function MatWrapper(varargin)
% Wrapper for the nw model 
    
    % Save library paths
MatlabPath = getenv('LD_LIBRARY_PATH');
Matpath = getenv('PATH');
% Make Matlab use system libraries
setenv('LD_LIBRARY_PATH','/usr/local/cuda-6.0/lib64:');
setenv('PATH', ['/usr/local/cuda-6.0/bin:/usr/lib/lightdm/lightdm:/' ...
                'usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games']);


N = 5000;
K = 3100:200:5000;
[fldr, contrast, muE, muI, inputTStart, inputTStop, inputStepSize] ...
    = DefaultArgs(varargin, {'~/Documents/cnrs/code/network_model/', [0.25 : 0.25 : 1], 0.1, 0.1, 0, 0, 0});
fldr = '/home/shrisha/Documents/cnrs/results/network_model_outFiles/';
%system(['make -f makecuda']);
for ik = 1 : length(K)
            % system(['./mysolver ', ...
            %         num2str(0), ' ', ...
            %         num2str(0.25), ' ', ...
            %         num2str(muE), ' ', ...
            %         num2str(muI), ' ', ...
            %         num2str(inputTStart), ' ', ...
            %         num2str(inputTStop), ' ', ...
            %         num2str(inputStepSize), ' ', ...
            %         num2str(K(ik))]);
            
            system(['./mysolver ', num2str(K(ik))]);
            system(['mv ', fldr, 'isynapEI ', fldr, 'isynapEI_K', num2str(K(ik)), '_N', num2str(N)]);
            system(['mv ', fldr, 'spkTimes ', fldr, 'spkTimes_K', num2str(K(ik)), '_N', num2str(N)]);
            
end
            %end
            

% Reassign old library paths
%setenv('LD_LIBRARY_PATH',MatlabPath);
%setenv('PATH', Matpath);            
end