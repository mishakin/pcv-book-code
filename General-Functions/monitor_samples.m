function monitor_samples(i,N_samples)
%% monitor_samples(i,N)
% i   =  current sample
% N   =  total number of samples
%
% wf 3/2018
d=ceil(N_samples/100);

        if mod(i,d)==0
            fprintf(num2str(i)),fprintf(',')
        end
        if mod(i,10*d)==0
            disp(' ')
        end

end

