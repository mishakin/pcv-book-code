%% init_rand_seed: initialize random number generation by specific seed
% to ensure repeatability of experiments
%
% seed = init_rand_seed(seed)
% seed = int, seed for initialization of rand and randn
%        default 42
%        if seed = 0 rand is initialized by timestamp
%
% Susane Wenzel 09/16
% wenzel@igg.uni-bonn.de

function seed = init_rand_seed(seed)

if nargin<1
    seed = 42;   
end
if seed == 0
    seed =  round( sum(10000*clock) );
end;

disp(['Random seed = ', num2str(seed)])

if verLessThan('matlab', '15')
    rand('state',seed); %#ok<*RAND>
    randn('state',seed);
else
    rng(seed,'v5uniform');
end




