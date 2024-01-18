%% rand_gauss
%
% ret = rand_gauss(Mu, C, n)
%
% Verrauschen eines Vektors mit entsprechender Kovarianzmatrix unter
% Beachtung der Korrelationen
% Erweiterung: Beachte Singularität von Cll
% 
% Mu  = N x 1 mean vector
% C   = N x N Covariance matrix
% n   = number of vectors to be sampled
%
% ret = N x n Matrix of n random vectors 
%
% @author Richard Steffen


function ret = sugr_rand_gauss(Mu, Cll, n)
    
% Zerlegung in Eigenvektoren und Eigenwerte
[R,S] = eig(Cll);   
% Eigenvectoren skalieren mit Wurzeln der Eigenwerte
A = R .* repmat(sqrt(abs(diag(S)))',length(S),1);
ret = repmat(Mu,1,n) + A*randn(length(Mu),n);
  
