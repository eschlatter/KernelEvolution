function []=fn_many_costben_sims(kernels,b,p,nbins,eflag,sx,sy,nmax,saveto_filepath)
% written by E Schlatter (eschlatt@bu.edu)
%
% Runs the function one_costben_sim many times, each for a different
% kernel, and outputs nonviable habitat and kin competition costs and 
% benefit for each kernel
%
% INPUTS:
%   kernels = matrix of kernels to use (each row a kernel, each column a distance
%   probability)
%   b = offspring produced per individual
%   p = probability of surviving dispersal
%   nbins = how many dispersal bins to actually use
%   eflag = which environment to use: 1=unbounded, 2=bounded, 3=reef,
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   nmax = maximum larval navigation distance (behavior)
%   saveto_filepath = where to save output files

out_fitness = zeros(size(kernels,1),1);
out_kincost = zeros(size(kernels,1),1);
out_mortcost = zeros(size(kernels,1),1);
dist_fitness = zeros(size(kernels,1),nbins);
dist_kincost = zeros(size(kernels,1),nbins);
dist_mortcost = zeros(size(kernels,1),nbins);

for i=1:size(kernels,1)
    v = kernels(i,:);
    [M,K]=fn_one_costben_sim(b,p,nbins,eflag,sx,sy,nmax,v);
    B = (1-M).*(1-K); % total benefit
    dist_fitness(i,:) = B;
    B(isnan(B))=0; % NaN entries mean no larvae displaced that distance,
    out_fitness(i) = v*B;
    dist_kincost(i,:) = K;
    dist_mortcost(i,:) = M;
    M(isnan(M))=0; % NaN entries mean no larvae displaced that distance,
    K(isnan(K))=0; % so no contribution to costs/benefits
    out_kincost(i) = v*K;
    out_mortcost(i) = v*M;
    clear K M B v
end

save(strcat([saveto_filepath,'.mat']),"out_fitness","out_kincost","out_mortcost","dist_fitness","dist_kincost","dist_mortcost","kernels")