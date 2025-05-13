function pop = fn_create_initial_pop(v,eflag,sx,sy,nbins_env,nbins,nmax,b)
% utility function to create an initial population array for use by
% IBM_dispersal function, given an initial displacement kernel to be shared
% by all members of the starting population
%
% INPUTS:
%   v = initial displacement kernel
%   eflag = which environment to use: 1=unbounded, 2=bounded, 3=reef,
%   4=heterogeneous
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   nbins_env = max number of dispersal bins the env can support (must be >= 2)
%   nbins = how many dispersal bins to actually use
%   nmax = maximum larval navigation distance (behavior)
%   b = offspring produced per individual
% 

K = 1;         % carrying capacity per patch

if nbins < 2; error('nbins_use must be at least 2'); end
if nbins > nbins_env; error('nbins_env must be bigger than nbins'); end

% check that the specified kernel is the right size
if nbins ~= size(v,2); error('displacement kernel must have dimensions 1 x nbins'); end

%-----LOAD-ENVIRONMENT----------------------------------------------------%
    if eflag==1
        load(strcat(['../output_environments/env_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '.mat']))
    elseif eflag==2
        load(strcat(['../output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '.mat']))
    elseif eflag==3
        load(strcat(['../output_environments/env_reef_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) 'dmap_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    elseif eflag==4
        load(strcat(['../output_environments/env_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    end
%-----LOAD-ENVIRONMENT----------------------------------------------------%


%-----INITALIZATION-------------------------------------------------------%
    N0 = round(K*S); % initial number of individuals - at capacity

    % create matrix to hold information for all individuals
    pop = zeros(N0,nbins+1);
    pop(:,1:nbins) = ones(N0,nbins).*v; % assign starting dispersal bin values (set equal prob of each distance)
    pop(:,nbins+1) = repmat(via_ID,1,K); % assign evenly to starting patches (patch named by absolute index (in xcoord,etc)


