function [] = fn_create_env_bounded(sx,sy,b)
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018
%
% USAGE: [] = create_env_bounded(sx,sy,nbins_env,nmax,b)
%
% Create the environment to use as a backdrop for the dispersal IBM.
% Here, the environment is a set of viable patches, surrounded by
%   nonviable patches.
% Only need to run this once - it creates a mat file that can be uploaded
%   later within IBM code.
%
% INPUTS:
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   b = number of offspring per site in the environment
%
% SAVED VARIABLES:
%   E = 2D matrix of patches, marked as viable (1) and nonviable (0)
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   S = number of actual viable sites in the environment
%   xcoord = list of x-coordinates for all patches in E
%   ycoord = list of y-coordinates for all patches in E
%   via_ID = list of viable patches in E, marked by their index number
%   bvec = vector of offspring capacity of each viable site, 1...S
%   dists = array of distances between each site

    stmp = 12+1; % make sure environment can hold all max dispersal distances (up to 12)

% create the environment
    E_via = ones(sy,sx);                 % viable patches in the environment
    E = zeros(sy+2*stmp,sx+2*stmp);      % full environment
    E(stmp+1:stmp+sy,stmp+1:stmp+sx) = E_via;
    via_ID = find(E);                  % ID of viable patches within E
    S = sum(sum(E));                   % total number of viable patches

% get dimensions and coordinates off the environment
    xdim = size(E,2);     % number of patches in x dimension
    ydim = size(E,1);     % number of patches in y dimension
    xcoord = repmat(1:xdim,ydim,1);
    xcoord = xcoord(:);                % vector of x-coordinates for each patch
    ycoord = repmat([1:ydim]',1,xdim);
    ycoord = ycoord(:);                % vector of y-coordinates for each patch

% %% create array of distances
    dists=zeros(length(xcoord),length(xcoord));
    for ii = 1:length(xcoord)
        for i = 1:length(xcoord)
            dists(ii,i) = abs(xcoord(ii)-xcoord(i)) + abs(ycoord(ii)-ycoord(i)); % vonNeumann
            %dists(ii,i)=((xcoord(ii)-xcoord(i))^2+(ycoord(ii)-ycoord(i))^2)^.5; % Euclidean
        end
    end
    %dists=floor(dists);

% create vector of the number of offspring produced at each site
bvec = b*ones(S,1);

% %% create vector patches, with birth rate of each site (viable), or 0
% (nonviable)
patches = zeros(length(xcoord),1);
patches(via_ID)=bvec; %added to create_env_bounded

clear E_via i ii ind stmp tmp xcoord_via ycoord_via xdim ydim y

save(strcat(['../output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_b=' num2str(b) '.mat']),'-v7.3')