function [mortality_cost, kincomp_cost, fitness_sitecapture,off] = fn_one_costben_sim(b,p,nbins,eflag,sx,sy,nmax,v)
% written by E Schlatter (eschlatt@bu.edu)
%
% For a specified displacement kernel, seascape, and set of parameters
% parameters, simulates one generation of dispersal and calculates 
% costs and fitness for each displacement distance.
%
% INPUTS:
%   b = offspring produced per individual
%   p = probability of surviving dispersal
%   nbins = how many dispersal bins to actually use
%   eflag = which environment to use: 1=unbounded, 2=bounded, 3=reef,
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   nmax = maximum larval navigation distance (behavior)
%   v = displacement kernel
%
% OUTPUTS:
% - mortality_cost: nonviable habitat cost (proportion of larvae that die because they can't reach
% suitable habitat)
% - kincomp_cost: kin competition cost (expected number of siblings encountered by a
% potentially-recruiting larva at its destination site, divided by its total number of siblings)
% - fitness_sitecapture: the proportion of larvae displacing each distance
% that survive through competition to successfully recruit


% check that kernel is correctly specified
if abs(sum(v)-1)>.01; error('kernel (v) must sum to 1'); end
if length(v)~=nbins; error('kernel (v) length must match nbins'); end

%-----LOAD-ENVIRONMENT----------------------------------------------------%
    if eflag==1
        load(strcat(['../output_environments/env_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '.mat']))
    elseif eflag==2
        load(strcat(['../output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_b=' num2str(b) '.mat']))
    elseif eflag==3
        load(strcat(['../output_environments/env_reef_nbins=' num2str(nbins) '_nmax=' num2str(nmax) 'dmap_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    elseif eflag==4
        load(strcat(['../output_environments/env_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    end
%-----LOAD-ENVIRONMENT----------------------------------------------------%

%-----INITALIZATION-------------------------------------------------------%
    K = 1;      % carrying capacity per patch
    N0 = round(K*S); % initial number of individuals - at capacity

    % create matrix to hold information for all individuals
    pop = zeros(N0,nbins+1);
    pop(:,1:nbins) = repmat(v,N0,1); % assign starting displacement kernel
    pop(:,nbins+1) = repmat(via_ID,1,K); % assign evenly to starting patches (patch named by absolute index (in xcoord,etc)

%-----INITALIZATION-------------------------------------------------------%



%-----SIMULATE------------------------------------------------------------%

%-----REPRODUCTION-----%
    % each individual produces b offspring with same patch and dispersal
    % parameters as parent
    off = [];
    for j = 1:size(pop,1)
        off = [off; repmat(pop(j,:),[patches(pop(j,nbins+1)),1])];  % !!! patches contains the number of larvae per patch (b=10)
    end

    off(1,nbins+2:nbins+5)=[0 0 0 0]; %add columns to off to store more info
    Noff = size(off,1); % total number of offspring

    % off: columns
    %  1:nbins > displacement kernel
    %  nbins+1 > starting patch
    %  nbins+2 > displacement patch (or 0, if died during displacement)
    %  nbins+3 > dispersal patch (or 0, if died during navigation)
    %  nbins+4 > recruitment patch (or 0, if died during competition)
    %  nbins+5 > displacement distance


%-----REPRODUCTION-----%


%-----DISPERSAL-----%
srand = rand(Noff,1); % random numbers to use for survival probabilities
surv = srand<p;        % which offspring survive (1) or not (0), if they attempt to disperse
drand = sum(v).*rand(Noff,1); % random numbers to use for dispersal probabilities -- scaled in case the kernel sums to slightly less than 1
d_ind = zeros(Noff,1); % to hold the index of dispersal distance traveled (1=stay, 2=distance 1, 3=distance 2, etc)

for j = 1:Noff

    % find dispersal bin corresponding to random dispersal distance
    % cumsum is the cumulative sum across the disperal bins -- i.e. the
    %   cumulative probability distribution
    % this finds the first bin where cumsum exceeds the random
    %   dispersal number generated
    % e.g. if off(j,1:nbins) = [0.1 0.8 0.1] this is 10%, 80%, and
    %      10% probability of traveling distances 0, 1, and 2
    %      and cumsum(off(j,1:nbins)) = [0.1 0.9 1]
    %      if drand(j) = 0.2 then d_ind(j) would be 2, that is the
    %      individual travels distance 1 during dispersal
    d_ind(j) = find(cumsum(off(j,1:nbins),2)>drand(j),1,'first');

end
clear ind val j

off(:,nbins+5) = d_ind; % store displacement distances in off
% d_ind starts at 1 for offspring who stay home

ind = find(d_ind==1); % which offspring didn't leave natal patch during displacement
off(ind,nbins+2) = off(ind,nbins+1); % their patch post-displacement is same as patch pre-displacement
clear ind

for j = 2:nbins     % for each of the possible displacement distances
    ind = find(d_ind==j); % find all offspring that drew j-1 displacement distance
    for i = ind'
        if surv(i)==1 % for offspring that survived displacement
            x = find(dists(:,off(i,nbins+1))==(j-1)); %possible patches to displace to
            y = randi(length(x)); % pick one at random
            off(i,nbins+2)=x(y); % save landing patch
        end
    end

end

% navigation
for i = 1:size(off,1) %for each larva
    if patches(off(i,nbins+2))==0 %if the larva has displaced to an uninhabitable patch
        x = find(dists(:,off(i,nbins+2))<=nmax); %find patches within navigation distance
        x_hab = x(patches(x)~=0); %restrict to habitable patches
        if ~isempty(x_hab) %if there are any habitable patches
            y = randi(length(x_hab)); %pick one at random
            off(i,nbins+3) = x_hab(y); %save settlement patch
        end
    else % if the larva has displaced to a habitable patch
        off(i,nbins+3) = off(i,nbins+2); %stay there
        % this also puts a 0 in the column nbins+3 for larvae that died
        % during displacement (i.e., off(i,nbins+2)=0), because
        % patches(0)==0 (the if condition above) is false
    end
end

clear srand surv drand d_ind ind i j
%-----DISPERSAL-----%


%-----COMPETITION-----%
% randomly reorder all offspring (to avoid spurious patterns in next
% step)
Noff = size(off,1); % total number of offspring
xind = randperm(Noff);
off = off(xind,:);
off(:,nbins+4) = off(:,nbins+3); % set all recruitment to dispersal patch

% only allow K offspring per patch to survive
Noffs = sum(off(:,nbins+3)==via_ID'); % number of offspring per patch
fullind = find(Noffs>K);  % index of overcrowded patches
for i = 1:length(fullind) % loop over each overcrowded patch
    patch_ind = find(off(:,nbins+3)==via_ID(fullind(i))); % index of offspring in patch
    off(patch_ind(K+1:end),nbins+4) = 0; % kill all but K of these
end
clear fullind patch_ind i

%-----COMPETITION-----%


%-----CALCULATE OUTPUTS-----%

    % site-capturing fitness, h(d)
    % probability a larva that displaces distance d captures its site

    fitness_sitecapture = zeros(nbins,1);

    for i = 1:nbins
        off_dist = off(off(:,nbins+5) == i,:); % pick just the ones that displaced the focal distance
        n_winners = sum(off_dist(:,nbins+4) ~= 0); % how many of them survived competition
        fitness_sitecapture(i) = n_winners/size(off_dist,1); % proportion
    end
    clear n_winners off_dist

    % direct cost (proportion of larvae that die because they can't reach
    % suitable habitat)

    mortality_cost = zeros(nbins,1);

    for i = 1:nbins
        off_dist = off(off(:,nbins+5) == i,:); % pick just the ones that displaced the focal distance
        died_in_nav = sum(off_dist(:,nbins+3)==0); % how many didn't find habitat
        mortality_cost(i) = died_in_nav/size(off_dist,1); % proportion
    end
    clear inds off_dist died_in_nav

    % indirect cost (expected number of siblings encountered by a
    % potentially-recruiting larva at its destination site, divided by total
    % number of competitors at that site)

    kincomp_cost = zeros(nbins,1);

    for i = 1:nbins
        % choose the rows that represent offspring that displaced the focal
        % distance AND survived to compete
        off_dist = off(off(:,nbins+5) == i,:);
        off_dist = off_dist(off_dist(:,nbins+3) ~= 0,:);

        kincost_i=[];
        % pick each individual in off_dist. Find its destination site. Then
        % count the total number of individuals (from off) in that
        % destination site, and the number that share the focal
        % individual's origin site.
        for j = 1:size(off_dist,1)
            origin = off_dist(j,nbins+1);
            destination = off_dist(j,nbins+3);
            %n_competitors = sum(off(:,nbins+3)==destination);
            sibs = sum((off(:,nbins+1)==origin).*(off(:,nbins+3)==destination))-1;
            %kincost_i = [kincost_i sibs/n_competitors];
            kincost_i = [kincost_i sibs/(b-1)];
        end
        clear origin destination

        kincomp_cost(i) = mean(kincost_i);
        clear off_dist
    end

%-----CALCULATE OUTPUTS-----%

end