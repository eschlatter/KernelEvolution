function [] = fn_IBM_dispersal(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del,pop_init,saveto_filepath)
% written by Allison K. Shaw (ashaw@umn.edu)
% updated December 2023 by E Schlatter
%
% IBM version of Hamilton & May model with spatially-explicit patches in
% 2-D
% Two steps to the dispersal process: individuals go from natal patch to
%   landing patch, then move from landing patch to final patch
%  > S sites in the environment, each supporting K individuals
%  > each adult produces b offspring per generation
%  > each individual is defined by set of dispersal strategies (probability
%      of traveling distance x)
%  > each dispersing individual has a probability p of surviving dispersal
%      (regardless of distance traveled)
%  > offspring inherit dispersal strategy from parent with small mutation
%  > individuals whose landing patch is not viable or who disperse beyond
%     the edge of the world die...except if they navigate back to a close
%     viable patch
%
% INPUTS:
%   gflag = whether (1) or not (0) to show plots during simulation
%   eflag = which environment to use: 1=unbounded, 2=bounded, 3=reef,
%   4=heterogeneous
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   nbins_env = max number of dispersal bins the env can support (must be >= 2)
%   nbins = how many dispersal bins to actually use
%   nmax = maximum larval navigation distance (behavior)
%   G = number of total generations to simulate
%   b = offspring produced per individual
%   p = probability of surviving dispersal
%   del = fraction of dispersal probability to move during mutation
%   pop_init = initial population array (e.g., from previous simulation)
%   saveto_filepath = where to save output files

rng('shuffle') % seed the random number generator from computer clock

K = 1;         % carrying capacity per patch

if nbins < 2; error('nbins_use must be at least 2'); end
if nbins > nbins_env; error('nbins_env must be bigger than nbins'); end

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

    % create initial population (incl. locations and kernels)
    % check if pop_init is the right dimensions and contains only viable patches
    if isequal(size(pop_init),[N0,nbins+1]) && sum(ismember(pop_init(:,nbins+1),via_ID))==N0 
        pop=pop_init;
    else
        pop_init=0;
        sprintf('pop_init not correctly specified; using default starting population')
    end

    if pop_init==0 % if no initial population matrix specified
        % create matrix to hold information for all individuals
        pop = zeros(N0,nbins+1);
        pop(:,1:nbins) = (1/nbins)*ones(N0,nbins); % assign starting dispersal bin values (set equal prob of each distance)
        pop(:,nbins+1) = repmat(via_ID,1,K); % assign evenly to starting patches (patch named by absolute index (in xcoord,etc)
    end


    % matrices to record population parameters over time
    dtime_avg = zeros(G,nbins); % population mean of dispersal parameters
    dtime_std = zeros(G,nbins); % pop. standard deviation of dispersal parammeters
    fitness_sitecapture = zeros(G,nbins);
    mortality_cost = zeros(G,nbins); % mortality cost at each distance
    kincomp_cost = zeros(G,nbins); % kin competition cost at each distance

    % if set to display graphics, display first figure and wait for keystrike
    if gflag==1
        figure(1); clf
        plot(pop(:,1:nbins)','o')
        ylim([0 1])
        xlabel('dispersal bin')
        ylabel('probability')
        pause(0.1)
    end
%-----INITALIZATION-------------------------------------------------------%

tic

%-----SIMULATE------------------------------------------------------------%
g=0; % this counts the number of generations that have passed

while g<G && size(pop,1)>0 % loop over generations (only while population not extinct)
    g = g+1 % step generation forward

    % record statistics
    dtime_avg(g,:) = mean(pop(:,1:nbins),1);   % average dispersal params
    dtime_std(g,:) = std(pop(:,1:nbins),[],1); % std. of dispersal params

    %-----REPRODUCTION-----%
        % each individual produces b offspring with same patch and dispersal
        % parameters as parent
        off = [];
        for j = 1:size(pop,1)
            off = [off; repmat(pop(j,:),[patches(pop(j,nbins+1)),1])];
        end

        off(1,nbins+2:nbins+5)=[0 0 0 0]; %add columns to off to store more info
        %  columns of off:
        %  1:nbins -> displacement kernel
        %  nbins+1 -> starting patch
        %  nbins+2 -> displacement patch (or 0, if died during displacement)
        %  nbins+3 -> dispersal patch (or 0, if died during navigation)
        %  nbins+4 -> recruitment patch (or 0, if died during competition)
        %  nbins+5 -> displacement distance (1=stay, 2=distance 1, 3=distance 2, etc)
        Noff = size(off,1); % total number of offspring

    %-----REPRODUCTION-----%


    %-----MUTATION-AND-DISPERSAL-----%
    srand = rand(Noff,1); % random numbers to use for survival probabilities
    surv = srand<p;        % which offspring survive (1) or not (0), if they attempt to disperse
    drand = rand(Noff,1); % random numbers to use for dispersal probabilities

    % assign displacement distances
    for j = 1:Noff

        % mutate dispersal strategy slightly by taking or adding del from dispersal bin
        % if del amount isn't remaining in dispersal bin, take everything
        % (i.e. bound dispersal probability at zero)
        ind = randperm(nbins,2);         % select two bins at random
        val = min(off(j,ind(1)),del);          % take del of prob, if available
        off(j,ind) = off(j,ind) + [-val +val]; % move del from first to second

        % find dispersal bin corresponding to random dispersal distance
        % cumsum is the cumulative sum across the disperal bins -- i.e. the
        %   cumulative probability distribution
        % this finds the first bin where cumsum exceeds the random
        %   dispersal number generated
        % e.g. if off(j,1:nbins) = [0.1 0.8 0.1] this is 10%, 80%, and
        %      10% probability of traveling distances 0, 1, and 2
        %      and cumsum(off(j,1:nbins)) = [0.1 0.9 1]
        %      if drand(j) = 0.2 then off(j,nbins+5) would be 2, that is the
        %      individual travels distance 1 during dispersal
        off(j,nbins+5) = find(cumsum(off(j,1:nbins),2)>drand(j),1,'first');

    end

    clear ind val j

    % assign displacement landing patches
    ind = find(off(:,nbins+5)==1); % which offspring didn't leave natal patch during displacement
    off(ind,nbins+2) = off(ind,nbins+1); % their patch post-displacement is same as patch pre-displacement
    clear ind

    for j = 2:nbins     % for each of the possible displacement distances
        ind = find(off(:,nbins+5)==j); % find all offspring that drew j-1 displacement distance
        for i = ind'
            if surv(i)==1 % for offspring that survived displacement
                x = find(dists(:,off(i,nbins+1))==(j-1)); %possible patches to displace to
                y = randi(length(x)); % pick one at random
                off(i,nbins+2)=x(y); % save landing patch
                clear x y
            end
        end
    end
    clear ind j

    % navigation
    for i = 1:size(off,1)
        if surv(i)==1 %only larvae that didn't die during displacement
            if patches(off(i,nbins+2))==0 %if the larva has displaced to an uninhabitable patch
                x = find(dists(:,off(i,nbins+2))<=nmax); %find patches within navigation distance
                x_hab = x(patches(x)~=0); %restrict to habitable patches
                if ~isempty(x_hab) %if there are any habitable patches
                    y = randi(length(x_hab)); %pick one at random
                    off(i,nbins+3) = x_hab(y); %save settlement patch
                end
            else % if the larva has displaced to a habitable patch
                off(i,nbins+3) = off(i,nbins+2); %stay there
            end
        end
    end

    clear srand surv drand i
    %-----MUTATION-AND-DISPERSAL-----%


    %-----COMPETITION-----%
    % randomly reorder all offspring (to avoid spurious patterns in next
    % step)
    Noff = size(off,1); % total number of offspring
    xind = randperm(Noff);
    off = off(xind,:);

    % only allow K offspring per patch to survive
    for i = 1:length(via_ID) % loop over each patch
        patch_ind = find(off(:,nbins+3)==via_ID(i)); % index of offspring currently in patch
        if(~isempty(patch_ind))
            n_settlers = min(K,size(patch_ind,1));
            off(patch_ind(1:n_settlers),nbins+4) = off(patch_ind(1:n_settlers),nbins+3); %choose K to survive (gets a patch in column nbins+4 = recruitment column)
        end
    end

    % save recruited offspring as new population
    remaining = find(off(:,nbins+4)~=0);
    pop = off(remaining,[1:nbins,nbins+3]);
    
    clear patch_ind n_settlers remaining i
    %-----COMPETITION-----%


    %-----OUTPUTS---------%

    % if set to display graphics, update figure 1
    if gflag==1
        N = size(pop,1);  % count number of individuals present
        figure(1); clf
        plot(pop(:,1:nbins)','o')
        ylim([0 1])
        xlabel('dispersal bin')
        ylabel('probability')
        title(strcat(['generation = ' num2str(g)]))
        pause(0.1)
    end
    
    % site-capturing fitness, h(d)
    % probability a larva that displaces distance d captures its site
    for i = 1:nbins
        off_dist = off(off(:,nbins+5) == i,:); % pick just the ones that displaced the focal distance
        n_winners = sum(off_dist(:,nbins+4) ~= 0); % how many of them survived competition
        fitness_sitecapture(g,i) = n_winners/size(off_dist,1); % proportion
    end
    clear n_winners off_dist

    % direct cost (proportion of larvae that die because they can't reach suitable habitat)
    for i = 1:nbins
        %%%%%%%%%%%%%%%%%%%%% should this be offspring with DISPERSAL distance i? %%%%%%%%%%%%%%%%%%%%
        %%% No. Because we're looking at what DISPLACEMENT kernel evolves
        %%% -- i.e., the fitness of each displacement distance.
        off_dist = off(off(:,nbins+5) == i,:); %offspring with displacement distance i
        died_in_nav = sum(off_dist(:,nbins+3)==0); %how many didn't reach suitable habitat
        mortality_cost(g,i) = died_in_nav/size(off_dist,1);
    end
    clear off_dist died_in_nav

    
    % indirect cost (expected number of siblings encountered by a
    % potentially-recruiting larva at its destination site, divided by its 
    % total number of siblings)
    for i = 1:nbins
        % choose the offspring that displaced the focal distance AND survived to compete
        off_dist = off(off(:,nbins+5) == i,:);
        off_dist = off_dist(off_dist(:,nbins+3) ~= 0,:);
        %vector to store the kin comp cost for each of them
        kincost_i=zeros(size(off_dist,1),1);
        
        % Pick each individual in off_dist. Find its destination site. 
        % Then count the number of individuals at the destination site 
        % that share the focal individual's origin site.
        for j = 1:size(off_dist,1)
            origin = off_dist(j,nbins+1);
            destination = off_dist(j,nbins+3);
            sibs = sum((off(:,nbins+1)==origin).*(off(:,nbins+3)==destination))-1; % subtract one for the focal individual
            kincost_i(j) = sibs/(b-1);
        end
        
        kincomp_cost(g,i) = mean(kincost_i);
        clear origin destination off_dist kincost_i
    end

    clear i j
    %-----OUTPUTS---------%

end % generation loop

toc

clear dists 

    % save output as mat file
    if eflag==1
        save(strcat([saveto_filepath '/IBM_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(p) '.mat']))
    elseif eflag==2
        save(strcat([saveto_filepath '/IBM_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(p) '.mat']))
    elseif eflag==3
        save(strcat([saveto_filepath '/IBM_reef_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '.mat']))
    elseif eflag==4
        save(strcat([saveto_filepath '/IBM_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '.mat']))
    end

    save(strcat([saveto_filepath '/dtime_avg.mat']),"dtime_avg")

%-----SIMULATE------------------------------------------------------------%