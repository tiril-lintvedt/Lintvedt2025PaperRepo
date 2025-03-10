function [simsamples, selected_candidates] = design_samples(species_composition,nslots, nsamples)
% 
% Simulates samples from blends of N species, with a given composition 
% , in a sample with Kcomponents/parts, and selects nsamples by using
% the Kennard-Stones algorithm.
% -------------------------------------------------------------------------
%
%       Input:
%                   species_composition - Saisir struct (d field : array K x N, i field:property name, v field: species name)
%                   nsamples            - number of samples to simulate
%
%
%       Output:
%                   simsamples          - Simulated samples (nsamples x K)
%                   selected_candidates - sample design (e.g [1 1 2 3 4] betyr 
%                                          at vi skal tilsette 2 deler av materiale 1, 
%                                           en del material 2, etc...)
%
%  
% -------------------------------------------------------------------------

Nspecies = size(species_composition.d, 2);
Kcomponents = size(species_composition.d, 1);

% Make candidate designs --------------------------------------------------

candidate_designs = nmultichoosek(1:Nspecies, nslots); % Make all possible compinations of sample with n slots from N species, with repitition 
ncandidates = size(candidate_designs,1);

candidate_compositions.d = zeros(ncandidates,Kcomponents);
candidate_compositions.v = species_composition.i;
candidate_compositions.i = [repelem('c',ncandidates,1),num2str([1:ncandidates]')];

for i = 1:ncandidates 
    design = candidate_designs(i,:);
    candidate_compositions.d(i,:) = [sum(species_composition.d(:,design),2)/nslots]' ;
    
end


% Stats of all candidate samples ------------------------------------------

cormap = corrcoef(candidate_compositions.d);
% Comp A vs comp B
scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/1.5 scrsz(4)-150])
hold on
set(gcf, 'Color', [1 1 1])
title('All candidate samples', 'Fontsize', 16)

c = 1;
for i = 1:Kcomponents          
    for j = 1:Kcomponents     
        if j <= i
            continue
        end
        subplot(3,ceil(Kcomponents/2),c)
        plot(candidate_compositions.d(:,i),candidate_compositions.d(:,j), '*' )
        text(0.05, 0.9, num2str(cormap(i,j)), 'Units', 'Normalized')
        xlabel(species_composition.i(i,:), 'Fontsize', 14)
        ylabel(species_composition.i(j,:), 'Fontsize', 14)
        c = c+1;
    end
end



% Select subset of candidate samples --------------------------------------

%ccompnorm = normalize(candidate_compositions.d);  % Normalized candidate compositions (avoid scaling effects?)
%[selected_samples,~] = KennardStone(ccompnorm,nsamples);
[selected_samples,~] = KennardStone(candidate_compositions.d,nsamples);
simsamples = candidate_compositions.d(selected_samples,:);
selected_candidates = candidate_designs(selected_samples,:);


% Stats of selected samples -----------------------------------------------

disp('Simulated component correlations: ---------------------------------')
cormap = corrcoef(simsamples)

% Comp A vs comp B
scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/1.5 scrsz(4)-150])
hold on
set(gcf, 'Color', [1 1 1])
title('Selected candidate samples', 'Fontsize', 16)
c = 1;
for i = 1:Kcomponents          
    for j = 1:Kcomponents     
        if j <= i
            continue
        end
        subplot(3,ceil(Kcomponents/2),c)
        plot(simsamples(:,i),simsamples(:,j), '*' )
        text(0.05, 0.9, num2str(cormap(i,j)), 'Units', 'Normalized')
        xlabel(species_composition.i(i,:), 'Fontsize', 14)
        ylabel(species_composition.i(j,:), 'Fontsize', 14)
        c = c+1;
    end
end


% Distribution of sample compositions
figure
set(gcf, 'color', [1 1 1])
hold on
for i = 1:Kcomponents
    subplot(2,ceil(Kcomponents/2),i)
    histogram(simsamples(:,i))
    title(species_composition.i(i,:))
end

% PCA on selected samples
[T,cum_var,var,mx,V,U] = pca(simsamples, 3);
X.d = simsamples;
X.i = num2str([1:size(simsamples,1)]');
plot_scores(T,X,1:size(X.i,2),1:size(X.i,2))

end