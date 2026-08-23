function [dI_sum,dIe,deg_Ie,ratio] = degrading_electrons_works3(ex_E,xs_i,E_bins,v_bin_mid,nN_i,num_e_bin_v)
% This function computes the inelastic collisions integral for electrons
% degrading in energy due to inelastic collisions
% Original code by Bjoern Gusavsson, modified by Margaretha Myrvang 
% INPUT:
% ex_E: [levels x 1] Excitation energy [J], levels = number of included
% cross sections
% xs_i: [energies x levels] Cross section as function of energy [m^2] 
% E_bins: [bin edges x 1] Energy bin edges [J], energies = bin edges-1
% v_bin_mid: [energies x 1] Speed in the middle of the bin [m/s]
% nN_i: [densities x 1] Density of neutral species [m^-3]
% num_e_bin_v: [energies x 1] Input histogram for speed [#e/m^3]
% OUTPUT: 
% dI_sum: [energies x 1] Sum over all different cross sections for gain 
% and loss for each bin [#e/(m^3*s)]
% dIe: [energies x levels] Number of electrons degrading from higher to lower
% energies [#e/(m^3*s)]
% deg_Ie: [energies x levels] Degraded electrons at lower energies
% distributed into their respective bins [#e/(m^3*s)]
% ratio: Ratio between electron loss and gain (should be 1)

%  Flux of electrons #e/(m^2*s)-> histogram of flux
Ie = v_bin_mid .* num_e_bin_v;  

% Matlab likes that you create empty matrices for for-loops
dIe           = zeros(length(v_bin_mid),length(ex_E)); 
deg_Ie        = zeros(length(v_bin_mid),length(ex_E)); 
dI_sum_per_xs = zeros(length(v_bin_mid),length(ex_E)); 
ratio         = zeros(1,length(ex_E)); 

% Width of energy bin, increases with energy
gradE = @(iE) E_bins(iE+1) - E_bins(iE); % [J]

% Loop over different excitation levels (ex_E) with corresponding cross sections (xs) and different energies in
% the middle of the bin 
for idE = 1:length(ex_E) 
  
  for i = length(E_bins)-1:-1:1 

    % Number of electrons degrading from higher to lower energies;
    % flux loss at higher energies to lower energies, i.e. how many electrons we loose
    % Units: m^-3 * #e/(m^2*s) * m^2 = #e/(m^3*s)
    dIe(i,idE) = nN_i(idE) * Ie(i) .* xs_i(i,idE); 
    
    % The flux of electrons degraded from energy bin [E(iE), E(iE)+gradE(iE)] to any lower energy bin 
    % by excitation of the idE-th state of the current species.
    Ie_degraded = (min(1,ex_E(idE)./gradE(i)))*dIe(i,idE); % Basically determines if dE (excitation energy) 
    % is larger/smaller than gradE and how much of Ie_degraded is outside Ie's bin (percent)
    
    % Find the energy bins that the electrons in the current energy bin will degrade to when losing dE_i(idE,1): 
    % Find index of energy (E_bins) between degraded energies: E(iE) - ex_E(1,idE) and E(iE) + gradE(i) - ex_E(1,idE) 
    i_upper_between = find( E_bins > E_bins(i) - ex_E(idE) & E_bins < E_bins(i) + gradE(i) - ex_E(idE) );
    
    % Find index of E_bins to the left of E(i) - ex_E(1,idE) (lower boundary)
    % First check if i_upper_between is non-empty or empty
    if isempty(i_upper_between) == 0 % If i_upper_between is not empty:
      
      if  i_upper_between(1) > 1 % If the index of i_upper_between = 2 or larger (i.e only important for the lowest energies)
        i_upper = zeros(length(i_upper_between)+1,1); % Create a new zero matrix that includes the new index
        i_upper(1) = i_upper_between(1)-1; % Now find the index of E_bins to the left of E(i) - ex_E(1,idE)
        i_upper(2:end) = i_upper_between; % Indexes between E_bins between E(i) - ex_E(1,idE) and E < E(i) + gradE(i) - ex_E(1,idE) 
      else 
        i_upper = i_upper_between; % If the index of the first element in i_upper_between = 1 (i.e only important for the lowest energies)
      end 
      
    else 
      i_upper = i_upper_between; % If i_upper_between is empty:
    end 
    
    % Fraction of electrons in a certain lower energy bin
    partition_fract = zeros(size(i_upper));
    
    if ( ~isempty(i_upper) && i_upper(1) < i ) 
      
      % Distribute the degrading electrons between those bins
      % Here partition_fract will contain all energy-bins electrons can
      % degrade into - including the current bin the electrons start at. 
      partition_fract(1) = min(1,(E_bins(i_upper(1)) + gradE(i_upper(1)) -  (E_bins(i) - ex_E(idE)))/gradE(i)); % gradE(i) correct?
      if length(i_upper)>2
        partition_fract(2:end-1) = min(1,gradE(i_upper(2:end-1))/gradE(i));
      end
      partition_fract(end) =  min(1,(E_bins(i) + gradE(i)-ex_E(idE) - E_bins(i_upper(end)) )/gradE(i));
      if i_upper(end) == i
        % Here we set the partition_fract for the current/starting bin to
        % zero - we are only interested in calculating the partitioning
	    % of the electrons degrading out of that bin
        partition_fract(end) = 0;
      end
      partition_fract = partition_fract/sum(partition_fract);
      
      % Distribute the degrading electrons in the correct bin 
      for i_u = find(partition_fract~=0) 
        deg_Ie(i_upper(i_u),idE) = ( deg_Ie(i_upper(i_u),idE) + Ie_degraded * partition_fract(i_u) ); % [#e/(m^3*s)]
      end
      
    end
    
  end 
  
  % Ratio between electron loss and gain (should be 1)
  ratio(idE) = abs(sum(dIe(:,idE))/sum(deg_Ie(:,idE)));
  
  % Set flux loss to be negative as it should [#e/(m^3*s)]
  dIe(:,idE) = -1*dIe(:,idE); 
  
  % Sum of electron loss and gain for each bin per cross section [#e/(m^3*s)]
  dI_sum_per_xs(:,idE) = dIe(:,idE) + deg_Ie(:,idE);  
  
end 

% Sum of electron loss and gain for each bin [#e/(m^3*s)]
dI_sum = sum(dI_sum_per_xs,2); 
 
