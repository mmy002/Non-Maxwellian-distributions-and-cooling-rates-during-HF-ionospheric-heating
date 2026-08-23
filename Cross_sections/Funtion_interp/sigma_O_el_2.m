function cross_sect = sigma_O_el_2(E,E_old,xs_old)

% Find values in cross sections energy vector E_old closest to E energy vector 
% (E might be out of range of cross section energy vector)
[~,closestIndex_1] = min(abs(E-E_old(1))); 
[~,closestIndex_end] = min(abs(E-E_old(end))); 
E_interp = E(closestIndex_1:closestIndex_end);

%----------------------------------------------------------------------------------------------------------------
% INTERPOLATION
%----------------------------------------------------------------------------------------------------------------
sigma_interp = exp(interp1(E_old,log(xs_old),E_interp,'pchip')); 

%----------------------------------------------------------------------------------------------------------------
% EXTRAPOLATION
%----------------------------------------------------------------------------------------------------------------
sigma_extrap = exp(interp1(E_old,log(xs_old),E,'linear','extrap')); 
sigma_extrap = sigma_extrap(1:closestIndex_1-1); 

%----------------------------------------------------------------------------------------------------------------
% CROSS SECTIONS
%----------------------------------------------------------------------------------------------------------------
cross_sect = [sigma_extrap; sigma_interp; zeros(length(E)-closestIndex_end,1)]; 

end 