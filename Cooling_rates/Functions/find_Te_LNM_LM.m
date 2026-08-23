function out = find_Te_LNM_LM(L_NM_tot,Te,Tn,nNe,nO,nN2,nO2)
% Compute electron temperature from Maxwellian cooling rates with
% the same total energy as non-Maxwellian
% L_NM_tot: Sum of non-Maxwellian cooling rates 
% L_M_tot: Maxwellian cooling rates, mosly from Pavlov analytical expressions

% Maxwellian cooling rates, mostly based on Pavlov analytical expression.
L_M_tot = cooling_rates_tables_pavlov_tot(Te,Tn,nNe,nO,nN2,nO2);  

% Output
out = L_NM_tot - L_M_tot;

end 