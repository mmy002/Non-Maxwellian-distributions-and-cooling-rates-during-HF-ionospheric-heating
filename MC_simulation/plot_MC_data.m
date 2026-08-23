%% Script for loading, processing and ploting data of MC runs 
% Format MC data: 10 000 electrons (row) x (301) n_colls+1 (columns)

%% Graphics setup
lftsz = 14; % font-size for labels
tftsz = 15; % font-size for titles
aftsz = 12; % font-size for axes and legends

%% Physical constants   
m_e   = 9.1093835611e-31;      % Electron rest mass [kg]
q_e   = 1.602176620898e-19;    % Elementary charge [C]
eVtoJ = 1.602176634e-19;       % Convert from [eV] to [J] 

%% Energy and velocity-distributions
% Here we reduce the particle energies and velocities down to energy and
% velocity-distributions - this is effectively histograms in 1 and 2
% dimensions. Depending on how many particles we trace and for how long the
% distributions will become increasingly smooth 

% Load data for the MC run (go first to the folder containing the data...)
dfiles = dir('*.mat');

% Bin-edges for the histcounts and histcounts2 function calls
E_bins = 0:0.025:4;

% Load the first file 
load(dfiles(1).name)

N_E = zeros(numel(E_bins)-1,size(E_all,2));
for it = size(E_all,2):-1:1
  N_E(:,it) = histcounts(E_all(:,it),E_bins)';
end

for i1 = 2:numel(dfiles)

  if rem(i1,10) == 0
    fprintf('%d at %s\n',i1,datestr(now,'HH:MM:SS'))
  end

  load(dfiles(i1).name)

  for it = 1:size(E_all,2)
    N_E(:,it) = N_E(:,it) + histcounts(E_all(:,it),E_bins)';
  end

end

%% Plot energy distribution and energy as a function of time 

% Energy/speed in the middle fo the bin
E_bin_mid = E_bins(1:end-1)/2 + E_bins(2:end)/2; 

% Pick certain times to plot
pick_ind_N_E = [1 2 8 26:length(avg_t)]; 
pick_ind_legend = [1 2 8 26:55:length(avg_t)]; 

[p1,p2] = fileparts(pwd); 
figure('name',p2)
subplot(2,1,1)
phNE = semilogy(E_bin_mid,N_E(:,pick_ind_N_E));
set(phNE,'linewidth',2)
xlabel('Energy [eV]','fontsize',lftsz)
ylabel('[#e]','fontsize',lftsz)
title('a)')
colormap(turbo)
legend(num2str(avg_t(pick_ind_legend)','%2.4f'))
set(gca,'fontsize',aftsz)
xlim([0 2.5])
subplot(2,1,2)
imagesc(avg_t,E_bins,log10(N_E))
xlabel('Time [s]','fontsize',lftsz)
ylabel('Energy [eV]','fontsize',lftsz)
title('b)')
set(gca,'fontsize',aftsz)



