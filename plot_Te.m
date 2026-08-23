%% Plot electron temperatures for non-Maxwellian distributions
% Plot figures for the following paper: non-Maxwellian electron distribution in the D
% region during artificial heating (Paper I): Model development and
% electron temperature

%% Graphical setup
font_a  = 12.0; % Font size for tics
font_lb = 14.5; % Font size for label
font_lg = 11.5; % Font size for legends
font_t  = 15;   % Font size for title

%% Sort data according to altitude and E-field strength
cd(cooling_rate_data_directory); 
dfiles = dir('*.mat');
load(dfiles(1).name)

% Extract height and electric field from files
get_z0_E0 = cell(numel(dfiles),2); 
for i1 = 1:numel(dfiles)
    load(dfiles(i1).name)
    my_filename = dfiles(i1).name; 
    [~,N,~] = fileparts(my_filename); 
    a1 = strsplit(N,'-');
    get_z0_E0{i1,1} = a1{3}; 
    get_z0_E0{i1,2} = a1{5}; 
end 

% Convert from cell to double
get_z0_E0_new = zeros(size(get_z0_E0)); 
for k = 1:size(get_z0_E0,1)
    get_z0_E0_new(k,1) = str2double(get_z0_E0{k,1}); 
    get_z0_E0_new(k,2) = str2double(get_z0_E0{k,2});
end 

% Sort files 
[sort_dfiles,ind_sort] = sortrows(get_z0_E0_new); 
dfiles_cell = struct2cell(dfiles)'; 
dfiles_new1 = cell(size(dfiles,1),6); 
for i1 = 1:numel(dfiles)
    index = find(i1 == ind_sort); 
    dfiles_new1{index,1} = dfiles_cell{i1,1}; 
    dfiles_new1{index,2} = dfiles_cell{i1,2}; 
    dfiles_new1{index,3} = dfiles_cell{i1,3}; 
    dfiles_new1{index,4} = dfiles_cell{i1,4}; 
    dfiles_new1{index,5} = dfiles_cell{i1,5}; 
    dfiles_new1{index,6} = dfiles_cell{i1,6};  
end 

% Convert from cell to struct
fields = ["name", "folder", "date", "bytes", "isdir", "datenum"];
for i1 = 1:numel(dfiles)
    dfiles_new2 = cell2struct(dfiles_new1,fields,2);
end 

%% Find height z0 and electric field strength E0 in data 
z0 = unique(sort_dfiles(:,1)); 
E0_field_height = zeros(length(z0),size(sort_dfiles,1)/length(z0)); 
for l = 1:length(z0)
    for k = 1:size(sort_dfiles,1)/length(z0)
        n = l-1; 
        m = k + size(sort_dfiles,1)/length(z0)*n; 
        E0_field_height(l,k) = sort_dfiles(m,2); 
    end 
end

%% Index at different heights 
% Rows are height index
% Columns are electric field index
ind = zeros(size(E0_field_height,2),length(z0))'; 
v = 0; 
for j = 1:length(z0)
    aa = 1+v; 
    ind(j,:) = aa:aa+size(E0_field_height,2)-1; 
    v = v+size(E0_field_height,2); 
end 

%% Load data (go first to the folder containing the data...)

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;  

% Non-Maxwellian electron temperature 
Te_2nd = zeros(1,length(z0)*size(E0_field_height,2));
Te_eff = zeros(length(Te_grad),length(z0)*size(E0_field_height,2));
for i1 = 1:numel(dfiles)

    load(dfiles_new2(i1).name)
    
    % Non-Maxwellian temperature (second moment) from electron distribution function 
    Te_2nd(i1) = Te_nm_2nd;

    % Gradient temperature
    Te_eff(:,i1) = Te_grad; 

    % Energy/speed for Te_grad: constant delta E 
    res = length(Te_grad)+1;  % Resolution
    E_bins_new = (linspace(0,eVtoJ*max_E,res))'; % Bin edges [J] 
    E_bin_mid_eV_height = (E_bins_new(1:end-1)/2 + E_bins_new(2:end)/2)*(1/eVtoJ); % [eV]

end 

%% Plot gradient temperature with seconed moment temperature

% Find electric field of 2.0 V/m for plotting 
[~,my_ind_E0_2] = find(E0_field_height==2); 

fig = figure; 
subplot(2,1,1)
plot(E0_field_height(1,:),Te_2nd(ind(1,:)),'--','LineWidth',2.0)
hold on; 
plot(E0_field_height(2,:),Te_2nd(ind(2,:)),'-.','LineWidth',2.0)
plot(E0_field_height(3,:),Te_2nd(ind(3,:)),':','LineWidth',2.0)
plot(E0_field_height(4,:),Te_2nd(ind(4,:)),'LineWidth',2.0)
hold off; 
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.008 0.005])
xlim([0.75 4])
legend('70 km','80 km','90 km','100 km',...
      'Location','northwest',...
      'FontSize',font_lg)
xlabel('E_0 [V/m]','FontSize',font_lb)
ylabel('T_{2nd} [K]','FontSize',font_lb)
title('(a)','FontSize',font_t)

subplot(2,1,2)
plot(E_bin_mid_eV_height,Te_eff(:,ind(1,my_ind_E0_2(1))),'--','LineWidth',2.0)
hold on; 
plot(E_bin_mid_eV_height,Te_eff(:,ind(2,my_ind_E0_2(2))),'-.','LineWidth',2.0)
plot(E_bin_mid_eV_height,Te_eff(:,ind(3,my_ind_E0_2(3))),':','LineWidth',2.0)
plot(E_bin_mid_eV_height,Te_eff(:,ind(4,my_ind_E0_2(4))),'LineWidth',2.0)
hold off; 
set(gca,...
    'fontsize',font_a,...
    'box','on',...
    'tickdir','both',...
    'TickLength',[0.008 0.005])
ylabel('T_{eff} [K]','FontSize',font_lb)
xlabel('Energy [eV]','FontSize',font_lb)
title('(b)','FontSize',font_t)
legend('70 km','80 km','90 km','100 km',...
       'Location','northeast',...
       'FontSize',font_lg)
ylim([0 12000])
yticks([0 2500 5000 7500 10000])
yticklabels({'0','2500','5000','7500','10000'})
xlim([0 2.0])
colororder('sail')

% fig.PaperPositionMode = 'manual';
% orient(fig,'landscape')
% print(fig,'fig06.pdf','-dpdf')

