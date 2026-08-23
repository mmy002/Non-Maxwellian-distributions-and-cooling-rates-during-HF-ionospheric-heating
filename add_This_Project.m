[S1,S2,S3] = fileparts(which('add_This_Project'));

addpath(fullfile(S1),                                         '-end')
addpath(fullfile(S1,'Cross_section'),                         '-end')
addpath(fullfile(S1,'Cross_section','Process_data'),          '-end')
addpath(fullfile(S1,'Functions'),                             '-end')
addpath(fullfile(S1,'Fun_plot'),                              '-end')
addpath(fullfile(S1,'Cross_section/Process_data','Functions'),'-end')
addpath(fullfile(S1,'Cross_section','Data'),                  '-end')
addpath(fullfile(S1,'Boltzmann_solver'),                      '-end')
addpath(fullfile(S1,'Cross_section','Funtion_interp'),        '-end')
addpath(fullfile(S1,'Fun_plot','cm_and_cb_utilities'),        '-end')
addpath(fullfile(S1,'MC_simulation','MC_data','Histrogram'),  '-end')
addpath(fullfile(S1,'Cross_section','Data','New_xs_rot'),     '-end')
addpath(fullfile(S1,'Cooling_rates','Functions'),             '-end')

cooling_rate_data_directory = [S1,'/../Data/Cooling-rate-run-UT-20260222-182632']; 
