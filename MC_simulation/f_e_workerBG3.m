function [t_all,E_all,v_pars,v_perps,v_x,v_y] = f_e_workerBG3(n_colls,B,Eo1_0,w_HF,theta2B,T_e,Nn,dE,ps_data)
%         1     2     3      4       5                              1 2     3    4       5   6   7     8  9 
% f_e_worker - calculate the evolution of electron-velocities in HF-heating
% of magnetized plasma. This version only takes idealized ellastic
% collisions into account with a somewhat naive approach to angular
% scattering.
% 
% Calling:
%  [t_all,E_all,v_pars,v_perps,v_x,v_y] = f_e_worker2(n_colls,B,Eo1_0,w_HF,theta2B,w_e,T_e.Sigma,dE)
% Input:
%  n_colls - number or collisions to run for, scalar integer
%  B       - magnetic field strength (T), double scalar
%  Eo1_0   - local E-field amplitude of the HF-radio-wave (V/m), double
%            scalar
%  w_HF    - angular frequency of the HF-radio-wave heating wave (rad/s),
%            double scalar
%  theta2B - Angle to the magnetic field of the HF-radio-wave propagation
%            (radians), double scalar.
%  w_e     - angular electron gyro-frequency (rad/s), double scalar
%  T_e     - electron temperature (K), double scalar
% Output:
%  t_all   - time at collisions (s), double array [1 x n_colls + 1]
%  E_all   - electron energy at collision (eV), double array [1 x n_colls + 1]
%  v_pars  - electron-velocity (m/s) component parallel to B, double array
%            [1 x n_colls + 1] 
%  v_perps - electron-velocity (m/s) component perpendicular to B, double
%            array [1 x n_colls + 1] 
%  v_x     - electron-velocity (m/s) component perp to B par E-W, double
%            array [1 x n_colls + 1] 
%  v_y     - electron-velocity (m/s) component perp to B par N-S, double
%            array [1 x n_colls + 1] 

n_p = 1;
iP = 1;

% Allocate memory for outputs 
t_all = zeros(n_p,n_colls+1);  
E_all = zeros(n_p,n_colls+1);  
v_pars = zeros(n_p,n_colls+1);  
v_perps = zeros(n_p,n_colls+1); 
v_x = zeros(n_p,n_colls+1);     
v_y = zeros(n_p,n_colls+1);     

% Physical constants
m_e = 9.1093835611e-31;   % electron rest mass [kg]
q_e = 1.602176620898e-19; % elementary charge [C]

% Initial velocity, from a T_e K thermal population
v0 =  v_thermal(T_e,m_e)/sqrt(3)*randn(3,1);
while isnan(v0==1) 
    v0 =  v_thermal(T_e,m_e)/sqrt(3)*randn(3,1); % In case v0 is NaN
end 
v_long = v0;

% Intial time 
t0 = 0;
t_long = t0;

% Corresponding electron energy 
E0 = m_e/2*sum(v0.^2)/q_e;
E = E0; 

% Convert from [eV] to [J] 
eVtoJ = 1.602176634e-19;  

for i1 = 1:n_colls

    % Cross sections 
    xs_e_N2elastic = sigma_N2_el(E*eVtoJ,ps_data{1}(:,1),ps_data{1}(:,2));
    xs_e_O2elastic = sigma_O2_el(E*eVtoJ,ps_data{2}(:,1),ps_data{2}(:,2));
    xs_e_Oelastic = sigma_O_el_2(E*eVtoJ,ps_data{3}(:,1),ps_data{3}(:,2));
    xs_e_N2vib0_1 = sigma_N2_vib_01(E*eVtoJ,ps_data{4}(:,1),ps_data{4}(:,2)); 
    xs_e_N2vib0_2 = sigma_N2_vib_02(E*eVtoJ,ps_data{5}(:,1),ps_data{5}(:,2));
    xs_e_N2vib0_3 = sigma_N2_vib_03(E*eVtoJ,ps_data{6}(:,1),ps_data{6}(:,2)); 
    xs_e_N2vib0_4 = sigma_N2_vib_04(E*eVtoJ,ps_data{7}(:,1),ps_data{7}(:,2));
    xs_e_N2vib0_5 = sigma_N2_vib_05(E*eVtoJ,ps_data{8}(:,1),ps_data{8}(:,2)); 
    xs_e_N2vib0_6 = sigma_N2_vib_06(E*eVtoJ,ps_data{9}(:,1),ps_data{9}(:,2));
    xs_e_N2vib0_7 = sigma_N2_vib_07(E*eVtoJ,ps_data{10}(:,1),ps_data{10}(:,2)); 
    xs_e_N2vib0_8 = sigma_N2_vib_08(E*eVtoJ,ps_data{11}(:,1),ps_data{11}(:,2)); 
    xs_e_N2rot0_2 = sigma_N2_rot_02(E*eVtoJ,ps_data{12}(:,1),ps_data{12}(:,2)); 
    xs_e_N2rot0_4 = sigma_N2_rot_04(E*eVtoJ,ps_data{13}(:,1),ps_data{13}(:,2));
    xs_e_N2rot0_6 = sigma_N2_rot_06(E*eVtoJ,ps_data{14}(:,1),ps_data{14}(:,2)); 
    xs_e_N2rot0_8 = sigma_N2_rot_08(E*eVtoJ,ps_data{15}(:,1),ps_data{15}(:,2)); 
    xs_e_O2a1Dg = sigma_O2_a1_delta_g(E*eVtoJ,ps_data{16}(:,1),ps_data{16}(:,2)); 
    xs_e_O2b1Sgp = sigma_O2_b1_delta_g(E*eVtoJ,ps_data{17}(:,1),ps_data{17}(:,2)); 
    xs_e_Ofine_1_0 = sigma_O_fine_str_10(E*eVtoJ,ps_data{18}(:,1),ps_data{18}(:,2));  
    xs_e_Ofine_2_0 = sigma_O_fine_str_20(E*eVtoJ,ps_data{19}(:,1),ps_data{19}(:,2)); 
    xs_e_Ofine_2_1 = sigma_O_fine_str_21(E*eVtoJ,ps_data{20}(:,1),ps_data{20}(:,2));  
    xs_e_O3P1D = sigma_O_3P_1D(E*eVtoJ,ps_data{21}(:,1),ps_data{21}(:,2)); 
    xs_e_O3P1S = sigma_O_3P_1S(E*eVtoJ,ps_data{22}(:,1),ps_data{22}(:,2)); 
    xs_e_O2_vib = sigma_O2_vib(E*eVtoJ,ps_data{23}(:,1),ps_data{23}(:,2)); 
    sigma = [xs_e_N2elastic; xs_e_O2elastic; xs_e_Oelastic; xs_e_N2vib0_1; xs_e_N2vib0_2; xs_e_N2vib0_3; xs_e_N2vib0_4; ...
    xs_e_N2vib0_5; xs_e_N2vib0_6; xs_e_N2vib0_7; xs_e_N2vib0_8; xs_e_N2rot0_2; xs_e_N2rot0_4; xs_e_N2rot0_6; xs_e_N2rot0_8; ...
    xs_e_O2a1Dg; xs_e_O2b1Sgp; xs_e_Ofine_1_0; xs_e_Ofine_2_0; xs_e_Ofine_2_1; xs_e_O3P1D; xs_e_O3P1S; xs_e_O2_vib];

    % Speed [m/s]
    v = sqrt(sum(v0.^2)); 

    % Mean time between collisions [s]
    T_colls = exprnd(1./(sigma.*Nn*v)); 
    
    % Choose the three collisions with the lowest mean time 
    [T_colls1,iDE] = sort(T_colls);
    t_curr = [0,T_colls1(1)];
    if ~isempty(t_long) 
    t_curr = t_curr + t_long(end);
    t0 = t_curr(1);
    end
    
    % Equation of motion for E-field and B-field. Velocity as function of time
    v0_before = v0; 
    v_curr = v_e_HF_L(B,Eo1_0,v0,m_e,q_e,t_curr([1 end])-t0,-t0,theta2B,w_HF); 
    v0 = v_curr(:,end);
    
    % Here we handle angular scattering of electrons after collisions
    e_post = randn(3,1); % Velocity-unit-vector (new direction of velocity)
    e_v = e_post/norm(e_post); 
  
    % Energy before and after collision 
    E0 = m_e/2*sum(v0.^2)/q_e; 
    E  = E0 - dE(iDE(1));
    
    % If energy is negative 
    if E < 0 
        t_curr = [0,T_colls1(2)];
        if ~isempty(t_long) 
        t_curr = t_curr + t_long(end);
        t0 = t_curr(1);
        end
        
        % Beveglessligning for E-felt og B-felt. Velocity as function of time
        v_curr = v_e_HF_L(B,Eo1_0,v0_before,m_e,q_e,t_curr([1 end])-t0,-t0,theta2B,w_HF); 
        v0 = v_curr(:,end); 

        E0 = m_e/2*sum(v0.^2)/q_e; 
        E  = E0 - dE(iDE(2)); 
    end 

    if E < 0 
        t_curr = [0,T_colls1(3)];
        
        if ~isempty(t_long) 
        t_curr = t_curr + t_long(end);
        t0 = t_curr(1);
        end
        
        % Beveglessligning for E-felt og B-felt. Velocity as function of time
        v_curr = v_e_HF_L(B,Eo1_0,v0_before,m_e,q_e,t_curr([1 end])-t0,-t0,theta2B,w_HF); 
        v0 = v_curr(:,end); 

        E0 = m_e/2*sum(v0.^2)/q_e; 
        E  = E0 - dE(iDE(3)); 
    end 

    % New electron velocity after collision
    v0 = e_v*sqrt(2*q_e*E/m_e);
    
    % Save time and velocity
    t_long = [t_long,t_curr(end)];
    v_long = [v_long,v_curr(:,end)];

end

% Calculate the electron energy and extract the velocity-components
E_all(iP,:) = m_e/2*sum(v_long.^2)/q_e;
v_pars(iP,:) = v_long(3,:);
v_perps(iP,:) = (v_long(1,:).^2+v_long(2,:).^2).^.5;
v_x(iP,:) = v_long(1,:);
v_y(iP,:) = v_long(2,:);
t_all(iP,:) = t_long;

end