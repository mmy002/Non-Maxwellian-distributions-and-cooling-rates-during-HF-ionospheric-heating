function rv = ode_boris_moverV(t,r0v0,q,m,Ein,Bin,wMax,nu,v_loss_fraction)
% ODE_BORIS_MOVER - charged particle equation-of-motion-solver in B and E-fields
%   The ode_boris_mover integrates the equations of motion in
%   magnetic and electrical fields with the Boris-mover scheme
%   [e.g. 1].
% 
% Calling:
%   rv = ode_boris_mover(t,r0v0,q,m,Ein,Bin[,wMax])
% Input:
%   t    - time (s), double array [1 x n_t] for desired output
%   r0v0 - initial particle position, double array [6 x 1] with the
%          first 3 components for the position (m) and the last 3
%          positions for the velocity (m/s). The velocities will be
%          treated as the velocity at half the time-step before t0,
%          that is no special treatment of the velocities at the
%          initiation.
%   q    - particle charge (C), scalar double
%   m    - particle mass (kg), scalar double
%   Ein  - electrical field (V/m), either a [3 x 1] array for a
%          static uniform E-field or a function handle to a
%          function f_E(t,r) returning a [3 x 1] electrical in
%          position r at time t.
%   Bin  - magnetic field (T), either a [3 x 1] array for a
%          static uniform B-field or a function handle to a
%          function f_B(t,r) returning a [3 x 1] electrical in
%          position r at time t.
%   wMax - maximum expected gyro-frequency, used to limit the
%          time-step of the solver. Optional argument, defaults to
%          zero. 
% Output:
%  rv    - particle phase-space coordinates, double array [6 x n_t].
%          The first 3 rows are the particle position (m), rows 4
%          to 6 are the particle velocity (m/s),
% 
% The equations of the Boris scheme are:
%   x_{k+1}   = x_k + Dt*v_{k+1/2}
%   v_{k+1/2} = u' + q' E_k
%
% with
%
%   u' = u + ( u + ( u � h ) ) � s
%   u  = v_{k-1/2} + q'*E_k
%   h  = q'*B_k
%   s  = 2*h/( 1 + h^2 )
%
%   q' = Dt�(q/2m)
% 
% Example:
%   % Physical constants and parameters
%   m_e = 9.1094e-31;        % electron mass (kg)
%   q_e = 1.6022e-19;        % electron charge (C)
%   B = 5e-5;                % Magnetic field strength (T)
%   w_e = (q_e*B/m_e);       % electron gyro-frequency
%   T_gyro = 1/(w_e/(2*pi)); % gyro-period
%   % Initial conditions and time-span
%   v_0    = (2*q_e*1/m_e).^(1/2); % velocity of 1 eV electron (m/s)
%   r_gyro = v_0/w_e;               % electron gyro-radius
%   t_out = linspace(0,100*T_gyro,10001); % requested time-steps 
%   % ODE-integration:
%   rv = ode_boris_mover(t_out,[0;-r_gyro;0;v_0(1);0;0],-q_e,m_e,[0;0;0],[0;0;B]);
%   % Digestion of solution
%   r_c = mean(rv(1:3,:),2);
%   subplot(2,2,1)
%   plot(rv(1,:),rv(2,:)),
%   ylabel('y-position (m)'),xlabel('x-position (m)'),
%   title('X-Y trajectory 100 gyrations')
%   subplot(2,2,2)
%   plot(rv(1,:),rv(2,:)),
%   ylabel('y-velocity (m/s)'),xlabel('x-velocity (m/s)'),
%   title('v_x-v_y 100 gyrations')
%   subplot(2,2,3)
%   plot(t_out,((rv(1,:)-r_c(1)).^2+(rv(2,:)-r_c(2)).^2).^.5);
%   r_max = max(((rv(1,:)-r_c(1)).^2+(rv(2,:)-r_c(2)).^2).^.5);
%   r_min = min(((rv(1,:)-r_c(1)).^2+(rv(2,:)-r_c(2)).^2).^.5);
%   r_g = mean(((rv(1,:)-r_c(1)).^2+(rv(2,:)-r_c(2)).^2).^.5);
%   ylabel('(r_{max} - r_{min})/r_g'),xlabel('time (s)'),
%   title(sprintf('relative variation of gyro-radius: %g',...
%                 (r_max-r_min)/r_gyro))
%   subplot(2,2,4)
%   plot(t_out,(rv(4,:).^2+rv(5,:).^2+rv(6,:).^2)*m_e/2/q_e)
%   K_max = max((rv(4,:).^2+rv(5,:).^2+rv(6,:).^2)*m_e/2/q_e);
%   K_min = min((rv(4,:).^2+rv(5,:).^2+rv(6,:).^2)*m_e/2/q_e);
%   ylabel('(K_{max} - K_{min})/K_0'),xlabel('time (s)'),
%   title(sprintf('relative variation of Kinetic energy: %g',...
%                 (K_max-K_min)/1))
% 
%   The current implementation adjusts the time-increment to 1/20
%   of the gyro-period at the current particle position. That might
%   lead to errors where the B-field varies abruptly over distances
%   smaller than |1/20*r_gyro|. The mover conserves kinetic energy
%   for particles in a uniform B-field to machine precision, and
%   keeps relative variation of gyro-radius to ~4e-4.

%  Copyright � Bjorn Gustavsson 20190118, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

% TODO: implement handling of w_max similar to E and B so that we
% can get the maximum gyro-frequency inside a sphere with radius
% r_gyro centred at the current position.

if nargin < 7 || isempty(wMax)
  wMax = 0;
end
if nargin < 8 || isempty(wMax)
  nu = 0;
end
if nargin < 8 || isempty(wMax)
  v_loss_fraction = 0;
end

n_p = size(r0v0,1);
%% initial speeds
for i_p = n_p:-1:1,
  v0_norm(i_p) = norm(r0v0(i_p,4:6));
end
%% Start assigning output variables:
rv(n_p,6,numel(t)) = 0;
rv(:,:,1) = r0v0;

% wbh = waitbar(0,'Patience Bjorn, patience...');
% Initialize the particle position and velocity:
r = rv(:,1:3);
v = rv(:,4:6);
% Starting time:
t_running = t(1);
i_next = 2;
while t_running < t(end)
  
  % Calculate the E and B-fields at current positions:
  if isempty(Ein) % "Special case for when the E-field depends on
                  % the local B-field (case: plasma-spheric
                  % corotation) - then the Bin function has to
                  % return both B and E fields
    [B,E] = Bin(t_running,r);
  else
    if isa(Ein,'function_handle')
      E = Ein(t_running,r);
    else
      E = Ein;
    end
    if isa(Bin,'function_handle')
      B = Bin(t_running,r);
    else
      B = Bin;
    end
  end
  
  %% Update the increment in time
  %  to ensure nice trajectories along the track use 1 20th of a
  %  gyro-period, or the time to the next requested time for output
  %  
  %  Start with the gyro-frequency:
  w_c = max(abs(q.*sum(B.^2,2).^(1/2)./m));
  Dt = min(0.05/(max(w_c,wMax)/2/pi),t(i_next)-t_running);
  
  for i_p = 1:n_p,
    %% Then the Boris-scheme:
    q_prime = Dt*q(i_p)/(2*m(i_p));
    h = q_prime*B(min(i_p,end),:);
    s = 2*h/(1+norm(h)^2);
    u  = v(i_p,:) + q_prime*E(min(i_p,end),:);
    u_prime = u + cross( u + cross( u, h ), s);
    
    v_halfstep = u_prime + q_prime*E(min(i_p,end),:);
    r(i_p,:) = r(i_p,:) + v_halfstep*Dt;
    v(i_p,:) = v_halfstep;
    if rand(1) < nu*Dt
      e_new = randn(3,1);
      e_new = e_new/norm(e_new);
      %v(i_p,:) = v0_norm(i_p) + (norm(v_halfstep)-v0_norm(i_p))*e_new*(1-v_loss_fraction);    
      v(i_p,:) = norm(v_halfstep)*e_new*(1-v_loss_fraction);    
    end
  end
  t_running = t_running + Dt;
  if t_running == t(i_next)
    rv(:,:,i_next) = [r,v];
    if i_p == n_p
      i_next = i_next+1;
    end
  end
end
