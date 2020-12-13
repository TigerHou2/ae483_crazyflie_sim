%% docontrol.m
function [U,e,int,m,xd,LPF] = ...
    docontrol(idx,des_pos,XN,E,INT,m0,tDelay,pidRates,LPF)
%DOCONTROL implements the cascaded PID controller from Crazyflie
%
% Inputs:
%       idx - integer, the index of the current timestep
%   des_pos - 3x1 matrix of desired position
%        XN - 12xN matrix of perturbed states over N timesteps
%         E - 12xN matrix of state errors
%       INT - 12x1 vector of integral errors
%        m0 - 12x1 vector of motor PWM values at the (N-1)th timestep
%     noise - 12x1 vector of sensor noise
%    tDelay - scalar, motor time delay from 0% to 100% thrust
%             we assume that the motor ramp-up time is linear; for exmaple,
%             going from 0% to 50% thrust will take half as long as
%             going from 0% to 100%.
%  pidRates - 5x1 vector of PID controller rates for
%             (1)position, (2)velocity, (3)attitude, (4)attitude rate
%             the last element (5) is the physics rate, which should be
%             an integer multiple of the highest controller rate
%       LPF - 12x1 cell of low pass filter data structures
%
% Outputs:
%         U - 4x1 vector of controller parameters
%               ([3x1] moment about x-, y-, and z-axis, and thrust)
%         e - 12xN matrix of state errors (updated for current timestep)
%       int - updated integral errors
%         m - 4x1 vector of motor outputs
%        xd - 12x1 vector of desired states
%       LPF - updated low pass filters
% 
% Note - the 12 states are as follows:
%   [3x1] position in fixed frame
%   [3x1] Euler angle (yaw, pitch, roll)
%   [3x1] velocity in body frame
%   [3x1] Eular angle rate (yaw, pitch, roll)

%% frame conversion
% convert between sim axis definition and Crazyflie axis definition
% because the sim's z-axis is pointed to ground, while the Crazyflie's
% z-axis is pointed upward
%
% note that we also adopt the ZYX Euler angle convention for the attitude
% whereas other states use the XYZ convention

XN = XN .* [  1, -1, -1, ... position
             -1, -1,  1, ... attitude
              1, -1, -1, ... velocity
              1, -1, -1, ... attitude rate
              ]';

%% calculate physics rate
dt = 1 / max(pidRates);

%% constants
thrustScale = 1000;
thrustBase  = 36000;
thrustMin   = 20000;
thrustCap   = 65535; % unsigned 16-bit int signal to motor, equals 60 grams
signed_int16 = 32767; % max value of 16-bit signed int
L = 0.092; % m, width of Crazyflie (between adjacent motors)
D = 0.045; % m, diameter of rotor
Cp = 0.10 * 1.5; % assume const. rotor power coeff. (0.1 from paper)
Ct = 0.15 * 1.5; % assume const. rotor thrust coeff. (0.15 from paper)
rho = 1.225; % kg/m^3, density of air
rpm = @(PWM) 0.2685 * PWM + 4070.3; % pwm to motor rpm
thr = @(PWM) Ct * rho * (rpm(PWM)/60)^2 * D^4; % thrust per motor
tau = @(PWM) Cp * rho * (rpm(PWM)/60)^2 * D^5 / (2*pi); % torque
idx = idx + 1;
yawDecay = 7.5; % deg/s, models the yaw decay of the Kalman filter

%% PID defaults
%   position_controller_pid.c for pos and vel gains
%   
%   12x3 matrix of proportional, integral, and derivative gains
PID = [  2.00    0.00    0.00 ; ... % x
         2.00    0.00    0.00 ; ... % y
         2.00    0.50    0.00 ; ... % z
         6.00    1.00    0.35 ; ... % yaw
         6.00    3.00    0.00 ; ... % pitch
         6.00    3.00    0.00 ; ... % roll
        25.00    1.00    0.00 ; ... % vx
        25.00    1.00    0.00 ; ... % vy
        25.00   15.00    0.00 ; ... % vz
       250.00  500.00    2.50 ; ... % roll rate
       250.00  500.00    2.50 ; ... % pitch rate
       120.00   16.70    0.00 ; ... % yaw rate
       ];

%% PID modifications
%   any changes to the default PID values should be made here
PID(1,:) = [2.2, 0, -0.05];
PID(2,:) = [2.2, 0, -0.05];
PID(3,:) = [2.6, 0.4, 0.05];
% PID(1,:) = [3.3, 0.1, 0.1];
% PID(2,:) = [3.6, 0.1, 0.1];
% PID(3,:) = [1.9, 0.2, 0.1];
% PID(4:6,:) = PID(4:6,:) * 2;
% PID(10:12,:) = PID(10:12,:) * 2;

%% integration limits
%   pid.h
%   12x1 vector of integration limits
ilim = [5000, 5000,  5000, ... position
         360,   20,    20, ... attitude
        5000, 5000,  5000, ... velocity
        33.3, 33.3, 166.7, ... attitude rate
       ]';

%% thrust limit
tlim = thrustCap/2/thrustScale;

%% control limits
%   position_controller_pid.c
%   12x1 vector of control limits
clim = [  Inf,  Inf,  Inf, ... position
          Inf,   22,   22, ... attitude
          1.1,  1.1,  1.1, ... velocity
          Inf,  Inf,  Inf, ... attitude rate
        ]';
% clim(7:9) = clim(7:9) * 0.7;

%% low pass filter cutoff frequencies
%   position_controller_pid.c
%   attitude_controller_pid.c
cutoff_freq = [ 20, 20, 20, ... position
                15, 15, 15, ... attitude
                20, 20, 20, ... velocity
                30, 30, 30, ... attitude rate
                ]';
            
%% initialize return values
xd = nan(12,1);  % desired states
int = nan(12,1); % integral gains

%% Position PID
sub = 1:3;
sub_next = 7:9;
xd(sub) = des_pos;
dT = 1/pidRates(1); % timestep for the position PID
IDX = idx - 1 - mod(idx-2,dT/dt);
IDXm = max(IDX-dT/dt,1);
E(sub,idx) = des_pos-XN(sub,IDX);
integral_terms = INT(sub) + E(sub,IDX).*dT * double(IDX == idx-1);
integral_terms = constrain(integral_terms,ilim(sub));
int(sub) = integral_terms;
deriv_term = (E(sub,IDX)-E(sub,IDXm)) ./ dT;
    [LPF{sub(1)},d1] = lpf(LPF{sub(1)}, deriv_term(1), pidRates(1), cutoff_freq(sub(1)));
    [LPF{sub(2)},d2] = lpf(LPF{sub(2)}, deriv_term(2), pidRates(1), cutoff_freq(sub(2)));
    [LPF{sub(3)},d3] = lpf(LPF{sub(3)}, deriv_term(3), pidRates(1), cutoff_freq(sub(3)));
    deriv_term = [d1;d2;d3];
des_vel = runPid( E(sub,IDX), integral_terms, deriv_term, PID(sub,:) );
des_vel = constrain(des_vel,clim(sub_next));

%% Velocity PID
sub = sub_next;
sub_next = 4:6;
xd(sub) = des_vel;
dT = 1/pidRates(2); % timestep for the velocity PID
IDX = idx - 1 - mod(idx-2,dT/dt);
IDXm = max(IDX-dT/dt,1);
E(sub,idx) = des_vel-XN(sub,IDX);
integral_terms = INT(sub) + E(sub,IDX).*dT * double(IDX == idx-1);
integral_terms = constrain(integral_terms,ilim(sub));
int(sub) = integral_terms;
deriv_term = (E(sub,IDX)-E(sub,IDXm)) ./ dT;
    [LPF{sub(1)},d1] = lpf(LPF{sub(1)}, deriv_term(1), pidRates(2), cutoff_freq(sub(1)));
    [LPF{sub(2)},d2] = lpf(LPF{sub(2)}, deriv_term(2), pidRates(2), cutoff_freq(sub(2)));
    [LPF{sub(3)},d3] = lpf(LPF{sub(3)}, deriv_term(3), pidRates(2), cutoff_freq(sub(3)));
    deriv_term = [d1;d2;d3];
des_ang = runPid( E(sub,IDX), integral_terms, deriv_term, PID(sub,:) );
des_ang = constrain(des_ang,flipud([tlim;clim(sub_next(2:3))]));
roll = des_ang(1); % note: this is unintuitive since roll doesn't control x
pitch = des_ang(2); %      the reason this is setup like this is because
                    %      these are raw values - they will now be
                    %      converted to the actual roll/pitch desired.
                    %      (As per Crazyflie documentation)
thrust = des_ang(3);
des_ang(1) =   constrain(XN(4,IDX),yawDecay*dT); % yaw
des_ang(2) =   roll *cos(XN(4,IDX)) + pitch*sin(XN(4,IDX));  % pitch
des_ang(3) = - pitch*cos(XN(4,IDX)) + roll *sin(XN(4,IDX));  % roll
des_ang = constrain(des_ang,clim(sub_next));

%% Thrust
thrust = constrain(thrust,tlim);
thrust = thrust * thrustScale + thrustBase;
thrust = max(thrust,thrustMin); % this is the pwm of U(4)

%% Attitude PID
sub = sub_next;
sub_next = 10:12;
xd(sub) = des_ang;
dT = 1/pidRates(3); % timestep for the attitude PID
IDX = idx - 1 - mod(idx-2,dT/dt);
IDXm = max(IDX-dT/dt,1);
E(sub,idx) = des_ang-rad2deg(XN(sub,IDX));
E(4,idx) = mod(E(4,idx)+180,360)-180;
integral_terms = INT(sub) + E(sub,IDX).*dT * double(IDX == idx-1);
integral_terms = constrain(integral_terms,ilim(sub));
int(sub) = integral_terms;
deriv_term = (E(sub,IDX)-E(sub,IDXm)) ./ dT;
%     [LPF{sub(1)},d1] = lpf(LPF{sub(1)}, deriv_term(1), pidRates(2), cutoff_freq(sub(1)));
%     [LPF{sub(2)},d2] = lpf(LPF{sub(2)}, deriv_term(2), pidRates(2), cutoff_freq(sub(2)));
%     [LPF{sub(3)},d3] = lpf(LPF{sub(3)}, deriv_term(3), pidRates(2), cutoff_freq(sub(3)));
%     deriv_term = [d1;d2;d3];
des_rate = runPid( E(sub,IDX), integral_terms, deriv_term, PID(sub,:) );
des_rate = flipud(des_rate); % change to roll, pitch, yaw (x,y,z)
des_rate = constrain(des_rate,clim(sub_next));

%% Attitude Rate PID
% note: in the state setup the order of angles has changed:
%   previously it was yaw, pitch, roll
%   now it is roll, pitch, yaw (x,y,z)
sub = sub_next;
xd(sub) = des_rate;
dT = 1/pidRates(4); % timestep for the attitude rate PID
IDX = idx - 1 - mod(idx-2,dT/dt);
IDXm = max(IDX-dT/dt,1);
E(sub,idx) = des_rate-rad2deg(XN(sub,IDX));
integral_terms = INT(sub) + E(sub,IDX).*dT * double(IDX == idx-1);
integral_terms = constrain(integral_terms,ilim(sub));
int(sub) = integral_terms;
deriv_term = (E(sub,IDX)-E(sub,IDXm)) ./ dT;
%     [LPF{sub(1)},d1] = lpf(LPF{sub(1)}, deriv_term(1), pidRates(2), cutoff_freq(sub(1)));
%     [LPF{sub(2)},d2] = lpf(LPF{sub(2)}, deriv_term(2), pidRates(2), cutoff_freq(sub(2)));
%     [LPF{sub(3)},d3] = lpf(LPF{sub(3)}, deriv_term(3), pidRates(2), cutoff_freq(sub(3)));
%     deriv_term = [d1;d2;d3];
output = runPid( E(sub,IDX), integral_terms, deriv_term, PID(sub,:) );
output = constrain(output,signed_int16);

%% Motor Commands
r = output(1) / 2;
p = output(2) / 2;
y = -output(3);
m1 = max(min(thrust-r+p+y,thrustCap),0);
m2 = max(min(thrust-r-p-y,thrustCap),0);
m3 = max(min(thrust+r-p+y,thrustCap),0);
m4 = max(min(thrust+r+p-y,thrustCap),0);

%% Motor Delay
dm = m0 - [m1;m2;m3;m4];
dmCap = thrustCap * dt / tDelay;
dm = constrain(dm,dmCap);
mOut = m0 - dm;
m1 = mOut(1); m2 = mOut(2); m3 = mOut(3); m4 = mOut(4);

%% Moments and Thrust
t1 = tau(m1); t2 = tau(m2); t3 = tau(m3); t4 = tau(m4);
f1 = thr(m1); f2 = thr(m2); f3 = thr(m3); f4 = thr(m4);
U = zeros(4,1);
U(4) = (f1+f2+f3+f4);
U(1) = L/2 * (f3+f4-f1-f2); % roll
U(2) = L/2 * (f2+f3-f1-f4); % pitch
U(3) = (t1-t2+t3-t4); % yaw

e = E(:,idx);
m = mOut;

end % docontrol.m
