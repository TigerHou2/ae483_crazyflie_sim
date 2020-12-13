function [data] = simulate(h_num,p,xe,x0,xTol,t0,tStab,tMax,data)
%
% Inputs:
%   h_num - matlabFunction of EOMs created by findeoms.m
%       p - structure containing drone parameters
%       m - [kg] scalar,         mass
%       g - [m/s^2] scalar,      gravity
%       J - [kg-m^2] 3x3 matrix, moment of inertia
%      td - [s] scalar,          time delay of motors
%      sd - 12x1 vector,         sensor standard deviations
%      xe - 12x1 vector,         target state
%      x0 - 12x1 vector,         initial state
%    xTol - [m] scalar,          position tolerance for convergence
%      t0 - [s] scalar,          initial time
%   tStab - [s] scalar,          hover time
%    tMax - [s] scalar,          max sim time
%  tDelay - [s] scalar,          time delay on motor response

% Choose a sample rate. This is the number of times per second that the
% controller will run. This is also the number of times per second that
% data will be collected (at least, in the simulation I've created here -
% as you know, things are different in the hardware experiments, where the
% onboard controller and the offboard data logger run at different rates).
% requires all lower rates to be a factor of the highest rate!
posRate = 100; % PID controller freq for position
velRate = 100; % PID controller freq for velocity
angRate = 500; % PID controller freq for attitude
rtnRate = 500; % PID controller freq for attitude rate
pidRates = [posRate,velRate,angRate,rtnRate];
physicsRate = 2; % simulate physics at a multiple of the max data rate
sampleRate = max(pidRates) * physicsRate;
pidRates(end+1) = sampleRate;

commandRate = 20; % the rate of Python commands during flight.
                  % the command rate determines how often position accuracy
                  % is checked for the purpose of waypoint hovering

% Compute the time step. This is the length of time between one sample and
% the next.
tStep = 1 / sampleRate;
sampleNumber = ceil(tMax*sampleRate);

% Define variables to keep track of the time, state, and input during the
% simulation. I like to put these variables all together in a struct - that
% makes it easy to store and pass the data around, and also to add more
% things to the data later on if I want (rotor spin rates, for example).
%
% I initialize the time and state with their initial values. I initialize
% the input as an empty matrix - no inputs have been chosen yet.
if isempty(data)
    data = struct('t', t0, ...
                  'x', x0, ...
                  'u', [0; 0; 0; p.m * p.g], ...
                  'm', zeros(4,1), ... % motor PWM
                  'e', zeros(12,1), ...
                  'int', zeros(12,1), ... % integral error
                  'xd', zeros(12,1), ...
                  'xn', x0, ...
                  'lpf', {cell(12,1)});
    idx0 = 1;
else
    idx0 = length(data.t);
end
          
% generate noise matrix
noiseMat = normrnd(0,repmat(p.sd,1,sampleNumber));

% get iterators for holding data at each timestep
int = data.int(:,end);
LPF = data.lpf;

% Loop through all time steps.
for i = idx0:idx0+sampleNumber-1
    
    % Grab the precomputed sensor noise for the current time.
    noise = noiseMat(:,i-idx0+1);
    
    % The "current time" is at time step "i". So, the current time and
    % state are in the i'th column of data.t and data.x.
    t = data.t(:, i);
    x = data.x(:, i);
    m = data.m(:, i);
    xn = data.x(:, i) + noise;
    
    % Choose input based on the current time and state.
    des_pos = xe(1:3);
    [u,e,int,m,xd,LPF] = ...
        docontrol( i, des_pos, data.xn, data.e, int, ...
                   m, p.td, pidRates, LPF);
    
    % append data to structure.
      data.e( :, i + 1) = e;
    data.int( :, i + 1) = int;
      data.u( :, i + 1) = u;
      data.m( :, i + 1) = m;
     data.xd( :, i + 1) = xd;
     data.xn( :, i + 1) = xn;
    data.lpf = LPF;
    
    % Numerically integrate equations of motion to compute what the next
    % state will be given the current state and input.
    %
    %   @(t, x) h_num(x, u)     This is a anonymous function that ode45
    %                           will call to find xdot as a function of t
    %                           and x. The syntax means that when ode45
    %                           calls this anonymous function with (t, x),
    %                           this function will return the value
    %                           h_num(x, u). Note that the "x" comes from
    %                           ode45 and is changing all the time, while
    %                           the "u" is what we computed above, from the
    %                           controller, and is constant over each time
    %                           step. For help on anonymous functions, do:
    %
    %                               doc anonymous
    %
    %   [t, t + tStep]          This is the time interval to simulate. The
    %                           state (i.e., current) time is t, the final
    %                           (i.e., next) time is t + tStep.
    %
    %   x                       This is the state at the state of the time
    %                           interval that we want to simulate. In other
    %                           words, it is the current state.
    %
    [t_sol, x_sol] = ode45(@(t, x) h_num(x, u), [t, t + tStep], x);
    
    % Parse the solution that is returned by ode45 and store the next time
    % and state in our data structure. These will become the *current* time
    % and state the next time through our loop.
    %
    % Here is what ode45 returns:
    %
    %   t_sol       m x 1 matrix of times
    %   x_sol       m x n matrix of states, where the state time t_sol(i)
    %               is given by the 1 x n matrix x_sol(i, :)
    %
    % So, the "next time" is the last element of t_sol, and the "next
    % state" is the last row of x_sol. I always prefer to represent the
    % state as an n x 1 matrix instead of as a 1 x n matrix, so I take the
    % transpose of the result before storing it in data.
    data.t(:, i + 1) = t_sol(end, :);
    data.x(:, i + 1) = x_sol(end, :)';
    
    % check if convergence is achieved
    if i >= sampleRate * tStab ...
    && mod(i, ceil(physicsRate/commandRate)) == 0 ...
    && max(vecnorm(data.e(1:3,i+1-sampleRate*tStab:i+1),2,1)) < xTol
        break
    end
    
end

end