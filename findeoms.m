function [h, x, u, h_num, params] = findeoms()

% Define all 12 states as symbolic variables. The "real" tag at the end
% tells MATLAB to assume that all of these variables are real-valued (as
% opposed to compled-valued, for example).
syms o1 o2 o3 hy hp hr v1 v2 v3 w1 w2 w3 real

% For convenience, define the position (a point), the linear velocity (a
% vector), and the angular velocity (a vector) in coordinates, in terms of
% the state variables. Note that we do NOT group the yaw, pitch, and roll
% angles into a 3x1 matrix, to emphasize that these angles are NOT the
% coordinates of any point or vector.
o = [o1; o2; o3];
v = [v1; v2; v3];
w = [w1; w2; w3];

% For convenience, define the set of all states as a 12x1 matrix "x". The
% order in which we choose to list the states is up to us - we just need to
% be consistent everywhere else.
x = [o; hy; hp; hr; v; w];

% Define all 4 inputs as real-valued symbolic variables.
syms u1 u2 u3 u4 real

% For convenience, define the set of all inputs as a 4x1 matrix "u".
u = [u1; u2; u3; u4];

% Define all parameter values. 
m = 0.0286;             % kg, mass
g = 9.81;               % m/s^2, acceleration of gravity
td = 0.18;              % s, motor time delay from zero to full thrust

% moment of inertia calculations
    % battery
    m_batt = 0.0071;
    x_batt = 0.03;
    y_batt = 0.02;
    z_batt = 0.007;
    J_batt = diag([ 1/12 * m_batt * (y_batt^2+z_batt^2), ...
                    1/12 * m_batt * (x_batt^2+z_batt^2), ...
                    1/12 * m_batt * (x_batt^2+y_batt^2) ]);
    % motors
    m_motor = 0.0027;
    x_motor = 0.092 / 2;
    y_motor = 0.092 / 2;
    z_motor = 0;
    J_motor = diag([ m_motor * (y_motor^2+z_motor^2), ...
                     m_motor * (x_motor^2+z_motor^2), ...
                     m_motor * (x_motor^2+y_motor^2) ]) * 4;
    % flowdeck v2
    m_flow = 0.0016;
    x_flow = 0.021;
    y_flow = 0.028;
    z_flow = 0.004;
    J_flow = diag([ 1/12 * m_flow * (y_flow^2+z_flow^2), ...
                    1/12 * m_flow * (x_flow^2+z_flow^2), ...
                    1/12 * m_flow * (x_flow^2+y_flow^2) ]);
    % main board
    m_board = m - m_batt - 4*m_motor - m_flow;
    x_board = 0.092;
    y_board = 0.092;
    z_board = 0.003;
    J_board = diag([ 1/12 * m_board * (y_board^2+z_board^2), ...
                     1/12 * m_board * (x_board^2+z_board^2), ...
                     1/12 * m_board * (x_board^2+y_board^2) ]);
    % total
    J = J_batt + J_motor + J_flow + J_board;

% model sensor noise using a zero-mean normal distribution
sd = [0.005, 0.005, 0.005, ...      position
      deg2rad([0.2, 0.2, 0.2]), ... attitude
      0.020, 0.020, 0.020, ...      velocity
      deg2rad([3, 3, 3]), ...       attitude rate
      ]';

% It is helpful to package all parameter values into a struct that can be
% stored for use later in control design.
params = struct('m', m, 'g', g, 'J', J, 'td', td, 'sd', sd);

% Find applied force (from gravity and rotors) in the room frame.
f = [0; 0; m * g] + GetR_ZYX(hy, hp, hr) * [0; 0; -u4];

% Find applied torque (from rotors) in the body frame.
tau = [u1; u2; u3];

% Find equations of motion in symbolic form. This is a 12x1 matrix h, each
% entry of which is a function of symbolic states and inputs. In our notes
% from class, we would write:
%
%   xdot = h(x, u)
%
h = [v;
     GetN_ZYX(hy, hp, hr) * w;
     (1 / m) * f;
     inv(J) * (tau - wedge(w) * J * w)];

% Find equations of motion in numeric form. This is a MATLAB function h_num
% that can called like any other function as h_num(x, u), with real-valued
% arguments x (a 12x1 matrix of numbers) and u (a 4x1 matrix of numbers).
h_num = matlabFunction(h, 'vars', {x, u});

end


% Returns the rotation matrix R that corresponds to yaw, pitch, and roll
% angles of hy, hp, and hr (assuming a ZYX Euler Angle sequence).
function R = GetR_ZYX(hy, hp, hr)
R = Rz(hy) * Ry(hp) * Rx(hr);
end

% Given yaw, pitch, and roll angles of hy, hp, and hr (assuming a ZYX Euler
% Angle sequence), returns the matrix N for which
%
%   [hydot; hpdot; hrdot] = N * w
%
% where w is the angular velocity and where hydot, hpdot, and hrdot are the
% angular rates.
function N = GetN_ZYX(hy, hp, hr)
R_1inB = Rx(hr);
R_BinA = Ry(hp);
N = inv([(R_BinA * R_1inB)'*[0; 0; 1] (R_1inB)'*[0; 1; 0] [1; 0; 0]]);
end

% Returns the matrix A for which the matrix multiplication A * b implements
% the cross product of two vectors a and b (both written in coordinates).
function A = wedge(a)
    A = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end

% Returns the rotation matrix corresponding to a rotation by an angle h
% about the x axis.
function R = Rx(h)
c = cos(h);
s = sin(h);
R = [ 1  0  0;
      0  c -s;
      0  s  c];
end

% Returns the rotation matrix corresponding to a rotation by an angle h
% about the y axis.
function R = Ry(h)
c = cos(h);
s = sin(h);
R = [ c  0  s;
      0  1  0;
     -s  0  c];
end

% Returns the rotation matrix corresponding to a rotation by an angle h
% about the z axis.
function R = Rz(h)
c = cos(h);
s = sin(h);
R = [ c -s  0;
      s  c  0;
      0  0  1];
end