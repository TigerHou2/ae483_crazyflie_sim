function out = runPid(P,I,D,pid)
%RUNPID calculates the output of a PID controller
%   P - Nx1 vector of proportional errors
%   I - Nx1 vector of integral errors
%   D - Nx1 vector of derivative errors
% pid - Nx3 matrix of gains

out = P .* pid(:,1) ...
    + I .* pid(:,2) ...
    + D .* pid(:,3) ;

end

