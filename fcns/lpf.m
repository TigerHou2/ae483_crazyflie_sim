function [LPF,out] = lpf(LPF,in,sample_freq,cutoff_freq)
%LPF Summary of this function goes here
%   Detailed explanation goes here

%% initialize the low pass filter if it does not exist
if isempty(LPF)
    LPF = struct('b0',0,'b1',0,'b2',0,'a1',0,'a2',0,'d1',0,'d2',0);
    fr = sample_freq/cutoff_freq;
    ohm = tan(pi/fr);
    c = 1 + 2*cos(pi/4)*ohm + ohm^2;
    LPF.b0 = ohm^2 / c;
    LPF.b1 = 2 * LPF.b0;
    LPF.b2 = LPF.b0;
    LPF.a1 = 2 * (ohm^2-1) / c;
    LPF.a2 = (1 - 2*cos(pi/4)*ohm + ohm^2) / c;
    LPF.d1 = 0;
    LPF.d2 = 0;
    out = 0;
%% run the low pass filter
else
    d0 = in - LPF.d1*LPF.a1 - LPF.d2*LPF.a2;
    if isnan(d0) || isinf(d0)
        d0 = 0;
    end
    out = d0*LPF.b0 + LPF.d1*LPF.b1 + LPF.d2*LPF.b2;
    LPF.d2 = LPF.d1;
    LPF.d1 = d0;
end

%% reset the low pass filter
% if reset
%     dval = sample / (LPF.b0+LPF.b1+LPF.b2);
%     LPF.d1 = dval;
%     LPF.d2 = dval;
%     out = 0;
% end

end

