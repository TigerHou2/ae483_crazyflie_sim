function out = constrain(in,limit)
%CONSTRAIN constrains the input to the bounds of +/- limit
%   the user is responsible for providing the correct data shape/format
out = max(min(in,limit),-limit);
end

