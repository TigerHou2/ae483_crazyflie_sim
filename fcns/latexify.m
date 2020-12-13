function latexify(w,h,fontSize)
%LATEXIFY Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    w = 18;
    h = 18;
end
if ~exist('fontSize','var')
    fontSize = 16;
end

set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)

set( gca, 'Color', [1 1 1] )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [w h])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0.2 1.2 w h])
set(gcf, 'PaperPosition', [0.2 1.2 w h])

end

