clear all
close all

% general graphics, this will apply to any figure you open
% (groot is the default figure object).
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 0.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 12, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1.5, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 8, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);
 
% set the tickdirs to go out - need this specific order
%set(groot, 'DefaultAxesTickDir', 'out');
%set(groot, 'DefaultAxesTickDirMode', 'manual');