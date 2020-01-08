 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME: magnifyOnFigure_examples
% 
% AUTHOR: David Fernandez Prim (david.fernandez.prim@gmail.com)
%
% PURPOSE: Shows the funcionality of 'magnifyOnFigure'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default interactive mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp( sprintf('This is the default interactive operation mode of ''magnifyOnFigure''') ) 
fig = figure;
hold on;
plot(rand(100,1), 'b'); plot(rand(300, 1), 'r'); 
grid on;
hold off;
magnifyOnFigure;
disp('Press a key...')
pause;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure handle passed as an input argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp( sprintf('The figure handle is here passed as an input argument.') ) 
fig = figure;
hold on;
plot(rand(100,1), 'b'); plot(rand(300, 1), 'r'); 
grid on;
hold off;
magnifyOnFigure(fig);
disp('Press a key...')
pause;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Properties (in interactive mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp( sprintf('Playing arround with the properties in interactive mode...') ) 
figHandler = figure;
hold on;
plot(rand(100,1), 'b'); plot(rand(300, 1), 'r'); 
grid on;
hold off; 
ylim([0 2]);
magnifyOnFigure(...
        figHandler,...
        'units', 'pixels',...
        'magnifierShape', 'ellipse',...
        'initialPositionSecondaryAxes', [326.933 259.189 164.941 102.65],...
        'initialPositionMagnifier',     [174.769 49.368 14.1164 174.627],...    
        'mode', 'interactive',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 2,...
        'edgeColor', 'black',...
        'secondaryAxesFaceColor', [0.91 0.91 0.91]... 
            ); 
disp('Press a key...')
pause;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Properties (in manual mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp( sprintf('Or in manual mode.') ) 
figHandler = figure;
hold on;
plot(rand(100,1), 'b'); plot(rand(300, 1), 'r'); 
grid on;
hold off; 
ylim([0 2]);
magnifyOnFigure(...
        figHandler,...
        'units', 'pixels',...
        'initialPositionSecondaryAxes', [326.933 259.189 164.941 102.65],...
        'initialPositionMagnifier',     [174.769 49.368 14.1164 174.627],...    
        'mode', 'manual',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 2,...
        'edgeColor', 'black',...
        'secondaryAxesFaceColor', [0.91 0.91 0.91]... 
            ); 
disp('Press a key...')
pause;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Working on images also
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp( sprintf('How the tool works on images.') ) 
h = figure; 
load clown; 
image(X); 
colormap(map) ; 
axis image
magnifyOnFigure(h, 'displayLinkStyle', 'straight',...
                    'EdgeColor', 'white',...
                    'magnifierShape', 'rectangle',...
                    'frozenZoomAspectratio', 'on',...
                    'edgeWidth', 2);
                
disp('Press a key...')
pause;
close all
                
             

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Working on contour plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp( sprintf('How the tool works on contour plots.') ) 
scrsz = get(0, 'ScreenSize');
h = figure('Position', [0.01*scrsz(3), 0.25*scrsz(4), 0.65*scrsz(3), 0.60*scrsz(4)]);
load clown;
hc1 = contour(X, 'LineWidth', 2);
axis image

magnifyOnFigure;

disp('Press a key...')
pause;
close all


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Properties (in interactive mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp( sprintf('Using multiple magnifiers on the same axis...') ) 
figHandler = figure;
hold on;
plot(rand(100,1), 'b'); plot(rand(300, 1), 'r'); 
grid on;
hold off; 
ylim([0 2]);

magnifyOnFigure(...
        figHandler,...
        'units', 'pixels',...
        'magnifierShape', 'rectangle',...
        'initialPositionSecondaryAxes', [105.6 278.81 130.2 102.69],...
        'initialPositionMagnifier',     [106.56 47.2 18.6823 186.519],...    
        'mode', 'interactive',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 2,...
        'edgeColor', 'black',...
        'secondaryAxesFaceColor', [0.91 0.91 0.91]... 
            );   

magnifyOnFigure(...
        figHandler,...
        'units', 'pixels',...
        'magnifierShape', 'rectangle',...
        'initialPositionSecondaryAxes', [365.6 275.81 130.2 102.69],...
        'initialPositionMagnifier',     [211.459 47.2 18.6823 186.519],...    
        'mode', 'interactive',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 2,...
        'edgeColor', 'black',...
        'secondaryAxesFaceColor', [0.91 0.91 0.91]... 
            );   
        
magnifyOnFigure(...
        figHandler,...
        'units', 'pixels',...
        'magnifierShape', 'rectangle',...
        'initialPositionSecondaryAxes', [364.6 78.81 130.2 102.69],...
        'initialPositionMagnifier',     [270.827 47.2 18.6823 186.519],...    
        'mode', 'interactive',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 2,...
        'edgeColor', 'black',...
        'secondaryAxesFaceColor', [0.91 0.91 0.91]... 
            ); 
        
     
        
