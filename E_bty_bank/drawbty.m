clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.
fontSize = 16;
axis on;
letsGo = 1;
while letsGo
    % set range limits
    rmin = input('Enter the lowest range (IN KM!! can be negative) to draw in:\n');
    rmax = input('Enter the highest range (IN KM!!) to draw in:\n');
    maxDepth = input('Enter the maximum depth\n');
    axis([rmin rmax -maxDepth 0])
    title('Bty', 'FontSize', fontSize);
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
    uiwait(msgbox(message));
    % User signs name here.
    hFH = imfreehand();
    % Get the xy coordinates of where they drew.
    xy = hFH.getPosition;
    % get rid of imfreehand remnant.
    delete(hFH);
    % Overlay what they drew onto the image.
    hold on; % Keep image, and direction of y axis.
    xCoordinates = xy(:, 1);
    yCoordinates = -xy(:, 2);
    plot(xCoordinates, -yCoordinates, 'ro-', 'LineWidth', 2);
    caption = sprintf('Points may not lie on adjacent pixels, depends on your speed of drawing!', 'FontSize', fontSize);
    title(caption, 'FontSize', fontSize);
    % Ask user if they want to save the bty file.
    promptMessage = sprintf('Do you want to save this bathymetry?');
    titleBarCaption = 'Continue?';
    button = questdlg(promptMessage, titleBarCaption, 'OK', 'Cancel', 'OK');
    if strcmpi(button, 'OK')
        letsGo = 0;
    end
    if strcmpi(button, 'Cancel')
        clf;
        close all;
        imtool close all
    end
    
end

% Save the data
filename = input('What do you want to name the file? (Include extension).\n','s');
fid = fopen(filename,'w');
fprintf(fid,'''L''\n'); %L for linear interp, C for curved
fprintf(fid,'%u\n',length(xCoordinates)); %number of points 
for i = 1:length(xCoordinates)
    fprintf(fid,'%1$-6.3f %2$-6.3f\n',xCoordinates(i),yCoordinates(i)); 
end