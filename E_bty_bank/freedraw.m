clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.
fontSize = 16;
% Read in a standard MATLAB gray scale demo image.
folder = fullfile(matlabroot, '\toolbox\images\imdemos');
baseFileName = 'cameraman.tif';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
	% File doesn't exist -- didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
grayImage = imread(fullFileName);
imshow(grayImage, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
message = sprintf('Sign your name.\nLeft click and hold to begin drawing.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
% User draws a 2d slice of the ocean bathymetry here.
hFH = imfreehand();
% Get the xy coordinates of where they drew.
xy = hFH.getPosition
% get rid of imfreehand remnant.
delete(hFH);
% Overlay what they drew onto the image.
hold on; % Keep image, and direction of y axis.
xCoordinates = xy(:, 1);
yCoordinates = xy(:, 2);
plot(xCoordinates, yCoordinates, 'ro-', 'LineWidth', 2);
caption = sprintf('Original Grayscale Image.\nPoints may not lie on adjacent pixels, depends on your speed of drawing!', 'FontSize', fontSize);
title(caption, 'FontSize', fontSize);
% Ask user if they want to burn the line into the image.
promptMessage = sprintf('Do you want to burn the line into the image?');
titleBarCaption = 'Continue?';
button = questdlg(promptMessage, titleBarCaption, 'OK', 'Cancel', 'OK');
if strcmpi(button, 'Cancel')
	return;
end
cla;
hold off;
for k = 1 : length(xCoordinates)
	row = int32(yCoordinates(k));
	column = int32(xCoordinates(k));
	grayImage(row, column) = 255;
end
imshow(grayImage, []);
axis on;
caption = sprintf('Grayscale Image with Burned In Curve.\nPoints may not lie on adjacent pixels, depends on your speed of drawing!', 'FontSize', fontSize);
title(caption, 'FontSize', fontSize);
