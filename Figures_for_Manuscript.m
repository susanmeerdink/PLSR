%This code will go through and create coefficient figures for the
%manuscript
%Susan Meerdink
%3/3/2016
%% INPUTS
saveDirectory = 'H:\users\meerdink\Dropbox\Code\Matlab\PLSR\Output_Results\'; %Set Directory
workspaceName = '2016_02_25_wSeasons'; %The date added on to the workspace name
topWave = [];

%% Cellulose Full
%Variables
traitName = 'Cellulose';
spectrumName = 'Full';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 3;
yStep = 1;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,           

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));

%% Cellulose HyspIRI
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'Cellulose';
spectrumName = 'HyspIRI';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 60;
yStep = 20;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,  

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));
%% Lignin Full
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'Lignin';
spectrumName = 'Full';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 1.2;
yStep = .4;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel, 

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));

%% Lignin AVIRIS
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'Lignin';
spectrumName = 'AVIRIS';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 60;
yStep = 20;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel, 

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));
%% LMA Full
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'LMA';
spectrumName = 'Full';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 30;
yStep = 10;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));
%% LMA HyspIRI
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'LMA';
spectrumName = 'HyspIRI';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 450;
yStep = 150;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));
%% Nitrogen TIR
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'Nitrogen';
spectrumName = 'TIR';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 0.3;
yStep = 0.1;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));
%% Nitrogen HyspIRI
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'Nitrogen';
spectrumName = 'HyspIRI';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 4;
yStep = 1;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));
%% Water Content Full
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'Water Content';
spectrumName = 'Full';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 0.012;
yStep = 0.004;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));
%% Water Content HyspIRI
% Clean up workspace
close all %Close the figures
clearvars -except trait traitName spec spectrumName functype workspaceName saveDirectory topWave

%Variables
traitName = 'Water Content';
spectrumName = 'HyspIRI';
funcName = 'All Samples';

% Open Workspace
filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
disp(filename)
load(filename);

% Display figure
y = 0.5;
yStep = 0.1;
betafigure(wavelengths,meanStandCoeff,stdStandCoeff,traitName,spec,y,yStep)

% PLSR Predicted versus Observed Plots
plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,

% find largest coefficients
[B,I] = findpeaks(meanStandCoeff,'SortStr','descend');
topWave = horzcat(topWave,wavelengths(I(1:10)));