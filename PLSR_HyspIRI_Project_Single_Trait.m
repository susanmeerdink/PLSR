%% PLS regression Code
%Susan Meerdink

%This PLSR code is being adjusted based on RSE reviewer comments
% 1. It splits total number of samples into calibration and validation
% datasets. splitValCalMore.m

% 2. Then determines the number of components/factors by using leave one out  
% cross validation. determinefactors.m

% 3. The number of factors are then used to calculate and formulate the
% PLSR regression. plsr.m

% 4. The product of the PLSR is then input into plsfigure.m to display a
% resulting figure. 
  
%% Select Trait to Analyze
%Trait Numbers go from 1 - 5, select the trait you want to analyze
% 1: Nitrogen (%)
% 2: Cellulose (%)
% 3: Lignin (%)
% 4: Water Content (%)
% 5: Leaf Mass per Area (LMA) (g/m^2)
trait = 2;

    %% Load Data
load Hyspiri_8bands_Workspace_2015_11_24.mat
% includes spectra and chem values
traitResults = zeros(90,12);
allCalRMSE = zeros(1000,5);
allCalRsq = zeros(1000,5);

%% Set up Files to Store Results
directory = 'E:\Meerdink\Dropbox\Code\Matlab\PLSR\Output_Results\'; %Set Directory
plsrResultFile = '2016_01_22_plsr_resultsCellulose.csv'; %File that will hold results/accuries from PLSR
betaResultFile = '2016_01_22_beta_resultsCellulose.csv'; %File that will hold the BETA coefficients derived from the PLSR

%PLSR FILE
if exist(plsrResultFile,'file') > 0 % IF the file already exists, just append to the end 
    %File to hold PLSR results
    fidPLSR=fopen(strcat(directory, plsrResultFile),'a+');
else %If the file does not exist create it and write to it
    %File to hold PLSR results
    fidPLSR=fopen(strcat(directory, plsrResultFile),'w');
    fprintf(fidPLSR, 'Trait,Spectra,FuncType, Number of Samples, Number of Validation, # of Factors, Mean Rsq of Calibration, SD Rsq of Calibration, Mean RMSE of calibration, SD RMSE of calibration, Rsq of validation, RMSE of validation\n'); %Add header to the new file
end

%BETA FILE
if exist(plsrResultFile,'file') > 0 % IF the file already exists, just append to the end 
    %File to hold Beta coefficients
    fidBeta=fopen(strcat(directory, betaResultFile),'a+');
else %If the file does not exist create it and write to it
    %File to hold Beta coefficients
    fidBeta=fopen(strcat(directory, betaResultFile),'w');
end

%% Setting Variables
if trait == 1
    traitName = 'Nitrogen';
elseif trait == 2
    traitName = 'Cellulose';
elseif trait == 3
    traitName = 'Lignin';
elseif trait == 4
    traitName = 'Water Content';
else
    traitName = 'LMA';
end
disp(['Running PLSR for '  traitName])
allChem = chemistry(:,trait); %pulls out the whole column of data

%% Determine Validation/Calibration Samples

[valIndexAll, valIndexBroad,valIndexNeedle,valIndexSpring,valIndexSummer,valIndexFall] = splitValCalBefore(allChem,sampleInfoNum.BroadFunction,sampleInfoNum.Season);
    
%% Loop through spectrums
for spec = 1:6
    %spectra = horzcat(AVIRISspectra,TIRhsypirispectra); %Pulls out and combines into HyspIRI data
    if spec == 1
        %Full Spectrum
        spectraAll = horzcat(ASDspectra, Nicoletspectra); %Pulls out and combines into HyspIRI data
        spectrumName = 'Full';
        wavelengths = waveFull;
        %ID = 18;
    elseif spec == 2
        %VSWIR Spectrum
        spectraAll = ASDspectra; %Pulls out and combines into HyspIRI data
        spectrumName = 'VSWIR';
        wavelengths = ASDWavelengths;
        %ID =14;
    elseif spec == 3
        %TIR Spectrum
        spectraAll = Nicoletspectra;
        spectrumName = 'TIR';
        wavelengths = NicoletWavelengths;
    elseif spec == 4
        %HyspIRI Spectrum
        spectraAll = horzcat(AVIRISspectra,HyspiriTIRspectra); %Pulls out and combines into HyspIRI data
        spectrumName = 'HyspIRI';
        wavelengths = HyspiriWavelengths;
        %ID = 15;
    elseif spec == 5
        %AVIRIS Spectrum
        spectraAll = AVIRISspectra;
        spectrumName = 'AVIRIS';
        wavelengths = AvirisWavelengths;
        %ID = 15;
    else 
        %HyTES Spectrum
        spectraAll = HyTESspectra; %Pulls out and combines into HyspIRI data
        spectrumName = 'HyTES';
        wavelengths = HyTESWavelengths;
    end %ends if, elseif, else choosing spectrum
        
    %% Loop through Functional Types
    for functype = 1:6
        allChem = chemistry(:,trait); %pulls out the whole column of data
        if functype == 1
            %Broadleaf
            funcName = 'Broadleaf';
            allChem = chemistry(:,trait); %pulls out the whole column of data
            allChem = allChem(find(sampleInfoNum.BroadFunction == 1 | sampleInfoNum.BroadFunction == 2));
            spectra = spectraAll((find(sampleInfoNum.BroadFunction == 1 | sampleInfoNum.BroadFunction == 2)),:);
            species = sampleInfoNum.Species(find(sampleInfoNum.BroadFunction == 1 | sampleInfoNum.BroadFunction == 2));
            season = sampleInfoNum.Season(find(sampleInfoNum.BroadFunction == 1 | sampleInfoNum.BroadFunction == 2));
            disp(['Running PLSR for broadleaf '  traitName ' ' spectrumName])
            valIndex = valIndexBroad;
            batchName = 'Broadleaf';
        elseif functype == 2
            %Needleaf
            funcName = 'Needleleaf';
            allChem = chemistry(:,trait); %pulls out the whole column of data
            allChem = allChem(find(sampleInfoNum.BroadFunction == 3));
            spectra = spectraAll((find(sampleInfoNum.BroadFunction == 3)),:);
            species = sampleInfoNum.Species(find(sampleInfoNum.BroadFunction == 3));
            season = sampleInfoNum.Season(find(sampleInfoNum.BroadFunction == 3));
            disp(['Running PLSR for needleaf '  traitName ' ' spectrumName])
            valIndex = valIndexNeedle;
            batchName = 'Needleleaf';
        elseif functype == 3
            %All Samples
            funcName = 'All Samples';
            allChem = chemistry(:,trait); %pulls out the whole column of data
            spectra = spectraAll;
            disp(['Running PLSR for general '  traitName ' ' spectrumName])
            batchName = 'General';
            species = sampleInfoNum.Species;
            season = sampleInfoNum.Season;
            valIndex = valIndexAll;
        elseif functype == 4
            %Spring
            funcName = 'Spring';
            allChem = chemistry(:,trait); %pulls out the whole column of data
            allChem = allChem(find(sampleInfoNum.Season == 1));
            spectra = spectraAll((find(sampleInfoNum.Season == 1)),:);
            species = sampleInfoNum.Species(find(sampleInfoNum.Season == 1));
            season = sampleInfoNum.Season(find(sampleInfoNum.Season == 1));
            disp(['Running PLSR for spring '  traitName ' ' spectrumName])
            valIndex = valIndexSpring;
            batchName = 'Spring';
        elseif functype == 5
            %Summer
            funcName = 'Summer';
            allChem = chemistry(:,trait); %pulls out the whole column of data
            allChem = allChem(find(sampleInfoNum.Season == 2));
            spectra = spectraAll((find(sampleInfoNum.Season == 2)),:);
            species = sampleInfoNum.Species(find(sampleInfoNum.Season == 2));
            season = sampleInfoNum.Season(find(sampleInfoNum.Season == 2));
            disp(['Running PLSR for summer '  traitName ' ' spectrumName])
            valIndex = valIndexSummer;
            batchName = 'Summer';
        else
            %Fall
            funcName = 'Fall';
            allChem = chemistry(:,trait); %pulls out the whole column of data
            allChem = allChem(find(sampleInfoNum.Season == 3));
            spectra = spectraAll((find(sampleInfoNum.Season == 3)),:);
            species = sampleInfoNum.Species(find(sampleInfoNum.Season == 3));
            season = sampleInfoNum.Season(find(sampleInfoNum.Season == 3));
            disp(['Running PLSR for Fall '  traitName ' ' spectrumName])
            valIndex = valIndexFall;
            batchName = 'Fall';
        end %end of if, elseif, and else choosing functional types

        %% Finding and removing NaN Values
        nanListFull = find(isnan(allChem)); %Finds all NaN values in trait dataset being analyzed
        allChem(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
        spectra([nanListFull],:) = []; %Removes all NaN values from spectra dataset being analyzed
        species(nanListFull) = [];%Removes all NaN values from species list
        season(nanListFull) = [];%Removes all NaN values from the season list

        %% Set calibration and validation
        %Call function splitValCalMore.m
        %[valChem,valspectra,valSeason, calChem,calspectra] = splitValCalMore(allChem,spectra,species,season);
        %Validation (20% of data)
        valTrait = allChem(valIndex);
        valSpectra = spectra([valIndex],:);
        valSeason = season(valIndex);
        
        %Calibration (80% of data)
        calTrait = allChem;
        calTrait(valIndex) = []; %Calibration with validation removed
        calSpectra = spectra;
        calSpectra([valIndex],:) = [];
        
        disp(['Total number of samples: '  num2str(size(allChem,1))])
        disp(['Number of samples set aside for validation: '  num2str(size(valTrait,1))])
        disp(['Number of samples set aside for calibration: '  num2str(size(calTrait,1))])

        %% Determine Components for PLSR
        %Calls function determinefactors.m
        [PRESSRMSEY, ID, Min, meanPCTVAR] = determinefactors(spectra,allChem,funcName,spectrumName);
        % ID: the number of factors to be used in PLSR
        % Min: The minimum PRESS statistic value
        % meanPCTVAR: the percent variation explained by each factor
        % PRESSRMSEY: the press statistic for each factor

        %% Determine PLSR Model
        %calls function plsr.m
        [results,stdBETA,meanBETA,valRegLine,valCoeff,valRsq,valRMSE,CIpred,CIsample,calRsq,calRMSE] = plsr_with_bins_prop(calSpectra,valSpectra,calTrait,valTrait,ID);
        %CImodel,

        %Save variables
        disp(['Done with PLSR for ' funcName ' ' spectrumName])
        traitResults(trait,:) = horzcat(spec, trait,functype,results); %Save results from Plsr into a table for displaying and printing
        allCalRsq(:,trait) = calRsq;
        allCalRMSE(:,trait) = calRMSE;

        %% PLSR Predicted versus Observed Plots
        %This calls function plsfigure.m
        plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)%CImodel,

        %% Output results to a csv
        %Prints to file that is specified at top of this code
        %Prints a line of text/numbers to file that go in this order:
        %'Trait,Spectra,FuncType, Trait #, Number of Samples, Number used for Validation, # of Factors, Rsq of Calibration, RMSE of calibration, std of calibration, Rsq of validation, RMSE of validation\n'
        fprintf(fidPLSR,[traitName ',' spectrumName ',' funcName ',' ]); %Prints the trait name, which spectrum was used, and results for what functional type 
        fprintf(fidPLSR,'%d , %d , %d , %f , %f , %f , %f , %f , %f',results); %Prints plsr results
        fprintf(fidPLSR,'\n'); %Prints a new line character

        %% Output results of BETA coefficients to csv
        %Prints to file that is specified at top of this code
        fprintf(fidBeta,[traitName ',' spectrumName ',' funcName '\n']);%Prints the trait name, which spectrum was used, and results for what functional type 
        fprintf(fidBeta,'Wavelengths,intercept,');%Prints row headers
        fprintf(fidBeta,'%f ,',wavelengths);%Prints the wavelengths of the spectrum used in PLSR model
        fprintf(fidBeta,'\n');%Prints a new line character
        fprintf(fidBeta,'Mean BETA,'); %Prints row headers
        fprintf(fidBeta,'%f ,',meanBETA); %Prints the mean BETA coefficients used in the PLSR model
        fprintf(fidBeta,'\n');%Prints a new line character
        fprintf(fidBeta,'SD BETA,');%Prints row header
        fprintf(fidBeta,'%f ,',stdBETA);%Prints the standard deviation of the BETA coefficients used in the PLSR model
        fprintf(fidBeta,'\n');%Prints a new line character
        fprintf(fidBeta,'\n');%Prints another new line character to add space between models
        
        %% Save workspace
        filename = strcat(directory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_2016_01_22'); %set workspace filename
        save(filename); %Save the workspace
        close all %Close all figures that have been opened.

    end %End of loop going through functional types     
end  %End of loop goieptng through spectrums
fclose('all');