%% Set up Files to Store Results
saveDirectory = 'H:\users\meerdink\Dropbox\Code\Matlab\PLSR\Output_Results\'; %Set Directory
plsrResultFile = '2016_02_25_plsr_results_complied.csv'; %File that will hold results/accuries from PLSR
workspaceName = '2016_02_25_wSeasons'; %The date added on to the workspace name
%PLSR FILE
if exist(plsrResultFile,'file') > 0 % IF the file already exists, just append to the end 
    %File to hold PLSR results
    fidPLSRAll=fopen(strcat(saveDirectory, plsrResultFile),'a+');
else %If the file does not exist create it and write to it
    %File to hold PLSR results
    fidPLSRAll=fopen(strcat(saveDirectory, plsrResultFile),'w');
    fprintf(fidPLSRAll, 'Trait,Spectra,FuncType, Number of Samples, Number of Validation, # of Factors, Mean Rsq of Calibration, SD Rsq of Calibration, Mean RMSE of calibration, SD RMSE of calibration, Rsq of validation, RMSE of validation\n'); %Add header to the new file
end


allTraitResults = []; %Create array to hold results from PLSR models
%% Setting Variables
for trait = 1:5
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


    %% Loop through spectrums
    for spec = 1:6
        %spectra = horzcat(AVIRISspectra,TIRhsypirispectra); %Pulls out and combines into HyspIRI data
        if spec == 1
            %Full Spectrum
            spectrumName = 'Full';
        elseif spec == 2
            %VSWIR Spectrum
            spectrumName = 'VSWIR';
        elseif spec == 3
            %TIR Spectrum
            spectrumName = 'TIR';
        elseif spec == 4
            %HyspIRI Spectrum
            spectrumName = 'HyspIRI';
        elseif spec == 5
            %AVIRIS Spectrum
            spectrumName = 'AVIRIS';
        else 
            %HyTES Spectrum
            spectrumName = 'HyTES';
        end %ends if, elseif, else choosing spectrum

        %% Loop through Functional Types
        for functype = 1:6

            if functype == 1
                %Broadleaf
                funcName = 'Broadleaf';
            elseif functype == 2
                %Needleaf
                funcName = 'Needleleaf';
            elseif functype == 3
                %All Samples
                funcName = 'All Samples';
            elseif  functype == 4
                funcName = 'Spring';
            elseif functype == 5
                funcName = 'Summer';
            else %Fall
                funcName = 'Fall';
            end %end of if, elseif, and else choosing functional types
            
            %% Open Workspace
            filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'_',workspaceName,'.mat'); %set workspace filename
            disp(filename)
            load(filename);
            
            %% SAve results to new array
            allTraitResults = vertcat(allTraitResults, results);
                        
            %% Output results to a csv
            %Prints to file that is specified at top of this code
            %Prints a line of text/numbers to file that go in this order:
            %'Trait,Spectra,FuncType, Trait #, Number of Samples, Number used for Validation, # of Factors, Rsq of Calibration, RMSE of calibration, std of calibration, Rsq of validation, RMSE of validation\n'
            fprintf(fidPLSRAll,[traitName ',' spectrumName ',' funcName ',' ]); %Prints the trait name, which spectrum was used, and results for what functional type 
            fprintf(fidPLSRAll,'%d , %d , %d , %f , %f , %f , %f , %f , %f',results); %Prints plsr results
            fprintf(fidPLSRAll,'\n'); %Prints a new line character
            
            %% Close workspace
            clearvars -except trait traitName spec spectrumName functype fidPLSRAll saveDirectory allTraitResults workspaceName
        end %end of for loop for functional type
    end %end of for loop for spectrum
end %end of for loop for trait
fclose('all');
disp('done');
