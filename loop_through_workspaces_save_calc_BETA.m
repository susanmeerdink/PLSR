%% INFO
%Susan Meerdink
%Last Updated: 11/24/2015

%% Loop through spectrums
for spec = 1:6
    
    saveDirectory = 'C:\Users\grad\Dropbox\Code\Matlab\PLSR\Output_Results\'; %Set this Directory to save files into 
    load Hyspiri_8bands_Workspace_2015_11_24.mat %Load workspace so that you have wavelength variables

    %spectra = horzcat(AVIRISspectra,TIRhsypirispectra); %Pulls out and combines into HyspIRI data
    if spec == 1 %Full Spectrum
        
        %Setting Variables
        spectrumName = 'Full';
        beta_ResultFile = strcat(saveDirectory,'beta_results_Full.csv'); %File that will hold regression coefficients from PLSR
        wavelengthIn = vertcat(ASDWavelengths, NicoletWavelengths);
        
    elseif spec == 2 %VSWIR Spectrum

        %Setting Variables
        spectrumName = 'VSWIR';
        beta_ResultFile = strcat(saveDirectory,'beta_results_ASD.csv'); %File that will hold regression coefficients from PLSR
        wavelengthIn = ASDWavelengths ;
        
    elseif spec == 3 %TIR Spectrum
       
        %Setting Variables
        spectrumName = 'TIR';
        beta_ResultFile = strcat(saveDirectory,'beta_results_Nicolet.csv'); %File that will hold regression coefficients from PLSR
        wavelengthIn = NicoletWavelengths;
        
    elseif spec == 4 %HyspIRI Spectrum
       %Setting Variables
        spectrumName = 'HyspIRI';
        beta_ResultFile = strcat(saveDirectory,'beta_results_HyspIRI.csv'); %File that will hold regression coefficients from PLSR
        wavelengthIn = HyspiriWavelengths;

    elseif spec == 5 %AVIRIS Spectrum
           
        %Setting Variables
        spectrumName = 'AVIRIS';
        beta_ResultFile = strcat(saveDirectory,'beta_results_AVIRIS.csv'); %File that will hold regression coefficients from PLSR
        wavelengthIn = AvirisWavelengths;
        
    else %HyTES Spectrum
        
        %Setting Variables
        spectrumName = 'HyTES';
        beta_ResultFile = strcat(saveDirectory,'beta_results_HyTES.csv'); %File that will hold regression coefficients from PLSR
        wavelengthIn = HyTESWavelengths;
                
    end %ends if, elseif, else choosing spectrum
    
    %Set up Files to Store Results
    if exist(beta_ResultFile,'file') > 0 % IF the file already exists, just append to the end 
        %File to hold PLSR results
        fidBetaFile=fopen(beta_ResultFile,'a+');
    else %If the file does not exist create it and write to it
        %File to hold PLSR results
        fidBetaFile=fopen(beta_ResultFile,'w');
        fprintf(fidBetaFile, 'Trait,Spectra,FuncType,Reg/Standard,intercept,'); %Add header to the new file
        fprintf(fidBetaFile, '%f ,',wavelengthIn); %Prints wavelengths of regression coefficients
        fprintf(fidBetaFile, '\n');%Prints a new line character
    end
    
    %% Loop through Traits
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
            elseif functype == 4
                %Spring
                funcName = 'Spring';
            elseif functype == 5
                %Summer
                funcName = 'Summer';
            else
                %Fall
                funcName = 'Fall';
            end %end of if, elseif, and else choosing functional types
            
            %% Open Workspace
            filename = strcat(saveDirectory, 'workspace_', traitName, '_',spectrumName,'_',funcName,'.mat'); %set workspace filename
            load(filename);
            
            %% Calculate Standard BETA 
            standardBETA = [];
            %length(meanBETA)
            
            %Loop through wavelengths
            for i = 2: (length(meanBETA))-1
                standardBETA(i) = (meanBETA(i)/((max(calspectra(:,i))*100)- (min(calspectra(:,i))*100)))*100;
            end
                                    
            %% Output mean and std BETA results to a csv
            %Prints to file that is specified at top of this code
            %Prints a line of text/numbers to file that go in this order:
            
            %Save mean BETA (regression coefficients) to file
            fprintf(fidBetaFile,[traitName ',' spectrumName ',' funcName ', mean,' ]); %Prints the trait name, which spectrum was used, results for what functional type, code (mean,std,or standard), and BETA
            fprintf(fidBetaFile, '%f ,',meanBETA); %Prints wavelengths of regression coefficients
            fprintf(fidBetaFile, '\n');%Prints a new line character
            
            %Save standard deviation of BETA (regression coefficients) to file
            fprintf(fidBetaFile,[traitName ',' spectrumName ',' funcName ', std,' ]); %Prints the trait name, which spectrum was used, results for what functional type, code (mean,std,or standard), and BETA
            fprintf(fidBetaFile, '%f ,',stdBETA); %Prints wavelengths of regression coefficients
            fprintf(fidBetaFile, '\n');%Prints a new line character
            
            %Save standard deviation of BETA (regression coefficients) to file
            fprintf(fidBetaFile,[traitName ',' spectrumName ',' funcName ', standard,']); %Prints the trait name, which spectrum was used, results for what functional type, code (mean,std,or standard), and BETA
            fprintf(fidBetaFile, '%f ,',standardBETA); %Prints wavelengths of regression coefficients
            fprintf(fidBetaFile, '\n');%Prints a new line character
            
            %% Close workspace
            %clearvars -except trait spec functype
        end %end of for loop for functional type
    end %end of for loop for trait
end %end of for loop for spectrum
fclose('all');
