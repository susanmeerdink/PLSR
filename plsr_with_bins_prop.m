%Susan Meerdink
%This function calculates the PLSR equation
%This function requires the calibration dataset for spectra and traits
%along with the number of factors generated using determinefactors.m
function [results,stdBETA,meanBETA,valRegLine,valCoeff,valRsq,valRMSE,CIpred,CIsample, calRsq,calRMSE] = plsr_with_bins_prop(calSpectra,valSpectra,calChem,valChem,ID)
%CImodel,
disp('Running PLSR Model')

BETA = zeros(1000,(size(calSpectra,2)+1)); %create variable that will hold BETA coefficients derived from PLSR models
calRsq = zeros(1000,1); %create variable that will hold r-squared values from calibration dataset 
calRMSE = zeros(1000,1); %create variable that will hold precent root mean square errors values from calibration dataset 
valAllCoeff = zeros(1000,length(valChem));%create variable that will hold the trait values derived using the BETA coefficients and calibration spectra
rng('default');

for iteration = 1:1000 %Create 1000 PLSR models 
    indexList = 1:1:(size(calChem,1));
    valIndex = [];
    %% Randomly select from each bin for a total of 30% of data
    if max(calChem)< 1
        step = .1;
    elseif max(calChem) > 100
        if size(calChem,1) < 100
            step = (size(calChem,1)/3);
        else
            step = (size(calChem,1)/10);
        end
    else
        step = (size(calChem,1)/100);
    end
    binranges = min(calChem):step:max(calChem);
    [bincounts,binindex] = histc(calChem,binranges);
    
    for num = 1:size(bincounts,1)
        if bincounts(num) < 0
            continue
        else
            inIndex = indexList(find(binindex == num))';
            random = randperm((bincounts(num)),round((.3*bincounts(num)))); %Get random numbers in this bin range
            valIndex = vertcat(valIndex, inIndex(random));%Add them to the validation list
        end
    end
        
    %% Look for duplicates and replace if necessary
    checkDup = unique(valIndex,'sorted');
    count = 1;
    while size(checkDup,1) < floor(.3*size(calChem,1))
        random = randperm((size(calChem,1)),(size(valIndex,1) - size(checkDup,1)));
        valIndex = checkDup; %Remove duplicate
        valIndex = vertcat(valIndex, inIndex(random));%Add them to the validation list
        checkDup = unique(valIndex);
        if count == 100 %If it hasn't been able to find a number that isn't a duplicate by 100 iterations quit and move on
            break
        end
        count = count +1;
    end
    
    %% Set final validation and calibration datasets
    %Validation (30% of data)
    tempValTrait = calChem(valIndex);
    tempValSpectra = calSpectra([valIndex],:);
    %valSeason = season(valIndex);
    %valSpeciesName = speciesName(val(:,2));
    
    %Calibration (70% of data)
    tempTrait = calChem;
    tempTrait(valIndex) = []; %Calibration with validation removed
    tempSpectra = calSpectra;
    tempSpectra([valIndex],:) = [];
            
    %Randomly select 70% of calibration dataset
%     random = randperm(size(calSpectra,1),round((.7*size(calSpectra,1)))); %Select random indices
%     tempSpectra = calSpectra([random],:); % Set the temporaray calibration dataset to only have spectra pulled with random set of indices 
%     tempTrait = calChem(random); % Set the temporaray calibration dataset to only have traits pulled with random set of indices
%     tempValTrait = calChem; %Set the temporary validation dataset equal to the total calibration dataset
%     tempValTrait([random],:) = []; %remove the random indices from trait dataset
%     tempValSpectra = calSpectra; %Set the temporary validation dataset equal to the total calibration dataset
%     tempValSpectra([random],:) = [];%remove the random indices from trait dataset
%     
    %Run PLSR and save BETA coefficients
    [~,~,~,~,BETA1,~,~,~] = plsregress(tempSpectra,tempTrait,ID); 
    
    %Calculate Results
    calCoeff1 = [ones(size(tempValSpectra,1),1),tempValSpectra]*BETA1; %Determine what the trait values will be based on the BETA derived from PLSR
    calMDL = LinearModel.fit(tempValTrait,calCoeff1,'linear'); %Calculate a regression equation between the predicted and observed
    calRMSE1 = (calMDL.RMSE/(max(tempValTrait)-min(tempValTrait)))*100; %Calculate the percent RMSE from the regression equation
    
    %Set Variables for later
    BETA([iteration],:) = BETA1'; %Store the BETA coefficients for future use
    calRsq(iteration) = calMDL.Rsquared.Ordinary;  %Store the regression model r-squared
    calRMSE(iteration) = calRMSE1; %Store the regression model RMSE
end

%% Math
stdBETA = std(BETA)'; %Calculate the Standard deviation of the BETA coefficients and transpose them
meanBETA = mean(BETA)'; %Calculate the mean of the Beta coefficients and transpose them

valCoeff = [ones(size(valSpectra,1),1),valSpectra]*meanBETA; %Calculate the predicted trait values using the mean BETA coefficients and the set aside validation dataset
valMDL = LinearModel.fit(valCoeff,valChem,'linear','RobustOpts','on'); %Calculate the regression results of predicted versus observed
valRegLine = [valMDL.Coefficients.Estimate(1), valMDL.Coefficients.Estimate(2)]; %Pull out the regression line
%Coefficients.Estimate(1) = intercept
%Coefficients.Estimte(2) = slope
valRMSE = (valMDL.RMSE/(max(valChem)-min(valChem)))*100; %Calculate the percent RMSE
valRsq = valMDL.Rsquared.Ordinary; %Pull out the r-squared of predicted versus observed
%valCI = coefCI(valMDL);
%valCIinter = valCI(1,:);
%valCIx = valCI(2,:);
figure
hold on
boxplot(calRMSE)
plot(1,valRMSE,'o')
hold off
figure
hold on
hist(calRMSE)
plot(valRMSE,1,'ro')
hold off
figure
hold on
boxplot(calRsq)
plot(1,valRsq,'ro')
hold off
figure
hold on
hist(calRsq)
plot(valRsq,1,'ro')
hold off

%% Confidence Intervals for Validation Dataset
%Loop through BETA to get range of values for sample validation
for num = 1:1000
    valAllCoeff([num],:) = [ones(size(valSpectra,1),1),valSpectra]*(BETA(num,:)'); %Calculate the predicted trait for each run
end
%Calculate the confidence Intervals of validation samples
for sampleNum = 1:length(valChem)
    %hist(valAllCoeff(:,sampleNum))
    SEM = std(valAllCoeff(:,sampleNum))/sqrt(length(valAllCoeff(:,sampleNum)));        % Standard Error
    ts = tinv([0.025  0.975],length(valAllCoeff(:,sampleNum))-1);      % T-Score
    CIsample(:,sampleNum) = mean(valAllCoeff(:,sampleNum)) + ts*SEM; % Confidence Intervals
end 

%% Confidence Intervals for Prediction
ci = coefCI(valMDL); %Calculate the 95% Confidence Interval 
CIpred = [ci(1,1) ci(1,2)]; %Add this confidence interval for intercepts to a new variable.

%% Output Table
results = zeros(1,10);
%Number of Samples, Number used for Validation, # of Factors, Rsq of Calibration, RMSE of calibration, std of calibration, Rsq of validation, RMSE of validation\n'
results = [(size(calChem,1)+size(valChem,1)), size(valChem,1), ID, mean(calRsq), std(calRsq), mean(calRMSE), std(calRMSE), valRsq, valRMSE];

end