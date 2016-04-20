%Susan Meerdink
%This function calculates the number of factors for PLSR
function [meanPRESSY, ID, Min, PCT, meanPCTVAR] = determinefactors(spectra,trait,funcName,spectrumName,traitName)

disp('Determining number of components')
ncomp = 20; %Setting the maximum number of factors to 30 (in order to save time). If time is not an issue, loop through the total number of Factors possible using code below
% nobs = size(spectra,1); %Determine the size of the dataset
% ncomp = nobs - 2; %Determine the total number of factors that could be used 
c = cvpartition((size(trait,1)),'leaveout'); % Parition the dataset into calibration and validation using leave one out cross validation
PRESSY = zeros(c.NumTestSets,ncomp+1); %create variable that will hold PRESS statistic
PCTVAR = zeros(c.NumTestSets,ncomp); %Create variaible that will hold the percent variation explained at each factor
for i = 1:c.NumTestSets %Loop through the test sets
    [~,~,~,~,~,PCTVAR1,MSE,~] = plsregress(spectra(c.training(i),:),trait(c.training(i)),ncomp);
    PRESSY(i,:) = MSE(2,:); %Pull out PRESS stat for y
    PCTVAR(i,:) = cumsum(100*PCTVAR1(2,:)); %Pull out percent variation explained for y
end
%Note: PLSR setting components or factors is run on full data set

%% Determine the Number of Factors using the minimum PRESS statistic Value
PRESSY(:,1) = []; %Drop column 1 as this is for 0 factors
meanPRESSY = mean(PRESSY); %Calculate the mean PRESS statistic from the above runs
meanPCTVAR = mean(PCTVAR); %Calcualte the mean percent variation explained from the above runs

%Go through the minimum PRESSRMSEY to figure out which factor to choose
localMinPRESS = [];
for k = 1:size(meanPRESSY,1) %Loop through PRESSRMSEY list to check and see if there is a local minimum
    if meanPRESSY(k) < meanPRESSY(k+1)  %If the current index has a slope less then the next 
        localMinPRESS = vertcat(localMinPRESS,(k));
    end 
end


if isempty(localMinPRESS) == 1 %if there were no values that fit the above requirements
    [Minm,IDm] = min(meanPRESSY); % Pulls lowest PRESS statistic to determine number of factors
    PCTm = meanPCTVAR(IDm); %Pulls out the percent variation explained at that factor
else
    %Determine number of factors for PRESSRMSEY
    [~,tempIndex] = max(localMinPRESS(:,2));%Set the minimum number determined by slope
    IDm = localMinPRESS(tempIndex,1);%Set the minimum number of factors determined by slope
    Minm = meanPRESSY(IDm);
    PCTm = meanPCTVAR(IDm);%Pulls out the percent variation explained at that factor
end 

%% Determine Components using Ttest
% %Variables Need for Analysis 
% tresults = zeros(1,ncomp); %Holds results of ttest
% sigFactor = []; %Holds the non significant factors
% 
% 
% %T-Test to determine if there are factors that are NOT significantly
% %different from each other
% for i = 1:(ncomp-1)
%     [h,tresults(i),~,~] = ttest2(PRESSY(:,i),PRESSY(:,(i+1))); %use two sample t-test to figure out if the current factor is different from Minimum
%     if h == 0 %If the null hypothesis (that the two populations are the same) is retained
%         %We want to find when there isn't a difference in factors any more
%         sigIndiv = [i tresults(i)];
%         sigFactor = vertcat(sigFactor, sigIndiv);
%     end
% end
% 
% %Determine Number of Factors for T-test
% if isempty(sigFactor) == 0 %If this array is not empty 
%     %If the ttest did find factors that were statistically not different
%     %check to see if fewer factors can be used
%     diff = IDm - sigFactor(:,1);
% %     for j = 1:length(sigFactor)
% %         diff = vertcat(diff, abs((sigFactor(j,1)-IDm)));%find the difference between the number of factors determined from the PRESS statistic and the number of factors from the ttes and save to array
% %     end
%     %disp(diff(1,:))
%     [~,I] = min(diff(:,1));%Find the minimum difference, this wil lbe the number of factors closest to the factors determined from teh PRESS statistic
%     if length(sigFactor) > 1 %If there was more than one factor that wasn't statistically significant
%         if abs(sigFactor(I,1) - sigFactor(I-1,1)) == 1 %check to make sure that minmum number of factors determined by ttest doesn't have the factor before it also chosen
%             IDt = sigFactor(I-1,1); %Set the minimum number of factors determined by ttest
%             Mint = meanPRESSY(IDt); %Set the minimum number determined by ttest
%             PCTt = meanPCTVAR(IDt); %Pulls out the percent variation explained at that factor
%         else
%             IDt = sigFactor(I,1);%Set the minimum number of factors determined by ttest
%             Mint = meanPRESSY(IDt); %Set the minimum number determined by ttest
%             PCTt = meanPCTVAR(IDt); %Pulls out the percent variation explained at that factor
%         end
%     elseif size(sigFactor,1) == 1 %If there is only one factor 
%         IDt = sigFactor(1,1);%Set the minimum number of factors determined by ttest
%         Mint = meanPRESSY(IDt); %Set the minimum number determined by ttest
%         PCTt = meanPCTVAR(IDt); %Pulls out the percent variation explained at that factor
%     else
%         disp(['All Factors under ' ncomp ' are p < 0.05'])
%         IDt = NaN(1);
%         Mint = NaN(1);
%         PCTt = NaN(1);
%     end
%     
%     %Display Results from T-test on screen
%     disp('Results from T-test')
%     disp(num2str(sigFactor))
% 
% else %If the array is empty, this means all factors were significantly different from each other
%     % If everything was statistically different, use the minimum PRESS
%     % statistic to determine the number of factors
%     IDt = NaN(1);
%     Mint = NaN(1);
%     PCTt = NaN(1);
% end

%% Slope
%Testing Slope
listSlope = []; %Holds the slope values
IDs = NaN(1); %Holds the number of factors determined by the slope
for i = 2:(ncomp-1) %Loop through factors
    slope = meanPRESSY(:,(i-1)) - meanPRESSY(:,(i)); %Find slope 
    totalSlope = meanPRESSY(1) - meanPRESSY(ncomp); %Determine the complete change in slope
    percentSlope = 100*(slope/totalSlope); %determine the percent change of this one factor compared to the whole thing
    indivSlope = [(i-1) percentSlope];
    listSlope = vertcat(listSlope,indivSlope);
    %disp(['This ' num2str(i) ' factor has a slope of ' num2str(percentSlope)])
%     if percentSlope < 5 %If the slope is less than 1% store this value to be used later
%         disp(['This factor has a slope less than 5% ' num2str(i)])
%         indivSlope = [i percentSlope];
%         listSlope = vertcat(listSlope,indivSlope);
%     end
end

%Go through the minimum slope to figure out which factor to choose
localMinSlope = [];
for j = 1:size(listSlope,1)-1 %Loop through Slope list to check and see if there is a local minimum
    if listSlope(j,2) < listSlope((j+1),2) && meanPCTVAR(j+1) > mean(meanPCTVAR)%If the current index has a slope less then the next  
       localMinSlope = vertcat(localMinSlope,listSlope((j+1),:));
    end 
end



% [Mins,minID] = min(listSlope(:,2));%Set the minimum number determined by slope
% IDs = listSlope(minID,1);%Set the minimum number of factors determined by slope
% PCTs = meanPCTVAR(IDs);%Pulls out the percent variation explained at that factor

% diffSlope = IDm - listSlope;
% [~,I] = min(diffSlope(:,1));%Find the minimum difference, this wil lbe the number of factors closest to the factors determined from teh PRESS statistic
% if listSlope(1,1) >= (IDm - 15) || size(listSlope,1) == 1 %If the slope found is within 15 factors of the factor determined by the min PrESS statistic
%     IDs = listSlope(1,1);%Set the minimum number of factors determined by slope
%     Mins = meanPRESSY(IDs); %Set the minimum number determined by slope
%     PCTs = meanPCTVAR(IDs); %Pulls out the percent variation explained at that factor
% else  % If the first slope is more than 15 factors away from min PRESS statistic..
%     for L = 1: length(listSlope) % Loop through listSlope to find the first slope that is within 15 factors
%         if listSlope(L,1) <= (IDm - 15) 
%             IDs = listSlope(L,1);%Set the minimum number of factors determined by slope
%             Mins = meanPRESSY(IDs); %Set the minimum number determined by slope
%             PCTs = meanPCTVAR(IDs); %Pulls out the percent variation explained at that factor 
%         end
%     end
% end

if isempty(localMinSlope) == 1 %if there were no slopes that fit the above requirements
    [Mins, IDs] = min(listSlope(:,2));
    %IDs %Set the minimum number of factors determined by slope
    %Mins; %Set the minimum number determined by slope
    PCTs = meanPCTVAR(IDs); %Pulls out the percent variation explained at that factor
else
    %Determine number of factors for slope
    [~,tempIndex] = max(localMinSlope(:,2));%Set the minimum number determined by slope
    IDs = localMinSlope(tempIndex,1);%Set the minimum number of factors determined by slope
    Mins = meanPRESSY(IDs);
    PCTs = meanPCTVAR(IDs);%Pulls out the percent variation explained at that factor
end 

%% Display the number of Factors as determined from PRESS, slope, & t-test
disp(['Number of Factors for: ' funcName ' ' spectrumName])
disp([ 'ID of minimum PRESS statistic for Y: ' num2str(IDm)])
disp([ 'Value of minimum PRESS statistic for Y: ' num2str(Minm)])
disp([ 'Percent Variance Explained at minimum PRESS statistic for Y: ' num2str(PCTm)])
% disp([ 'ID from ttest for Y: ' num2str(IDt)])
% disp([ 'Value of minimum using ttest for Y: ' num2str(Mint)])
% disp([ 'Percent Variance Explained at minimum ttest for Y: ' num2str(PCTt)])
disp([ 'ID from Slope Analysis for Y: ' num2str(IDs)])
disp([ 'Value of minimum using Slope Analysis for Y: ' num2str(Mins)])
disp([ 'Percent Variance Explained at minimum Slope Analysis for Y: ' num2str(PCTs)])
%% Set final number of factors
% if isnan(IDt) == 1 %If IDt was set to NaN because all factors were determine significantly different
%     IDt = 20; %Too be graphed correctly.
%     PCTt = meanPCTVAR(IDt);
%     Mint = meanPRESSY(IDt);
%     %Just compare factors determined from PRESS and slope
%     if IDs < IDm
%         ID = IDs;
%         Min = Mins;
%     else
%         ID = IDm;
%         Min = Minm;
%     end
% else
%     %compare factors determined from PRESS, slope, and ttest
%     if IDs < IDm && IDs < IDt
%         ID = IDs;
%         Min = Mins;
%     elseif IDm < IDs && IDm < IDt
%         ID = IDm;
%         Min = Minm;
%     else
%         ID = IDt;
%         Min = Mint;        
%     end
% end

if IDs < IDm %If the # of Factors from slope is less then PRESSRMSEY set overall Factors to slope factors
    ID = IDs;
    Min = Mins;
    PCT = PCTs;
else %Otherwise set them to PRESSRMSEY factors
    ID = IDm;
    Min = Minm;
    PCT = PCTm;
end
% if isnan(IDt) == 1 %If IDt was set to NaN because all factors were determine significantly different
%    ID = IDm;
%    Min = Minm;
% elseif IDm < IDt %If the minimum determined from PRESS statistic is lower than factors determined by Ttest
%    ID = IDm;
%    Min = Minm;
% else %If min factors determined from ttest is lower than factors determined by PRESS statistic
%     ID = IDt;
%     Min = Mint;
% end
       
disp(['The number of factors that will be used: ' num2str(ID)])
%% Figure for % Variance Explained

%Figure of PRESS statistic for y
%Percent Variance Figure
figure('units','normalized','outerposition',[0.5 0 0.5 1])
hold on
plot(1:(ncomp),meanPCTVAR,'-bo'); axis square;
grid on
xlabel('Number of Factors');
ylabel ('Percent of Variance Explained');
axis([0 ncomp 0 100])
scatter(IDm,PCTm,'r','filled'); %Adds minimum onto the graph
%scatter(IDt,PCTt,'b','filled'); %Adds minimum onto the graph
scatter(IDs,PCTs,'g','filled'); %Adds minimum onto the graph
%ttText = ['ttest = ' sprintf('%0.2f',PCTt) ' at ' num2str(IDt)];
pressText = ['PRESS = ' sprintf('%0.2f',PCTm) ' at ' num2str(IDm)];
slopeText = ['Slope = ' sprintf('%0.2f',PCTs) ' at ' num2str(IDs)];
text(0.95,1,{pressText,slopeText},'Units','normalized','VerticalAlignment','top','HorizontalAlignment','right','FontSize',14)
legend('% Var Explained','PRESS statistic','Slope','Location', 'BestOutside');
title([traitName ' ' funcName ' ' spectrumName]);
hold off

%% Figure for PRESS stat

%Figure of PRESS statistic for y
%PRESS Statistic Variable
figure('units','normalized','outerposition',[0 0 0.5 1])
hold on
plot(1:(ncomp),meanPRESSY,'-bo'); axis square;
grid on
xlabel('Number of Factors');
ylabel ('PRESS Statistic');
axis([0 ncomp 0 meanPRESSY(1)])
scatter(IDm,Minm,'r','filled'); %Adds minimum onto the graph
%scatter(IDt,Mint,'b','filled'); %Adds minimum onto the graph
scatter(IDs,Mins,'g','filled'); %Adds minimum onto the graph
%ttText = ['ttest = ' sprintf('%0.2f',Mint) ' at ' num2str(IDt)];
pressText = ['PRESS = ' sprintf('%0.2f',Minm) ' at ' num2str(IDm)];
slopeText = ['Slope = ' sprintf('%0.2f',Mins) ' at ' num2str(IDs)];
text(0.95,1,{pressText,slopeText},'Units','normalized','VerticalAlignment','top','HorizontalAlignment','right','FontSize',14)
legend('PRESS Statistic','Min PRESS','Slope','Location', 'BestOutside');
title([traitName ' ' funcName ' ' spectrumName]);
hold off

%% Prompt For Input
prompt = 'Number of factors to be used?';
ID = input(prompt)
%% Old Way of Figuring out Factors
% Figure out the lowest number of components using PRESS 
% meanPRESSY = mean(PRESSY);
% meanPCTVAR = mean(PCTVAR);
% IDm = [];
% IDt = [];
% 
% [Minm,IDm] = min(meanPRESSY(10:20)); % Pulls lowest PRESS statistic to determine number of factors
% IDm = IDm + 9;
% PCTm = meanPCTVAR(IDm);
%         
% disp([ 'ID of minimum PRESS statistic for Y: ' num2str(IDm)])
% disp([ 'Value of minimum PRESS statistic for Y: ' num2str(Minm)])
% disp([ 'Percent Variance Explained at minimum PRESS statistic for Y: ' num2str(PCTm)])
% tresults = zeros(1,ncomp);
% temp = [];
% for i = 10:20
%     [~,tresults(i),~,~] = ttest2(PRESSY(i,:),PRESSY((i+1),:));
%     if tresults(i) > 0.05
%         temp = vertcat(temp,[i,tresults(i)]);
%     end
% end
% 
% disp('Results from T-test')
% disp(num2str(temp))
% 
% diffTP = [];
% for j = 1:size(temp,1)
%     diffTP = vertcat(diffTP,abs((temp(j,1)-IDm)));
% end
% 
% [~,I] = min(diffTP);
% if size(temp,1) > 1
%     if abs(temp(I,1) - temp(I-1,1)) == 1
%         IDt = temp(I-1,1);
%     else
%         IDt = temp(I,1);
%     end
% elseif size(temp,1) == 1
%     IDt = temp(1,1);
% else
%     disp('All Factors under 20 are p < 0.05')
%     IDt = IDm;
% end
% 
% Mint = meanPRESSY(IDt);
% PCTt = meanPCTVAR(IDt);
% 
% disp([ 'ID from ttest for Y: ' num2str(IDt)])
% disp([ 'Value of minimum using ttest for Y: ' num2str(Mint)])
% disp([ 'Percent Variance Explained at minimum ttest for Y: ' num2str(PCTt)])

% Set final number of factors
% if 10<=IDt && IDt<=20 || 10<=IDm && IDm<=20
%     if IDm < IDt
%         ID = IDm;
%         Min = Minm;
%     else
%         ID = IDt;
%         Min = Mint;
%     end
% elseif isnan(IDt) == 1
%     ID = IDm;
%     Min = Minm;
% else
%     disp('Number of Factors are not between 10 and 20')
%     ID = IDt;
%     Min = Mint;
% end
% disp(['The number of factors that will be used: ' num2str(ID)])
end
