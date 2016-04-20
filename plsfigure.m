function plsfigure(valTrait, valRegLine,valCoeff,valRsq,valRMSE,traitName,CIpred,CIsample,valSeason,funcName,spectrumName)
%CImodel,
%% Setting up Variables

% Set the upper limit of figures
if roundn(max(valTrait),2) > roundn(max(valTrait),0) %If the maximum value of trait to the hundreth is larger than the maximum trait value to the integer
    upperLim = roundn(max(valTrait),2);%Set the upper limit to the maximum trait value rounded to the hundreth
else % If not...
    upperLim = roundn(max(valTrait),0); %Set the upper limit to the maximum trait value rounded to the nearest integer
end

if roundn(min(valTrait),2) < roundn(min(valTrait),0)
    lowerLim = roundn(min(valTrait),2);
else
    lowerLim = roundn(min(valTrait),0);
end

if strcmp(traitName,'LMA') == 1
    units = '(g/m^2)';
else
    units = '(%)';
end

%% Valibration figure gray scale
figure; 
axis square;
hold on
%Add regression line
hf = @(x) valRegLine(1) + valRegLine(2)*x;
he = ezplot(hf,[0 upperLim]);
set(he,'Color','k','LineWidth',1.5)

%CI for PLSR models
hciup = @(x) CIpred(1) + valRegLine(2)*x;
hciupf = ezplot(hciup,[0 upperLim]);
set(hciupf,'Color',[119/256 136/256 153/256],'LineWidth',1.5)
hcidw = @(x) CIpred(2) + valRegLine(2)*x;
hcidwf = ezplot(hcidw,[0 upperLim]);
set(hcidwf,'Color',[119/256 136/256 153/256],'LineWidth',1.5)

%Add 1:1 reference line
hRefLine = refline(1,0);
set(hRefLine,'Color','k','LineStyle',':','LineWidth',1.5);

%CI for validation samples
% CIsampleInput = ((CIsample(1,:) - CIsample(2,:))/2)';
% errorbar(valCoeff, valTrait, CIsampleInput, 'k.')
    
%Add scatter plot of valibration samples
for season = 1:3
    indexSeason = find(valSeason == season);
  
    %Scatterplot
    s = scatter(valCoeff(indexSeason),valTrait(indexSeason),60,'filled');
    if season == 1
        color = [255/256 255/256 255/256]; %Spring
        outline = [0/256 0/256 0/256];
    elseif season == 2
        color = [119/256 136/256 153/256]; %Summer
        outline = [112/256 128/256 144/256];
    else
        color = [0/256 0/256 0/256]; %Fall
        outline = [0/256 0/256 0/256];
    end
    set(s, 'MarkerFaceColor',color,'MarkerEdgeColor',outline);
end

%Add additional features of graph
set(gca,'FontSize',14)
set(gca,'YTick',[lowerLim:(upperLim*.25):upperLim])
set(gca,'XTick',[lowerLim:(upperLim*.25):upperLim])
set(gca,'YLim',[lowerLim upperLim])
set(gca,'XLim',[lowerLim upperLim])
ylabel(['Observed ' traitName ' ' units],'FontSize',16);
xlabel(['Predicted ' traitName ' ' units],'FontSize',16);
rsqText = ['R^2 = ' sprintf('%0.2f',valRsq)];
rmseText = ['RMSE = ' sprintf('%0.2f',valRMSE)];
biasText = ['Bias = ' sprintf('%0.2f',valRegLine(1))]; %Model intercept
nText = ['N_v_a_l = ' num2str(size(valTrait,1))];
text(0.05,1,{rsqText,rmseText,biasText, nText},'Units','normalized','VerticalAlignment','top','FontSize',14)
%title([funcName ' ' spectrumName]);
title('');
hold off
end