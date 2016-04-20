%Susan Meerdink
%This function calculates the number of factors for PLSR
function [valIndexAll, valIndexBroad,valIndexNeedle,valIndexSpring,valIndexSummer,valIndexFall] = splitValCalBefore(allTrait,FuncType,Season)

%% Loop through Functional Types
    for functype = 1:6
        valIndex = [];
        if functype == 1
            %Broadleaf
            trait = allTrait(find(FuncType == 1 | FuncType == 2));
            nanListFull = find(isnan(trait)); %Finds all NaN values in trait dataset being analyzed
            trait(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
            indexList = 1:1:(size(trait,1));
        elseif functype == 2
            %Needleaf
            trait = allTrait(find(FuncType == 3));
            nanListFull = find(isnan(trait)); %Finds all NaN values in trait dataset being analyzed
            trait(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
            indexList = 1:1:size(trait,1);
        elseif functype == 3
            %Spring
            trait = allTrait(find(Season == 1));
            nanListFull = find(isnan(trait)); %Finds all NaN values in trait dataset being analyzed
            trait(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
            indexList = 1:1:(size(trait,1));
        elseif functype == 4
            %Summer
            trait = allTrait(find(Season == 2));
            nanListFull = find(isnan(trait)); %Finds all NaN values in trait dataset being analyzed
            trait(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
            indexList = 1:1:(size(trait,1));
        elseif functype == 5
            %Fall
            trait = allTrait(find(Season == 3));
            nanListFull = find(isnan(trait)); %Finds all NaN values in trait dataset being analyzed
            trait(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
            indexList = 1:1:(size(trait,1));
        else
            %All Samples
            trait = allTrait;
            nanListFull = find(isnan(trait)); %Finds all NaN values in trait dataset being analyzed
            trait(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
            indexList = 1:1:(size(trait,1));
        end %end of if, elseif, and else choosing functional types
        
        %% Finding and removing NaN Values
%         nanListFull = find(isnan(trait)); %Finds all NaN values in trait dataset being analyzed
%         trait(nanListFull) = []; %Removes all NaN values from trait dataset being analyzed
%         indexList(nanListFull) = []; %Removes all NaN values from dataset being analyzed
       
        %% Randomly select from each bin for a total of 20% of data
        if max(trait)< 1
            step = .1;
        elseif max(trait) > 100
            if size(trait,1) < 100
                step = (size(trait,1)/3);
            else
                step = (size(trait,1)/10);
            end
        else
            step = (size(trait,1)/100);
        end
        binranges = min(trait):step:max(trait);
        [bincounts,binindex] = histc(trait,binranges);
        rng('default')
        
        for num = 1:size(bincounts,1)
            if bincounts(num) < 0
                continue
            else
                inIndex = indexList(find(binindex == num))';
                random = randperm((bincounts(num)),round((.2*bincounts(num)))); %Get random numbers in this bin range
                valIndex = vertcat(valIndex, inIndex(random));%Add them to the validation list
            end
        end
        
        %% Look for duplicates and replace if necessary
        checkDup = unique(valIndex,'sorted');
        count = 1;
        while size(checkDup,1) < floor(.2*size(trait,1))
            random = randperm((size(trait,1)),(size(valIndex,1) - size(checkDup,1)));
            valIndex = checkDup; %Remove duplicate
            valIndex = vertcat(valIndex, inIndex(random));%Add them to the validation list
            checkDup = unique(valIndex);
            if count == 100 %If it hasn't been able to find a number that isn't a duplicate by 100 iterations quit and move on
                break
            end
            count = count +1;
        end
        
        if functype == 1
            %Broadleaf
            valIndexBroad = sort(valIndex);
        elseif functype == 2
            %Needleaf
            valIndexNeedle = sort(valIndex);
        elseif functype == 3
            %Spring
            valIndexSpring = sort(valIndex);
        elseif functype == 4
            %Summer
            valIndexSummer = sort(valIndex);
        elseif functype == 5
            %Fall
            valIndexFall = sort(valIndex);
        else
            %All Samples
            valIndexAll = sort(valIndex);
        end %end of if, elseif, and else choosing functional types
        
    end
end