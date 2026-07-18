function [xdata, ydata] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset, tweezerNum)
    
    if analyVar.UseTweezer
        tweezerUse = 1;
    else 
        tweezerUse = 0;
    end
    if nargin < 6
        tweezerNum = 1; % Default value if TweezerNum
        tweezerUse = 0;
    end

    xdata = cell(analyVar.numBasenamesAtom,1);                              % unfiltered data X initialized
    ydata = cell(analyVar.numBasenamesAtom,1);                              % unfiltered data Y initialized
    
    if tweezerUse
        %% Extracting x,y data for each BasenamesAtom
        for i = 1:analyVar.numBasenamesAtom
        
            xdata{i} = indivDataset{i}.(indVarField) + analyVar.DropTimeOffset;
            ydata{i} = indivDataset{i}.(depVarField)(:,tweezerNum);
    
        end  
    else
    
        %% Extracting x,y data for each BasenamesAtom
        for i = 1:analyVar.numBasenamesAtom
        
            xdata{i} = indivDataset{i}.(indVarField) + analyVar.DropTimeOffset;
            ydata{i} = indivDataset{i}.(depVarField);
    
        end
    end
    
end