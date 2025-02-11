function [xdata, ydata] = getxy(indVarField, depVarField, analyVar, indivDataset, avgDataset)

    xdata = cell(analyVar.numBasenamesAtom,1);                              % unfiltered data X initialized
    ydata = cell(analyVar.numBasenamesAtom,1);                              % unfiltered data Y initialized
    

    
    %% Extracting x,y data for each BasenamesAtom
    for i = 1:analyVar.numBasenamesAtom
        
        xdata{i} = indivDataset{i}.(indVarField);
        ydata{i} = indivDataset{i}.(depVarField);

    end
    
end