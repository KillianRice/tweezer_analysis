function [xdata_clean, ydata_clean] = getxy_filtered(indVarField, depVarField, analyVar, indivDataset, avgDataset)

    xdata = cell(analyVar.numBasenamesAtom,1);                              % unfiltered data X initialized
    ydata = cell(analyVar.numBasenamesAtom,1);                              % unfiltered data Y initialized
    xdata_clean = cell(analyVar.numBasenamesAtom,1);                        % filtered data X initialized
    ydata_clean = cell(analyVar.numBasenamesAtom,1);                        % filtered data X initialized

    %% Adding DAQ voltages to the indivDataset struct.
    numchannels = 8;
    time_axis = 1;                                                          % if 1 plots against time axis, else plots against imagevcoatom
    use_channels = [1 1 1 1 1 0 1 0];                                       % which channels to plot
    channel_names = {'922nmMOTCavityPD (V)',...                             % AI 0
                    '461nmZeemanPD (V)',...                                 % AI 1
                    '413nm_monPD',...                                       % AI 2
                    '461nm_MOTPD',...                                       % AI 3
                    '461nm_MOTcavityPD',...                                 % AI 4
                    '',...                                                  % AI 5
                    '826nm_TransmissionPD (V)',...                          % AI 6
                    '',...                                                  % AI 7
                    };
    daq_line_format = '%{yyyy.MM.dd-HH:mm:ss}D%f%f%f%f%f%f%f%f%f';
    legend_logic = zeros(size(analyVar.basenamevectorAtom));
    % For loop over filenames in batch file.
    for i = 1:analyVar.numBasenamesAtom
        
        daq_file = [analyVar.dataDir char(analyVar.basenamevectorAtom{i}) '_daq_voltages.dat'];
        
        if exist(daq_file, 'file') == 2
            df = fopen(daq_file);
            indivDataset{i}.daq_voltages = textscan(df, daq_line_format, 'headerLines',1);
            legend_logic(i) = 1;
        else
            disp(strcat(analyVar.basenamevectorAtom{i}, ' daq file not found.'));
            indivDataset{i}.daq_voltages = {};
        end
        
    end

    %% Extracting x,y data for each BasenamesAtom filtered by the conditions set by find. 
    % Mainly use PD voltages to look for unlocked lasers and cavities.
    for i = 1:analyVar.numBasenamesAtom
        
        xdata{i} = indivDataset{i}.(indVarField);
        ydata{i} = indivDataset{i}.(depVarField);
        % Filter Conditions using 'find'
        filtered_idx = find(indivDataset{i}.daq_voltages{3}>analyVar.IR_MOTCavPD(1) & indivDataset{i}.daq_voltages{3}<analyVar.IR_MOTCavPD(2)... % 922nm power drift filter.
            & indivDataset{i}.numberAtom < analyVar.AtomNumLims(2) & indivDataset{i}.numberAtom > analyVar.AtomNumLims(1)... % Atom number filter.
            & indivDataset{i}.daq_voltages{4}>analyVar.ZeemanPD(1) & indivDataset{i}.daq_voltages{4}<analyVar.ZeemanPD(2)... % ZeemanBeam PD filter.
            & indivDataset{i}.daq_voltages{5}>analyVar.SpecBeam2Threshold); % 826 nm ULE lock filter from Transmission PD.
        indivDataset{i}.filtered_idx = filtered_idx;
        xdata_clean{i} = xdata{i}(filtered_idx);
        ydata_clean{i} = ydata{i}(filtered_idx);
    end
    
end