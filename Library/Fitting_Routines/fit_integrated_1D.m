function PCell = fit_integrated_1D( ...
    analyVar, ...
    OD_Fit_ImageCell, ...
    InitGuess, ...
    lowBndVec, ...
    upBndVec, ...
    optimOpt)

% Fit X and Y integrated profiles independently.
%
% The resulting parameters are converted back into the normal
% 2D PCell parameter convention so that all downstream code
% remains unchanged.
%
% Temporary X-profile parameters:
%   [AmpX, xCntr, sigX, OffsetX, SlopeX]
%
% Temporary Y-profile parameters:
%   [AmpY, yCntr, sigY, OffsetY, SlopeY]
%
% Final PCell remains:
%   [Amp, sigX, sigY, xCntr, yCntr, Offset, SlopeX, SlopeY]

PCell = cell(size(OD_Fit_ImageCell));


for n = 1:numel(OD_Fit_ImageCell)

    %% ============================================================
    %  GET IMAGE
    % ============================================================

    OD = ...
        double(OD_Fit_ImageCell{n});


    nY = size(OD,1);
    nX = size(OD,2);


    %% ============================================================
    %  INTEGRATE IMAGE
    % ============================================================

    % Integrate over Y -> X profile.

    xData = ...
        sum(OD,1)';


    % Integrate over X -> Y profile.

    yData = ...
        sum(OD,2);


    % Pixel coordinates.

    xAxis = ...
        (1:nX)';

    yAxis = ...
        (1:nY)';


    %% ============================================================
    %  GET EXISTING INITIAL GUESS
    % ============================================================

    p2D = ...
        InitGuess{n};


    parameterNames = ...
        analyVar.InitCondList;


    %% ============================================================
    %  FIND 2D PARAMETER INDICES
    % ============================================================

    ampIndex = ...
        find(strcmpi(parameterNames,'Amp'),1);

    sigXIndex = ...
        find(strcmpi(parameterNames,'sigX'),1);

    sigYIndex = ...
        find(strcmpi(parameterNames,'sigY'),1);

    xCntrIndex = ...
        find(strcmpi(parameterNames,'xCntr'),1);

    yCntrIndex = ...
        find(strcmpi(parameterNames,'yCntr'),1);

    offsetIndex = ...
        find(strcmpi(parameterNames,'Offset'),1);

    slopeXIndex = ...
        find(strcmpi(parameterNames,'SlopeX'),1);

    slopeYIndex = ...
        find(strcmpi(parameterNames,'SlopeY'),1);


    %% ============================================================
    %  X PROFILE INITIAL GUESS
    % ============================================================

    sigmaX0 = ...
        p2D(sigXIndex);

    xCntr0 = ...
        p2D(xCntrIndex);


    % The integrated amplitude is approximately the peak 2D
    % amplitude multiplied by sqrt(2*pi)*sigmaY.

    ampX0 = ...
        p2D(ampIndex) .* ...
        sqrt(2*pi) .* ...
        p2D(sigYIndex);


    % Integrated constant background.

    offsetX0 = ...
        nY .* ...
        p2D(offsetIndex);


    % Integrated X slope.

    slopeX0 = ...
        nY .* ...
        p2D(slopeXIndex);


    pX0 = [
        ampX0
        xCntr0
        sigmaX0
        offsetX0
        slopeX0
        ];


    %% ============================================================
    %  Y PROFILE INITIAL GUESS
    % ============================================================

    sigmaY0 = ...
        p2D(sigYIndex);

    yCntr0 = ...
        p2D(yCntrIndex);


    % Integrated amplitude.

    ampY0 = ...
        p2D(ampIndex) .* ...
        sqrt(2*pi) .* ...
        p2D(sigXIndex);


    % Integrated constant background.

    offsetY0 = ...
        nX .* ...
        p2D(offsetIndex);


    % Integrated Y slope.

    slopeY0 = ...
        nX .* ...
        p2D(slopeYIndex);


    pY0 = [
        ampY0
        yCntr0
        sigmaY0
        offsetY0
        slopeY0
        ];


    %% ============================================================
    %  PROFILE BOUNDS
    % ============================================================

    % Only the amplitude, center, and sigma need physical
    % constraints. Background terms remain unconstrained.

    xLower = [
        0
        1
        0
        -Inf
        -Inf
        ];

    xUpper = [
        Inf
        nX
        Inf
        Inf
        Inf
        ];


    yLower = [
        0
        1
        0
        -Inf
        -Inf
        ];

    yUpper = [
        Inf
        nY
        Inf
        Inf
        Inf
        ];


    %% ============================================================
    %  FIT X PROFILE
    % ============================================================

    pX = ...
        lsqcurvefit( ...
            @Gaussian1D, ...
            pX0, ...
            xAxis, ...
            xData, ...
            xLower, ...
            xUpper, ...
            optimOpt);


    %% ============================================================
    %  FIT Y PROFILE
    % ============================================================

    pY = ...
        lsqcurvefit( ...
            @Gaussian1D, ...
            pY0, ...
            yAxis, ...
            yData, ...
            yLower, ...
            yUpper, ...
            optimOpt);


    %% ============================================================
    %  EXTRACT FITTED PARAMETERS
    % ============================================================

    fittedSigX = ...
        pX(3);

    fittedSigY = ...
        pY(3);

    fittedXCntr = ...
        pX(2);

    fittedYCntr = ...
        pY(2);


    %% ============================================================
    %  CONVERT PROFILE AMPLITUDES BACK TO 2D AMPLITUDE
    % ============================================================

    ampFromX = ...
        pX(1) ./ ...
        (sqrt(2*pi) .* fittedSigY);


    ampFromY = ...
        pY(1) ./ ...
        (sqrt(2*pi) .* fittedSigX);


    % Both profiles should give approximately the same 2D
    % amplitude. Their average is a simple symmetric estimate.

    fittedAmp = ...
        mean([ampFromX ampFromY]);


    %% ============================================================
    %  CONVERT PROFILE BACKGROUNDS TO 2D BACKGROUND
    % ============================================================

    fittedOffset = ...
        mean([ ...
            pX(4)./nY
            pY(4)./nX
            ]);


    fittedSlopeX = ...
        pX(5)./nY;


    fittedSlopeY = ...
        pY(5)./nX;


    %% ============================================================
    %  INSERT RESULTS INTO NORMAL 2D PARAMETER VECTOR
    % ============================================================

    p2D(ampIndex) = ...
        fittedAmp;

    p2D(sigXIndex) = ...
        fittedSigX;

    p2D(sigYIndex) = ...
        fittedSigY;

    p2D(xCntrIndex) = ...
        fittedXCntr;

    p2D(yCntrIndex) = ...
        fittedYCntr;

    p2D(offsetIndex) = ...
        fittedOffset;

    p2D(slopeXIndex) = ...
        fittedSlopeX;

    p2D(slopeYIndex) = ...
        fittedSlopeY;


    %% ============================================================
    %  STORE NORMAL 2D PCell
    % ============================================================

    PCell{n} = ...
        p2D;

end

end

function z = Gaussian1D(p,x)

z = ...
    p(1) .* ...
    exp( ...
        -(x-p(2)).^2 ./ ...
        (2*p(3)^2)) ...
    + p(4) ...
    + p(5).*x;

end