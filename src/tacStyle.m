function S = tacStyle()
%TACSTYLE  Single source of truth for TAC paper figure styling.
%   Usage (top of any figure script):
%       S = tacStyle;                 % constants + LaTeX (CM) interpreters
%       ... plot ...
%       tacAxes(ax, S);               % apply shared axis styling
%
%   Sizes are tuned for a 3.5 in (single IEEE column) wide export. Include
%   each PDF in LaTeX at its natural size (no width scaling) so the fonts
%   match across every figure.
%
%   Also flips the groot defaults to the LaTeX interpreter so tick labels
%   render in Computer Modern, matching the axis/legend text.

    % --- interpreters: Computer Modern everywhere ------------------------
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    % --- palette (Okabe–Ito, colorblind-safe) ----------------------------
    S.cE    = [0   114 178]/255;   % blue         (East)
    S.cN    = [213  94   0]/255;   % vermillion   (North)
    S.cU    = [0   158 115]/255;   % bluish green (Up)
    S.blue  = S.cE;                % accent blue reused in schematics
    S.green = [0.85 1.00 0.85];    % light region fill (target set A)
    S.grey1 = 1 - [0.2 0.2 0.2];   % lighter hatch band
    S.grey2 = 1 - [0.1 0.1 0.1];   % darker  hatch band

    % --- sizing ----------------------------------------------------------
    S.colW   = 3.5;                % IEEE single-column width [in]
    S.panelW = 0.48*S.colW;        % one cell of a 2x2 grid within a column

    % --- line widths -----------------------------------------------------
    S.lwMain = 0.75;               % primary traces / strokes
    S.lwThin = 0.5;                % boundaries / secondary lines
    S.lwLead = 0.35;               % thin leader lines to labels
    S.axLW   = 0.5;                % axis rule width

    % --- font sizes [pt] (for a 3.5 in wide export) ----------------------
    S.fsTick  = 8;                 % tick labels
    S.fsLab   = 9;                 % axis labels
    S.fsLeg   = 8;                 % legends
    S.fsAnnot = 8;                 % in-figure math labels (half-column schematic panels)
    S.fsSmall = 7;                 % units / band captions

    % --- highlight band --------------------------------------------------
    S.ffColor = [0.80 0.84 0.90];
    S.ffAlpha = 0.35;
end