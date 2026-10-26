function tacAxes(ax, S)
%TACAXES  Apply the shared TAC axis styling to axes AX using struct S.
    grid(ax,'on');
    box(ax,'on');
    ax.GridAlpha  = 0.15;
    ax.LineWidth  = S.axLW;
    ax.FontSize   = S.fsTick;
    ax.TickLabelInterpreter = 'latex';
    ax.XMinorGrid = 'off';
end