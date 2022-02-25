function MAF_Beautify_Ax(hAx, font_size)

if nargin < 1
    if	isempty(get(0,'currentfigure'))
        disp('ESN_Beautify_Plot :: ERROR, no figure');
        return;
    end
    hAx = gcf;
end

if nargin < 2
    font_size = 12;
end


font_name = 'Arial'; % 'Courier New'; %

hAx.FontName      = font_name;
hAx.FontUnits     = 'points';
hAx.FontSize      = font_size;
hAx.Box           = 'off';
hAx.TickDir       = 'out';
hAx.TickLength    = [.02 .02];
try
    hAx.XMinorTick    = 'off';
    hAx.YMinorTick    = 'off';
    hAx.XGrid         = 'off';
    hAx.YGrid         = 'off';
    hAx.GridLineStyle = '--';
    hAx.XLabel.FontName  = font_name;
    hAx.XLabel.FontUnits = 'points';
    hAx.XLabel.FontSize  = font_size;
    hAx.YLabel.FontName  = font_name;
    hAx.YLabel.FontUnits = 'points';
    hAx.YLabel.FontSize  = font_size;
catch
end
%     hAx.XColor        = [0 0 0];
%     hAx.YColor        = [0 0 0];
hAx.LineWidth     = 1;
hAx.Clipping      = 'off';

hAx.Title.FontName   = font_name;
hAx.Title.FontUnits  = 'points';
hAx.Title.FontSize   = font_size;
hAx.Title.FontWeight = 'bold';
hold(hAx, 'off');

end