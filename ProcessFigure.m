function ProcessFigure(hFigure, fileName, varargin)
% ProcessFigure(hFigure, fileName, varargin)
% [fontSize, paperSize]
% History:
% Dmytro Velychko - Created
% Shrisha - Modified 
    
    [fontSize, paperSize, IF_SAVE] = DefaultArgs(varargin, {8, [4.4, 4.4], 1});
    allAxesInFigure = findall(hFigure, 'type', 'axes');
    for hAxis = allAxesInFigure'
        set(hAxis, 'box', 'off');
        set(hAxis, 'TickDir', 'out');
        set(hAxis,'FontSize', fontSize);
        set(get(hAxis,'XLabel'), 'FontSize', fontSize);
        set(get(hAxis,'YLabel'), 'FontSize', fontSize);
        txtHdl = findall(hAxis, 'type', 'text');
        for kTxtHdl = txtHdl
            set(kTxtHdl, 'fontSize', fontSize);
        end
    end
    set(hFigure, 'PaperPosition', [0, 0, paperSize]);
    set(hFigure, 'PaperSize', paperSize);
    set(hFigure, 'Renderer', 'Painters');
    saveas(hFigure, [fileName '.pdf'], 'pdf');
    if IF_SAVE
        saveas(hFigure, [fileName, '.fig'], 'fig');
    end
end