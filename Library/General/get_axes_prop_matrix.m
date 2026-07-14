function limMat = get_axes_prop_matrix(axH, propName)

axH = axH(:);
axH = axH(isgraphics(axH));

if isempty(axH)
    limMat = [];
    return;
end

raw = get(axH, propName);

if iscell(raw)
    limMat = cell2mat(raw);
else
    limMat = raw;
end

end