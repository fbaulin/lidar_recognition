% -----------------------------
function datatipTxt = gscatterDatatipCallback(obj,evt,xnam,ynam)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

group = getappdata(target,'group');
groupname = getappdata(target,'groupname');
gind = getappdata(target,'gind');

if isempty(xnam)
    xnam = 'x';
end
if isempty(ynam)
    ynam = 'y';
end
if isempty (gind)
    % One group
    % Leave group name alone, it may be empty
    % Line index number is the same as the original row
    obsind = ind;
else
    % Multiple groups
    % If group name not given, assign it its number
    if isempty(groupname)
        groupname = num2str(group);
    end
    % Map line index to the original row
    obsind = gind(ind);
end

[xVal,yVal] = matlab.graphics.internal.makeNonNumeric(obj,pos(1),pos(2));
xVal = convertToString(xVal);
yVal = convertToString(yVal);

datatipTxt = {...
    [xnam ': ' xVal]...
    [ynam ': ' yVal]...
    ''...
    getString(message('stats:gscatter:Observation',num2str(obsind)))
    };

if ~isempty(groupname)
    datatipTxt{end+1} = getString(message('stats:gscatter:Group',groupname));
end

function str = convertToString(val)
if(isnumeric(val))
    str = num2str(val);
else
    str = char(string(val));
end