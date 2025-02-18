% Example Usage:
myCellArray = {1, 'hello', 3.14, true; [1 2 3], {'inner' 'cell'}, NaN, pi};
myStringMatrix = cellArrayToStringMatrix(myCellArray);
disp(myStringMatrix);
disp(class(myStringMatrix)); % Display the class (char array)
disp(' ')

emptyCellArray = {};
emptyStringMatrix = cellArrayToStringMatrix(emptyCellArray);
disp(emptyStringMatrix);
disp(class(emptyStringMatrix)); % Display the class (char array)
disp(' ')


notACell = 5;
notACellString = cellArrayToStringMatrix(notACell);
disp(notACellString);
disp(class(notACellString)); % Display the class (char array)
disp(' ')


function stringMatrix = cellArrayToStringMatrix(cellArray)
% Converts a cell array to a matrix of strings.
%
% Args:
%   cellArray: The input cell array.
%
% Returns:
%   A matrix of strings, where each element corresponds to the string 
%   representation of the corresponding cell in the input.  Returns an
%   empty character array if the input is not a cell array. Handles
%   different data types within the cells.

if ~iscell(cellArray)
    warning('Input is not a cell array. Returning empty string.');
    stringMatrix = char.empty(); % Return empty char array
    return;
end

if isempty(cellArray)
    stringMatrix = char.empty();
    return;
end


[rows, cols] = size(cellArray);  % Get dimensions of the cell array
stringMatrix = cell(rows, cols); % Preallocate the cell array of strings

for i = 1:rows
    for j = 1:cols
        currentValue = cellArray{i, j};

        if ischar(currentValue)
            stringMatrix{i, j} = currentValue;
        elseif isnumeric(currentValue)
            stringMatrix{i, j} = num2str(currentValue);
        elseif islogical(currentValue)
            stringMatrix{i, j} = num2str(currentValue);
        else
            warning(['Unsupported data type at (' num2str(i) ', ' num2str(j) '). Converting to string ''Unsupported''.']);
            stringMatrix{i, j} = 'Unsupported';  % Or handle differently
        end
    end
end

% stringMatrix = char(stringMatrix); % Convert cell array of strings to char array

end % end of function