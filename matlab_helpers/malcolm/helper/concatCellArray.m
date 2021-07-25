function y = concatCellArray(x)
%concatCellArray.m Concatenates cell array into an array

y = [];
for i = 1:numel(x)
    y = [y; x{i}];
end

end

