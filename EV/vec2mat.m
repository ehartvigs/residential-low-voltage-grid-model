function mat = vec2mat(vec, n)
%VEC2MAT Replaces vec2mat(vec, n) from the Communications Toolbox
%   mat = vec2mat(vec,n) converts the vector vec into a matrix with
%   n columns, creating one row at a time. If the length of vec is not
%   a multiple of n, then extra zeros are placed in the last row of
%   mat. The matrix mat has ceil(length(vec)/n) rows.
numPadded = mod(numel(vec),n);
if numPadded > 0
    numPadded = n - numPadded
    mat = reshape([vec.' padding(1:numPadded)], n, []).'
else
    numPadded % No padding required
    mat = reshape(vec.', n, []).'
end

end

