%roots finder in interval

function root = select_root(r, lower_bound, upper_bound)

    tol = 1e-3;

    % Keep roots that are approximately real
    valid_roots = r(abs(imag(r)) < tol);

    % Extract real part
    valid_roots = real(valid_roots);

    % Keep roots inside the desired interval
    valid_roots = valid_roots(valid_roots >= lower_bound - tol & valid_roots <= upper_bound + tol);

    % If no suitable root is found
    if isempty(valid_roots)
    error('No valid root found within the specified interval.');
    end

    % If multiple roots exist, select the smallest one
     root = valid_roots;

end