
function [v, t] = find_poly(oyll)
%[Description]

%[Parameters]
%oyll - This is the '*', which was previously introduced in task-detail
%This function is to find the value of the characteristic polynomial by
%recursive method. In 1st step, I read matrix A from file. And then in
%step 2, I ought to judge if matrix A is a sqare matrix and a Triple
%diagonal matrix. In step 3,
%[Returns]
%v - vector of coefficients
%t - run time

%document cputime t1
t1 = cputime;

%********task 3 (a)********
%********step 1********
file_name = sprintf("data/Amat%d.m", oyll);
run(file_name) %run the resource file and load data by the way

%********step 2********
[m, n] = size(A);
%judge if it is square matrix
if m ~= n
    error('The given matrix is not square');
end

%judge if it is tridiagonal matrix
if ~isbanded(A,1,1)
    error('The given matrix is not tridiagonal');
end

%********step 3********
% Calculate the characteristic polynomial of a tridiagonal matrix
% represented by the matrix T.
T = A;
matrix_size = size(T, 1);

% Extract the subdiagonal and diagonal of the matrix T.
subdiag = diag(T(2:end, 1:end-1));
main_diag = diag(T);

% Initialize the coefficients of the characteristic polynomial.
characteristic_poly = zeros(matrix_size+1, matrix_size+1);
characteristic_poly(1, 1) = 1;
characteristic_poly(2, 1:2) = [-main_diag(1), 1];

% Iterate over the diagonal of the matrix T to compute the remaining
% coefficients of the characteristic polynomial.
for iteration = 2:matrix_size
    % Compute the reversed cumulative products of the subdiagonal entries
    % up to the current iteration.
    reversed_subdiag_products = flipud(cumprod(flipud(subdiag(1 : iteration-1))));
    
    % Extract the subdiagonal values and previous polynomial coefficients
    % up to the current iteration.
    subdiag_values = (T(1:iteration-1, iteration));
    prev_poly = characteristic_poly(1:iteration-1, :);
    
    % Compute the recursive term and convolve the previous polynomial
    % with the latest subdiagonal and diagonal entries.
    recursive_term = sum(bsxfun(@times, reversed_subdiag_products.*subdiag_values, prev_poly), 1);
    conv_poly = my_conv(characteristic_poly(iteration-1+1, :), [-main_diag(iteration), 1]);
    
    % Subtract the recursive term from the convolved polynomial to
    % obtain the coefficients of the latest iteration's polynomial.
    characteristic_poly(iteration+1, :) = conv_poly(1:end-1) - recursive_term;
end

% Reverse the order of the coefficients to obtain the characteristic
% polynomial in the correct order.
v = fliplr(characteristic_poly(end, :));

%********step 4********
v_file = sprintf('result/vvec%d.m', oyll);
fid = fopen(v_file, 'w');
%write v by following steps
if isreal(v)
    %if it is a real matrix
    fprintf(fid,['v = [',repmat('%f;',1,n+1),'];'], v');
else
    %if it is not
    fprintf(fid,['v = complex([',repmat('%f;',1,n+1),'],'], real(v)');
    fprintf(fid,['[',repmat('%f;',1,n+1),']);'], imag(v)');
end
fclose(fid);

%document cputime t2
t2 = cputime;
%compute running time t
t = t2 - t1;

%%********task 3 (b)********
if oyll == 11
    task3b_file = 'result/task_3b.txt';
    fid = fopen(task3b_file, 'w');
    fprintf(fid, "run time = %8.4f s \n\n", t);
    %write v by following steps
    if isreal(v)
        %if it is a real matrix
        fprintf(fid,['answer:\nv = [',repmat('%f;',1,n+1),'];'], v');
    else
        %if it is not
        fprintf(fid,['answer:\nv = complex([',repmat('%f;',1,n+1),'],'], real(v)');
        fprintf(fid,['[',repmat('%f;',1,n+1),']);'], imag(v)');
    end
    fclose(fid);
end

end


function out = my_conv(p, q)
    m = length(p);
    n = length(q);
    out = zeros(1, m + n - 1);
    for i = 1:m
        for j = 1:n
            out(i+j-1) = out(i+j-1) + p(i) * q(j);
        end
    end
end

