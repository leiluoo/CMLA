%§¡§Ó§ä§à§â§ã§Ü§à§Ö §á§â§Ñ§Ó§à ? 2022 §£§®§¬-3 §°§å§ñ§ß §­§ï§Û§Ý§à
function [x,t] = lu_guass(oyll)
%[Description]
%This function is for labotory task 1(b), and it is to slove linear
%system Ax = b. In the 1st step and 2cd step I read matirx A and vector b
%from given files. In the 3rd step I ought to find out if columns of A
%match elements of b. In the 4th and 5th step, I try to get existing matrix
%L and U. In the step 6, I solve the equation Ly = b and Ux = y and then
%find out vector x. In the last step, I write the result vector x in a new
%file.
%[Parameters]
%oyll - This is the '*', which was previously introduced in task-detail
%file
%[Returns]
%x - x vector
%t - t run time


%document cputime t1
t1 = cputime;

%********task 1 (b)********
%********step 1********
file_name_A = sprintf("data\\Amat%d.m", oyll);
run(file_name_A) %run the resource file and load data by the way

%********step 2********
file_name_b = sprintf("data\\bvec%d.m", oyll);
run(file_name_b) %run the resource file and load data by the way

%********step 3********
[column_A, ~] = size(A); %get size of matrix A
[column_b, ~] = size(b); %get size of matrix b

%judge if column_A equal column_b
if column_A ~= column_b
    error("column_A and column_b are unequal");
end

%********step 4&step 5********
try 
    %if it can be found in existing file
    file_name_L = sprintf("data\\Lmat%d.m", oyll);
    file_name_U = sprintf("data\\Umat%d.m", oyll);
    run(file_name_L);
    run(file_name_U);
catch
    %if it can not
    file_name_L = sprintf("result\\Lmat%d.m", oyll);
    file_name_U = sprintf("result\\Umat%d.m", oyll);
    make_lu(oyll);
    run(file_name_L);
    run(file_name_U);
end

%********step 6********
n = column_b;
y = zeros(n,1);
y(1) = b(1);

%solve Ly = b
for i = 2:n
    for j = 1:i-1
        b(i) = b(i) - L(i,j) * y(j);
    end
    y(i) = b(i);
end

%solve Ux = y
x(n) = y(n) / U(n,n);
for i = n-1:-1:1
    for j = n:-1:i+1
        y(i) = y(i) - U(i,j) * x(j);
    end
    x(i) = y(i) / U(i,i);
end

%********step 7********
x_file = sprintf('result\\xvec%d.m', oyll);
fid = fopen(x_file, 'w');
%write x by following steps
if isreal(x)
    %if it is a real matrix
    fprintf(fid,['x = [',repmat('%8.4f;',1,column_b),'];'], x');
else
    %if it is not
    fprintf(fid,['x = complex([',repmat('%8.4f;',1,column_b),'],'], real(x)');
    fprintf(fid,['[',repmat('%8.4f;',1,column_b),']);'], imag(x)');
end
fclose(fid);

%document cputime t2
t2 = cputime;
%compute running time t
t = t2 - t1;

%%********task 1 (c)********
if oyll == 4
    task1c_file = 'result\\task_1c.txt';
    fid = fopen(task1c_file, 'w');
    fprintf(fid, "run time = %f s \n\n", t);
    %write x by following steps
    if isreal(x)
        %if it is a real matrix
        fprintf(fid,['answer:\nx = [',repmat('%8.4f;',1,column_b),'];'], x');
    else
        %if it is not
        fprintf(fid,['answer:\nx = complex([',repmat('%8.4f;',1,column_b),'],'], real(x)');
        fprintf(fid,['[',repmat('%8.4f;',1,column_b),']);'], imag(x)');
    end
    fclose(fid);
end
end

