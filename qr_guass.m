%§¡§Ó§ä§à§â§ã§Ü§à§Ö §á§â§Ñ§Ó§à ? 2022 §£§®§¬-3 §°§å§ñ§ß §­§ï§Û§Ý§à
function [x,t] = qr_guass(oyll)
%[Description]
%This function is for labotory task 2(b), and it is to slove linear
%system Ax = b. In the 1st step and 2cd step I read matirx A and vector b
%from given files. In the 3rd step I ought to find out if columns of A
%match elements of b, meanwhile check if matirx A is square and non-singular.
%In the 4th and 5th step, I try to get existing matrix Q and R. In the step
%6, I solve the equation Qy = b and Rx = y and then find out vector x. In 
%the last step, I write the result vector x in a new file.
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
[column_A, row_A] = size(A); %get size of matrix A
[column_b, ~] = size(b); %get size of matrix b

%judge if column_A equal column_b
if column_A ~= column_b
    error("column_A and column_b are unequal");
end

%judge if matrix A is square and non-sigular
if row_A ~= column_A
    error("Matrix A is not square");
end

if det(A) == 0
    error("Matrix A is singular");
end

%********step 4&step 5********
try 
    %if it can be found in existing file
    file_name_Q = sprintf("data\\Qmat%d.m", oyll);
    file_name_R = sprintf("data\\Rmat%d.m", oyll);
    run(file_name_Q);
    run(file_name_R);
catch
    %if it can not
    file_name_Q = sprintf("result\\Qmat%d.m", oyll);
    file_name_R = sprintf("result\\Rmat%d.m", oyll);
    make_qr(oyll);
    run(file_name_Q);
    run(file_name_R);
end

%********step 6********0
n = column_b;
x = zeros(n,1);
%solve Qy = b
y = Q' * b;

%solve Rx = y
for i = n:-1:1
    x(i) = (y(i) - R(i, i+1:n) * x(i+1:n))/ R(i,i);
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
if oyll == 9
    task1c_file = 'result\\task_2c.txt';
    fid = fopen(task1c_file, 'w');
    fprintf(fid, "run time = %f s \n\n", t);
    %write x by following steps
    if isreal(x)
        %if it is a real matrix
        fprintf(fid,['answer:\nx = [',repmat('%10.4f;',1,column_b),'];'], x');
    else
        %if it is not
        fprintf(fid,['answer:\nx = complex([',repmat('%10.4f;',1,column_b),'],'], real(x)');
        fprintf(fid,['[',repmat('%10.4f;',1,column_b),']);'], imag(x)');
    end
    fclose(fid);
end
end


