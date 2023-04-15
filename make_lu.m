%§¡§Ó§ä§à§â§ã§Ü§à§Ö §á§â§Ñ§Ó§à ? 2022 §£§®§¬-3 §°§å§ñ§ß §­§ï§Û§Ý§à
function [L, U] = make_lu(oyll)
%[Description]
%This function is for labotory task 1(a), and it is to perform a LU
%decomposition of matrix A. In the 1st and 2cd step I ought to initially
%judge if matrix A is square and regular, cause it can not be decomposed
%if it is not.In the 3rd step, I perform the formula of LU decomposition
%In the last step, I write the result matrices L and U in a new  file.
%[Parameters]
%oyll - This is the '*', which was previously introduced in task-detail
%file
%[Return]
%L - L matrix
%U - U matrix

%********task 1 (a)********
%********step 1********
file_name = sprintf("data\\Amat%d.m", oyll);
run(file_name) %run the resource file and load data by the way

%********step 2********
[a, b] = size(A);
%judge if it is square matrix
if a ~= b
    error('The given matrix is not square'); 
end

n = size(A,1); %number of columns

%judge if it is regular matrix
for i = 1:n
    if abs(det(A(1:i,1:i))) < 1e-300
        error('The given matrix is not regular');
    end
end

%********step 3********
Ak = A;
for k = 1:n
    % at this point AofK is A^(k-1)
    for j = k:n
        U(k,j) = Ak(k,j);
    end
    % check that we don't divide by zero...
    if U(k,k) < 1e-300
        error('** A^(k-1)_{k,k}==0 in LU decomp')
    end
    for i = k:n
        L(i,k) = Ak(i,k)/U(k,k);
    end
    % now modify AofK so that we can use it in the
    % next iteration
    for i = k:n
        for j = k:n
            Ak(i,j) = Ak(i,j) - L(i,k)*U(k,j);
        end
    end
end


%********step 4********
%create a new file to save matrix L
L_file = sprintf('result\\Lmat%d.m', oyll);
fid = fopen(L_file, 'w');

%write L by following steps
if isreal(L)
    %if it is a real matrix
    fprintf(fid,'L = ...\n');
    fprintf(fid,['[',repmat('%8.4f',1,n),';\n'], L(1,:)');
    fprintf(fid,[repmat('%8.4f',1,n),';\n'], L(2:n-1,:)');
    fprintf(fid,[repmat('%8.4f',1,n),'];\n'], L(n,:)');
else
    %if it is not
    fprintf(fid,'L = ...\ncomplex(');
    fprintf(fid,['[',repmat('%8.4f',1,n),';\n'], real(L(1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),';\n'], real(L(2:n-1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),'],'], real(L(n,:))');
    fprintf(fid,['[',repmat('%8.4f',1,n),';\n'], imag(L(1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),';\n'], imag(L(2:n-1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),']);\n'], imag(L(n,:))');
end
fclose(fid);

%create a new file to save matrix U
U_file = sprintf('result\\Umat%d.m', oyll);
fid = fopen(U_file, 'w');
%write L by following steps
if isreal(U)
    %if it is a real matrix
    fprintf(fid,'U = ...\n');
    fprintf(fid,['[',repmat('%8.4f',1,n),';\n'], U(1,:)');
    fprintf(fid,[repmat('%8.4f',1,n),';\n'], U(2:n-1,:)');
    fprintf(fid,[repmat('%8.4f',1,n),'];\n'], U(n,:)');
else
    %if it is not
    fprintf(fid,'U = ...\ncomplex(');
    fprintf(fid,['[',repmat('%8.4f',1,n),';\n'], real(U(1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),';\n'], real(U(2:n-1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),'],'], real(U(n,:))');
    fprintf(fid,['[',repmat('%8.4f',1,n),';\n'], imag(U(1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),';\n'], imag(U(2:n-1,:))');
    fprintf(fid,[repmat('%8.4f',1,n),']);\n'], imag(U(n,:))');
end
fclose(fid);

end

