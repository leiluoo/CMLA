
function [Q,R] = make_qr(oyll)
%[Description]
%This function is for labotory task 2(a), and it is to perform a QR
%decomposition of matrix A. In the 1st and 2cd step I ought to initially
%judge if matrix A is square, cause it can not be decomposed if it is not.
%In the 3rd step, I perform 3 methods of QR decomposition
%In the last step, I write the result matrices L and U in a new  file.
%[Parameters]
%oyll - This is the '*', which was previously introduced in task-detail
%file
%[Return]
%Q - Q matrix
%R - R matrix

%********task 1 (a)********
%********step 1********
file_name = sprintf("data/Amat%d.m", oyll);
run(file_name) %run the resource file and load data by the way

%********step 2********
[m,n] = size(A);
%judge if it is square matrix
if m ~= n
    disp('The given matrix is not square'); 
end

n = size(A,1); %number of columns

%********step 3********
%--------Method 1 => The Householder QR method--------
if Method == 1
    R = A;
    Q = eye(n);
    I = eye(n);
    for j = 1:n-1
        x = R(j:n,j);
        v = -sign(x(1)) * norm(x) * eye(n-j+1,1) - x;
        if norm(v) > 0
            v = v / norm(v);
            P = I;
            P(j:n,j:n) = P(j:n,j:n) - 2 * v * v';
            R = P * R;
            Q = Q * P;
        end
    end
%-----------------------------------------------------
%----------Method 2 => The Givens QR method-----------
elseif Method == 2
    Q = eye(n);
    R = A;
    for j = 1:n
        for i = n:(-1):j+1
            x = R(:,j);
            if norm([x(i-1),x(i)]) > 0
                c = x(i-1) / norm([x(i-1),x(i)]);
                s = -x(i) / norm([x(i-1),x(i)]);
                G = eye(n); 
                G([i-1,i],[i-1,i]) = -[c,s;-s,c];
                R = G' * R;
                Q = Q * G;
            end
        end
    end
%-----------------------------------------------------
%--------Method 3 => The Gram-Schmidt process---------
% else
%     Q = zeros(n); 
%     R = zeros(n);
%     for j = 1:n
%         v = A(:,j);
%         for i = 1:j-1
%             R(i,j) = Q(:,i)' * v;
%             v = v - R(i,j) * Q(:,i);
%         end
%         R(j,j)= norm(v);
%         Q(:,j)= v / R(j,j);
%     end
%     
% end
else
    Q = zeros(m, n);
    R = zeros(n, n);
    max_iter = 20;
    for j = 1:n
        v = A(:,j);
        
        % Orthogonalize against previous vectors
        for i = 1:j-1
            R(i,j) = Q(:,i)' * v;
            v = v - R(i,j) * Q(:,i);
        end
        
        % Improve orthogonality by reorthogonalization
        for k = 1:max_iter
            for i = 1:j-1
                R(i,j) = R(i,j) + Q(:,i)' * v;
                v = v - (Q(:,i)' * v) * Q(:,i);
            end
        end
        
        % Normalize v to obtain q_j
        R(j,j) = norm(v);
        Q(:,j) = v / R(j,j);
    end
end
%-----------------------------------------------------

%********step 4********
%create a new file to save matrix L
Q_file = sprintf('result/Qmat%d.m', oyll);
fid = fopen(Q_file, 'w');

%write L by following steps
if isreal(Q)
    %if it is a real matrix
    fprintf(fid,'Q = ...\n');
    fprintf(fid,['[',repmat('%10.4f',1,n),';\n'], Q(1,:)');
    fprintf(fid,[repmat('%10.4f',1,n),';\n'], Q(2:n-1,:)');
    fprintf(fid,[repmat('%10.4f',1,n),'];\n'], Q(n,:)');
else
    %if it is not
    fprintf(fid,'L = ...\ncomplex(');
    fprintf(fid,['[',repmat('%10.4f',1,n),';\n'], real(Q(1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),';\n'], real(Q(2:n-1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),'],'], real(Q(n,:))');
    fprintf(fid,['[',repmat('%10.4f',1,n),';\n'], imag(Q(1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),';\n'], imag(Q(2:n-1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),']);\n'], imag(Q(n,:))');
end
fclose(fid);

%create a new file to save matrix U
R_file = sprintf('result/Rmat%d.m', oyll);
fid = fopen(R_file, 'w');
%write L by following steps
if isreal(R)
    %if it is a real matrix
    fprintf(fid,'R = ...\n');
    fprintf(fid,['[',repmat('%10.4f',1,n),';\n'], R(1,:)');
    fprintf(fid,[repmat('%10.4f',1,n),';\n'], R(2:n-1,:)');
    fprintf(fid,[repmat('%10.4f',1,n),'];\n'], R(n,:)');
else
    %if it is not
    fprintf(fid,'R = ...\ncomplex(');
    fprintf(fid,['[',repmat('%10.4f',1,n),';\n'], real(R(1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),';\n'], real(R(2:n-1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),'],'], real(R(n,:))');
    fprintf(fid,['[',repmat('%10.4f',1,n),';\n'], imag(R(1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),';\n'], imag(R(2:n-1,:))');
    fprintf(fid,[repmat('%10.4f',1,n),']);\n'], imag(R(n,:))');
end
fclose(fid);


end

