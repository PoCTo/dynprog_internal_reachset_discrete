function [approximations, centers, ellipses, ls] = ...
                good_curves(l_1, l_2, k_0, k_1, directions_count)
syms k;
%A_s(k)=sym(A);
%B_s(k)=sym(B);
%Q_s(k)=sym(Q);
%Q_c_s(k)=sym(Q_c);

h = waitbar(0, 'Calculating...');

A_sym = @(t)subs(A,k,t);
B_sym = @(t)subs(B,k,t);
Q_sym = @(t)subs(Q,k,t);
Q_c_sym = @(t)subs(Q_c);

to_array=@(k)(k-k_0+1);
            
some_A_value = A_sym(k_0);
n = size(some_A_value, 1);

% assuming k_1 >= k_0
t_count = k_1 - k_0 + 1;

steps = t_count*directions_count;
step = 0;

As = cell(t_count,1);
As{1} = some_A_value;
for (i = (k_0+1):k_1)
    %[i - k_0 + 1]
    As{to_array(i)} = A_sym(i);
end
    
Fundamentals = cell(t_count); % fund(k_0,k)
Fundamentals{to_array(k_0)} = eye(n,n);
for i=k_0+1:k_1
    Fundamentals{to_array(i)} = A_sym(i-1)*Fundamentals{to_array(i-1)};
end

% Normalized random directions
directions = cell(directions_count);
for i=1:directions_count
    directions{i} = rand(n,1);
    directions{i} = directions{i}/norm(directions{i});
end
directions{1} = l_1;
directions{2} = l_2;

centers = cell(t_count);
centers{1} = X_0_c;

fund = @(y)Fundamentals{to_array(y)};


% Centers
for (k=(k_0 + 1):k_1)
    centers{to_array(k)} = A_sym(k-1)*centers{to_array(k-1)}+...
                           B_sym(k-1)*Q_c_sym(k-1);    
end

ellipses = cell(t_count, directions_count);
ls = cell(t_count);
for i=1:t_count
    ls{i} = directions;
end
for (i=1:directions_count)
    ellipses{1,i} = X_0;
end

approximations = cell(t_count,1);
approximations{1} = ellipsoidalProjection(X_0_c,X_0,l_1,l_2,100);

M = eye(n,n);
S = eye(n,n);
x0s = sqrtm(X_0+eye(n,n)*1e-6)*S';
border_set = cell(t_count, directions_count);

for dir=1:directions_count
    for k=(k_0+1):k_1
        ls{to_array(k)}{dir} = fund(k)*directions{dir};
    end
end
step = step+t_count;
waitbar(step/steps,h);


for dir=1:directions_count
    l = directions{dir};
    M_k = sqrtm(X_0+eye(n,n)*1e-6);
    S = eye(n,n);
    S_k = S;
    
    for k=(k_0):(k_1-1)
        
        ellipses{to_array(k+1),dir} = A_sym(k)*ellipses{to_array(k),dir}*(A_sym(k)') + ...
            ((B_sym(k))*Q_sym(to_array(k))*(B_sym(k)'))...
               + ...
            fund(k+1)*(M_k*S_k*sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*(pinv(fund(k+1))')+...
            sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*(pinv(fund(k+1))')*(S_k')*(M_k'))*...
            (fund(k+1)');
        
        M_k = M_k + sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*(fund(k+1)')*S_k;
        S_k = generate_rotation_matrix(S*sqrtm(X_0+eye(n,n)*1e-6)*l,...
            sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*pinv(fund(k+1)')*l);
      
        
        l1 = ls{to_array(k)+1}{1}/norm(ls{to_array(k)+1}{1});
        l2 = ls{to_array(k)+1}{2}/norm(ls{to_array(k)+1}{2});
        border_set{to_array(k)+1,dir} = ...
            ellipsoidalProjection(centers{to_array(k)+1},ellipses{to_array(k)+1,dir},...
                l1, l2, 100);
        step = step+1;
        waitbar(step/steps,h);

    end
end
for (k=(k_0+1):k_1)
    xs = border_set{to_array(k),1}(1, :);
    ys = border_set{to_array(k),1}(2, :);
    [xs,ys] = poly2cw(xs,ys);
    %xs = xs(end:-1:1);
    %ys = ys(end:-1:1);
    for r = 2:directions_count
        [border_set{to_array(k),r}(1,:),border_set{to_array(k),r}(2,:)] = ...
            poly2cw(border_set{to_array(k),r}(1,:),border_set{to_array(k),r}(2,:));
        [xs,ys] = polybool('union',xs,ys,...
            border_set{to_array(k),r}(1,:),border_set{to_array(k),r}(2,:));
    end
    approximations{to_array(k)} = [xs; ys];
end

close(h);

% figure();
% % colormap copper
% % shading interp
% hold on;
% for k = k_0:k_1
%   k_sh = k-k_0+1;
%   len = size(approximations{k_sh}, 2);
%   leftk=k-1/2;
%   rightk=k+1/2;
%   if (k==k_0)
%       leftk=k;
%   elseif (k==k_1)
%       rightk=k;
%   end
%      [leftk,rightk] 
%   mesh(linspace(leftk, rightk, len)'*ones(1, len), ...
%        ones(1, len)'*approximations{k_sh}(1, :),...
%        ones(1, len)'*approximations{k_sh}(2, :));
%   grid on
% end
% xlabel('$t$','interpreter', 'latex');
% ylabel('$e_1$','interpreter', 'latex');
% zlabel('$e_2$','interpreter', 'latex');
