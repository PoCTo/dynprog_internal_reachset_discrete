function [approximations, centers, ellipses] = ...
                no_good_curves(l_1, l_2, k_0, k_1, directions_count)

syms k;
%A_s(k)=sym(A);
%B_s(k)=sym(B);
%Q_s(k)=sym(Q);
%Q_c_s(k)=sym(Q_c);

h = waitbar(0, 'Calculating...');

plot_count = 30;

A_sym = @(t)subs(A,k,t);
B_sym = @(t)subs(B,k,t);
Q_sym = @(t)subs(Q,k,t);
Q_c_sym = @(t)subs(Q_c);

to_array=@(k)(k-k_0+1);
            
some_A_value = A_sym(k_0);
n = size(some_A_value, 1);

% assuming k_1 >= k_0
t_count = k_1 - k_0 + 1;

steps = t_count*directions_count+3;
step = 0;

As = cell(t_count,1);
As{1} = some_A_value;
for (i = (k_0+1):k_1)
    %[i - k_0 + 1]
    As{to_array(i)} = A_sym(i);
end
    
Fundamentals = cell(t_count,t_count);
for (i = 1:t_count)
    Fundamentals{i,i} = eye(n,n);
    %[i,i]
end
for (i = k_0:k_1)
    for (j = (i+1):k_1)
        Fundamentals{to_array(i),to_array(j)} = ...
            As{to_array(j-1)} * Fundamentals{to_array(i),to_array(j-1)};
        %[i-k_0+1,j-k_0+1]
    end
end

step = step+1;
waitbar(step/steps,h);


% Normalized random directions
directions = cell(directions_count);
for i=1:directions_count
    directions{i} = rand(n,1)*2 - 1;
    directions{i} = directions{i}/norm(directions{i});
end

centers = cell(t_count);
centers{1} = X_0_c;

fund = @(x,y)Fundamentals{to_array(y),to_array(x)};

for (k=(k_0 + 1):k_1)
    centers{to_array(k)} = A_sym(k-1)*centers{to_array(k-1)}+...
                           B_sym(k-1)*Q_c_sym(k-1);    
end


step = step+1;
waitbar(step/steps,h);

ellipses = cell(t_count, directions_count);
for (i=1:directions_count)
    ellipses{1,i} = X_0;    
    
    step = step+1;
    waitbar(step/steps,h);
end

approximations = cell(t_count,1);
approximations{1} = ellipsoidalProjection(X_0_c,X_0,l_1,l_2,100);

for (k=(k_0+1):k_1) % t_count
    border_set = cell(1, directions_count);
    S = eye(n,n);
    for (j=1:directions_count)
        ellipses{to_array(k),j} = S*sqrtm(X_0+eye(n,n)*1e-6)*(fund(k,k_0)');
        l = directions{j};
        for (i=k_0:(k-1))
            ellipses{to_array(k),j} = ellipses{to_array(k),j}+...
                generate_rotation_matrix(S*sqrtm(X_0+eye(n,n)*1e-6)*(fund(k,k_0)')*l,...
                    sqrtm((B_sym(i))*Q_sym(i)*(B_sym(i)')+eye(n,n)*1e-6)*(fund(k,i+1)')*l)*...
                sqrtm((B_sym(i))*Q_sym(i)*(B_sym(i)')+eye(n,n)*1e-6)*((fund(k,i+1))');
%             lambda_i = sqrt(dot(l,... 
%                 fund(k,i+1)*B_sym(i)*Q_sym(i)*...
%                     (B_sym(i)')*(fund(k,i+1)')*l));
%             lambda_i = lambda_i / sqrt(dot(l,fund(k,k_0)*X_0*(fund(k,k_0)')*l));
%             ellipses{to_array(k),j} = ellipses{to_array(k),j} + ...
%                 lambda_i;
        end
        %ellipses{to_array(k),j} = ellipses{to_array(k),j}*S*sqrtm(X_0)*(fund(k,k_0)');
        ellipses{to_array(k),j} = (ellipses{to_array(k),j}')*ellipses{to_array(k),j};
        %ellipses{to_array(k),j};
        
        border_set{j} = ...
            ellipsoidalProjection(centers{to_array(k)},ellipses{to_array(k),j},...
                l_1, l_2, plot_count);
            
        step = step+1;
        waitbar(step/steps,h);
    end
    xs = border_set{1}(1, :);
    ys = border_set{1}(2, :);
    [xs,ys] = poly2cw(xs,ys);
    %xs = xs(end:-1:1);
    %ys = ys(end:-1:1);
    for r = 2:directions_count
        [border_set{r}(1,:),border_set{r}(2,:)] = ...
            poly2cw(border_set{r}(1,:),border_set{r}(2,:));
        [xs,ys] = polybool('union',xs,ys,...
            border_set{r}(1,:),border_set{r}(2,:));
    end
    approximations{to_array(k)} = [xs; ys];
    %[xs;ys]
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









