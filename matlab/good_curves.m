function [approximations, centers, ellipses, points] = ...
                good_curves(l_1, l_2, k_0, k_1, directions_count)

%A_s(k)=sym(A);
%B_s(k)=sym(B);
%Q_s(k)=sym(Q);
%Q_c_s(k)=sym(Q_c);

h = waitbar(0, 'Calculating...');

syms k;
A_sym1 = @(t)subs(A,k,t);
%A_sym = @(t)(A_sym1(t) + A_svd(A_sym1(t))*const_svd);
B_sym = @(t)subs(B,k,t);
Q_sym = @(t)subs(Q,k,t);
Q_c_sym = @(t)subs(Q_c);

%const_svd = 1e-6;
eps = 1e-4;
sing = arrayfun(@(s)max(svd(A_sym1(s))),[k_0,k_1]);
si = max(sing);
ressing = 1;
for i = 1:k_1-k_0-1
    ressing = ressing + si^i;
end
const_sqrt = (eps / 2 / ressing)/2;

sing = arrayfun(@(s)max(svd(B_sym(s)*Q_sym(s)*B_sym(s)')),[k_0,k_1]);
si = max(sing);

syms ee1;
e1 = eval(solve(ee1+sqrt(ee1^2+2*si*ee1)-2*const_sqrt - 1e-5, ee1));

maxsinga = -100500;
for i=k_0:k_1
    Matx = A_sym1(i)*(A_svd(A_sym1(i))');
    Matx = Matx + Matx';
    svds = svd(Matx);
    maxsinga = max(max(svds),maxsinga);
end
const_svd = (-maxsinga+sqrt(maxsinga^2 + 4*eps^2))/3;
%A_sym = @(t)(A_sym1(t) + A_svd(A_sym1(t))*const_svd);

A_sym = @(t)regularize_nonsymm(A_sym1(t),1e-4);

R_sym = @(t)(regularize_symm(B_sym(t)*Q_sym(t)*(B_sym(t)'),1e-4));

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

directions_count = directions_count + 1;

% Normalized random directions
directions = cell(directions_count);
angles = linspace(0,2*pi,directions_count);
  e_1 = l_1 / norm(l_1);
  e_2 = l_2 - e_1 * (e_1' * l_2);
  e_2 = e_2 / norm(e_2);
  %e_1;
  %e_2;
for i=1:directions_count-1
    directions{i} = rand(n,1)*2 - 1;
    
    directions{i} = (e_1*cos(angles(i))+e_2*sin(angles(i)));
    %directions{i} = directions{i}/norm(directions{i});
    %directions{i};
end

directions{directions_count} = e_2;

% % Normalized random directions
% directions = cell(directions_count);
% for i=1:directions_count
%     directions{i} = rand(n,1);
%     directions{i} = directions{i}/norm(directions{i});
% end
% directions{1} = l_1;
% directions{2} = l_2;

centers = cell(t_count);
centers{1} = X_0_c;

fund = @(y)Fundamentals{to_array(y)};


% Centers
for (k=(k_0 + 1):k_1)
    centers{to_array(k)} = A_sym(k-1)*centers{to_array(k-1)}+...
                           B_sym(k-1)*Q_c_sym(k-1);    
end

ellipses = cell(t_count, directions_count);
% ls = cell(t_count);
% for i=1:t_count
%     ls{i} = directions;
% end
for (i=1:directions_count)
    ellipses{1,i} = X_0;
end

approximations = cell(t_count,1);
approximations{1} = ellipsoidalProjection(X_0_c,X_0,l_1,l_2,100);

%M = eye(n,n);
%S = eye(n,n);
%x0s = sqrtm(X_0+eye(n,n)*1e-6)*S';
border_set = cell(t_count, directions_count);

% for dir=1:directions_count
%     for k=(k_0+1):k_1
%         ls{to_array(k)}{dir} = fund(k)*directions{dir};
%     end
% end
step = step+t_count;
waitbar(step/steps,h);

points = zeros(2,directions_count, t_count);
find_support_point = @(q,Q,ell) q + Q * ell / sqrt(ell' * Q * ell);

for dir=1:directions_count
    l = directions{dir};
    M_k = sqrtm(X_0+eye(n,n)*const_sqrt);
    S = eye(n,n);
    %S_k = S;
    
    for p=(k_0):(k_1-1)
            ps = to_array(p);
        k = p;
        fkp1 = fund(k+1);
        fkp1r = pinv(fund(k+1));
        S_k = generate_rotation_matrix( ...
                S*sqrtm(X_0+eye(n,n)*const_sqrt)*l,...
...%                    sqrtm(fkp1r *... %fund(p+1)*...
...%                         (B_sym(p)*Q_sym(p)*(B_sym(p)')+eye(n,n)*const_sqrt)*... %(fund(p+1)')
...%                    (fkp1r'))*l ...
                ...%sqrtm((B_sym(p)*Q_sym(p)*(B_sym(p)')+eye(n,n)*const_sqrt))*...
                sqrtm(R_sym(p))*...
                (fkp1r')*l...
            );
        
        newsqrt = sqrtm(...
                        ...%fkp1r*...
                        ...%(B_sym(p)*Q_sym(p)*(B_sym(p)')+eye(n,n)*const_sqrt)...%*...
                        R_sym(p)...
                        ...%(fkp1r')...
                        )*(fkp1r');
        
        ellipses{ps+1,dir} =(...
            A_sym(k)*ellipses{ps,dir}*(A_sym(k)') + ...
            R_sym(p)+...%B_sym(p)*Q_sym(p)*(B_sym(p)') + ...
            fkp1*(...
                M_k * S_k * newsqrt...
                +... 
                newsqrt' * (S_k') * (M_k')...                        
            )*(fkp1'));        
        M_k = M_k + ...
            ...%sqrtm(...
                ...%fkp1r*(B_sym(p)*Q_sym(p)*(B_sym(p)')+eye(n,n)*const_sqrt)*(fkp1r'))...
                ...%sqrtm(B_sym(p)*Q_sym(p)*(B_sym(p)')+eye(n,n)*const_sqrt)*(fkp1r')...
                sqrtm(R_sym(p))*(fkp1r')...
                *(S_k');
        
%         ellipses{to_array(k+1),dir} = A_sym(k)*ellipses{to_array(k),dir}*(A_sym(k)') + ...
%             ((B_sym(k))*Q_sym(to_array(k))*(B_sym(k)'))...
%                + ...
%             fund(k+1)*(M_k*S_k*sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*(pinv(fund(k+1))')+...
%             sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*(pinv(fund(k+1))')*(S_k')*(M_k'))*...
%             (fund(k+1)');
%         
%         M_k = M_k + sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*(fund(k+1)')*S_k;
%         S_k = generate_rotation_matrix(S*sqrtm(X_0+eye(n,n)*1e-6)*l,...
%             sqrtm(B_sym(k)*Q_sym(k)*(B_sym(k)')+eye(n,n)*1e-6)*pinv(fund(k+1)')*l);
%       
        k = p;
        %l1 = fund(k+1)'*l_1;%/norm(fund(k+1)*l_1); %ls{to_array(k)+1}{1}/norm(ls{to_array(k)+1}{1});
        %l2 = fund(k+1)'*l_2;%/norm(fund(k+1)*l_2); %ls{to_array(k)+1}{directions_count}/norm(ls{to_array(k)+1}{directions_count});
        l1 = fkp1r'*l_1;
        l2 = fkp1r'*l_2;
        border_set{to_array(k)+1,dir} = real(...
            ellipsoidalProjection(centers{to_array(k)+1},ellipses{to_array(k)+1,dir},...
                l1, l2, 100));
        step = step+1;
        waitbar(step/steps,h);
        
        point = find_support_point(centers{to_array(p+1)},...
            ellipses{to_array(p+1),dir},fkp1r'*directions{dir});
        points(:,dir,ps+1) = projection2(point, l1, l2);
    end
end



% for i = 1:directions_count
%     point = find_support_point(centers{to_array(k_1)},...
%         ellipses{to_array(k_1),i},fund(k_1)*directions{i})
%     points(i,:) = projection(point, fund(k_1)*l_1, fund(k_1)*l_2);
% end
% points

for (k=(k_0+1):k_1)
    xs = border_set{to_array(k),1}(1, :);
    ys = border_set{to_array(k),1}(2, :);
    k
    [xs,ys] = poly2cw(xs,ys);
    
    %xs = xs(end:-1:1);
    %ys = ys(end:-1:1);
    for r = 2:directions_count
        k,r
        [border_set{to_array(k),r}(1,:),border_set{to_array(k),r}(2,:)] = ...
            poly2cw((border_set{to_array(k),r}(1,:)),(border_set{to_array(k),r}(2,:)));
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
