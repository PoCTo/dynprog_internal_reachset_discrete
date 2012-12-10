function [reach_tube, reach_ellipses] = ...
                no_good_curves(l_1, l_2, k_0, k_1, directions_count)
some_A_value = A(k_0);
n = size(some_A_value, 1);

% assuming k_1 >= k_0
t_count = k_1 - k_0 + 1;

As = zeros(n, n, t_count);
As(:,:,1) = some_A_value;
for (i = (k_0+1):k_1)
    %[i - k_0 + 1]
    As(:,:,i - k_0 + 1) = A(i);
end
    
Fundamentals = zeros(n,n,t_count,t_count);
for (i = 1:t_count)
    Fundamentals(:,:,i,i) = eye(n,n);
    %[i,i]
end
for (i = k_0:k_1)
    for (j = (i+1):k_1)
        Fundamentals(:,:,i-k_0+1,j-k_0+1) = ...
            As(:,:,j-k_0) * Fundamentals(:,:,i-k_0+1,j-k_0);
        %[i-k_0+1,j-k_0+1]
    end
end

% Normalized random directions
directions = rand(n, directions_count);
norms = arrayfun(@(i)norm(directions(:,i)),1:directions_count);
directions = directions ./ repmat(norms,n,1);

centers = zeros(n, t_count);
centers(:,1) = X_0_c;

for (k=(k_0 + 1):k_1)
    i = k - k_0 + 1;
    centers(:,i) = Fundamentals(1, i)*X_0_c;
    for (j = k_0:k-1)
        %Fundamentals(:,:,j+1-k_0+1,i)
        % B(j-k_0+1)
        % Q_c(j-k_0+1)
        centers(:,i) = centers(:,i) + Fundamentals(:,:,j+1-k_0+1,i)*...
            B(j-k_0+1)*Q_c(j-k_0+1);
    end
end

ellipses = zeros(n, n, t_count, directions_count);
for (i=1:directions_count)
    ellipses(:,:,1,i) = X_0;
end

approximations = cell(t_count);
approximations{1} = ellipsoidalProjection(X_0_c,X_0,l_1,l_2,200);

for (k=(k_0+1):k_1) % t_count
    border_set = cell(1, directions_count);
    i = k - k_0 + 1;
    S = eye(n,n);
    for (j=1:directions_count)
        ellipses(:,:,i,j) = S*sqrtm(X_0);
        for (k_star=k_0:(k-1))
            i_star = k_star-k_0+1; % i_star
            l = directions(:,j);
            l_k_star = sqrt(dot(l,...
                Fundamentals(:,:,1,i_star)*B(k_star)*Q(k_star)*...
                    B(k_star)'*Fundamentals(:,:,1,i_star)'*directions(:,j))/...
                        dot(l,Fundamentals(:,:,1,i_star)*X_0*Fundamentals(:,:,1,i_star)'*l));
            ellipses(:,:,i,j) = ellipses(:,:,i,j) + ...
                l_k_star*S*sqrtm(X_0);
        end
        ellipses(:,:,i,j) = ellipses(:,:,i,j)*Fundamentals(:,:,1,i)';
        ellipses(:,:,i,j) = ellipses(:,:,i,j)'*ellipses(:,:,i,j);
        ellipses(:,:,i,j)
        border_set{1,j} = ...
            ellipsoidalProjection(centers(:,i),ellipses(:,:,i,j),...
                l_1, l_2, 200);
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
    approximations{i} = [xs; ys];
    %[xs;ys]
end

figure();
% colormap copper
% shading interp
hold on;
for k = k_0:k_1
  k_sh = k-k_0+1;
  len = size(approximations{k_sh}, 2);

  mesh(linspace(k-1, k, len)'*ones(1, len), ...
       ones(1, len)'*approximations{k_sh}(1, :),...
       ones(1, len)'*approximations{k_sh}(2, :));
  grid on
end
xlabel('$t$','interpreter', 'latex');
ylabel('$e_1$','interpreter', 'latex');
zlabel('$e_2$','interpreter', 'latex');







