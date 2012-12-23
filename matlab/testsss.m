% A = @(t) [-1 1 0 1; 1 1 1 0; 1 0 -1 0; -1 1 1 -1];
% B = @(t) [0.1 0 0.2; 0.2 0 0.1; 0.1 0 -0.2; -0.2 0 0.1];
% Q = @(t) [1 0 1; 0 1 1; 1 1 3];
% q = @(t) [0 1 1]';
% m = [0 0 0 0]';
% M = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
% l1 = [0 3 1 0]';
% l2 = [3 0 2 1]';
% k1 = 5;
% k0 = -2;

% A = @(t) [0.8 0 3 -2; 0 0.8 1 0; 0 0 2 0; 0 0 0 2];
% A = @(t) [0.5 0.5 0 0.5; 1 0.5 0.5 0; 0 1 1 1; 1 1 1 1];
% B = @(t) [1 0 1; 0.5 0.5 1; 0 0 1; 0 1 0];
% Q = @(t) [0.5 1 0; 1 4 0; 0 0 16];
% q = @(t) [1 1 1]';
% m = [0.5 0.5 0.5 0.5]';
% M = [1 1 0 0; 1 4 0 0; 0 0 3 2; 0 0 2 3];
% l1 = [1 3 2 0]';
% l2 = [0 1 1 4]';
% k1 = 6;
% k0 = 1;

k1 = 5;
N_dim = [5 3 5 10 20 50 100 150 200 250 300];

time2 = [];

% N_dim = [3:10:500];
% N_dim = [5 5 10 20 50 100 150 200 250];
% % N_dim = [20 50 250]
% 
% for n_dim = N_dim(1:end)
%   n_dim
%   tic
% 
% 
% % A = @(t) (rand(n_dim) + ones(n_dim)) * matrix;
% A =  eye(n_dim);
% B =  eye(n_dim);
% Q =  eye(n_dim);
% q =  zeros(1, n_dim)';
% m = zeros(1, n_dim)';
% M = eye(n_dim);
% l1 = zeros(1, n_dim)';
% l1(1) = 1;
% l2 = zeros(1, n_dim)';
% l2(2) = 1;
% k0 = 1;
% 
% A_tool = mat2cell(A,linspace(1,1,n_dim),linspace(1,1,n_dim));
% B_tool = mat2cell(B,linspace(1,1,n_dim),linspace(1,1,n_dim));
% Q_tool = mat2cell(Q,linspace(1,1,n_dim),linspace(1,1,n_dim));
% q_tool = mat2cell(q,linspace(1,1,n_dim));
% 
% for i=1:n_dim
%     for j=1:n_dim
%         A_tool{i,j}=num2str(A_tool{i,j});
%         B_tool{i,j}=num2str(B_tool{i,j});
%         Q_tool{i,j}=num2str(Q_tool{i,j});
%         q_tool{i}=num2str(q_tool{i});
%     end
% end
% m_tool = m;
% M_tool = M;
% 
% % dim = size(A(0), 1);
% % N_directions = 1;
% % Ell = cell(1, N_directions);
% % for i = 1:N_directions
% %   Ell{i} = rand(1, dim)';
% %   Ell{i} = Ell{i} / norm(Ell{i});
% % end
% 
% % A_tool = {'sin(k)' '0' '0'; '1/2*cos(k)' '0' '1/2*cos(k)'; '0' 'sin(k)' '0'}
% 
% dim = n_dim;%size(A_tool, 1);
% 
% % e1 = l1 / norm(l1);
% % e2 = l2 - e1 * (e1' * l2);
% % e2 = e2 / norm(e2);
% % Ell = cell(1, N_directions);
% % for i = 1:N_directions
% %   Ell{i} = rand(1, dim);
% %   Ell{i} = Ell{i} / norm(Ell{i});
% % end
% no_good_curves2(A_tool,B_tool,Q_tool,q_tool,M_tool,m_tool,l1,l2,k0,k1,1);
% 
% % plot_ea(tube);
% 
%   time2 = [time2 toc]
% 
% end

time1 = [];

% N_dim = [3:10:500];

% N_dim = [20 50 250]

for n_dim = N_dim(1:end)
  n_dim
  tic


% A = @(t) (rand(n_dim) + ones(n_dim)) * matrix;
A = @(t) eye(n_dim);
B = @(t) eye(n_dim);
Q = @(t) eye(n_dim);
q = @(t) zeros(1, n_dim)';
m = zeros(1, n_dim)';
M = eye(n_dim);
l1 = zeros(1, n_dim)';
l1(1) = 1;
l2 = zeros(1, n_dim)';
l2(2) = 1;

k0 = 1;

A_tool = A(0);
B_tool = B(0);
Q_tool = Q(0);
q_tool = q(0);
m_tool = m;
M_tool = M;

dim = size(A(0), 1);
N_directions = 1;
Ell = cell(1, N_directions);
for i = 1:N_directions
  Ell{i} = rand(1, dim)';
  Ell{i} = Ell{i} / norm(Ell{i});
end

% A_tool = {'sin(k)' '0' '0'; '1/2*cos(k)' '0' '1/2*cos(k)'; '0' 'sin(k)' '0'}

dim = size(A_tool, 1);

e1 = l1 / norm(l1);
e2 = l2 - e1 * (e1' * l2);
e2 = e2 / norm(e2);
% Ell = cell(1, N_directions);
% for i = 1:N_directions
%   Ell{i} = rand(1, dim);
%   Ell{i} = Ell{i} / norm(Ell{i});
% end
Ell = cell2mat(Ell);
ls = linsys(A_tool, B_tool, ellipsoid(q_tool, Q_tool), [], [], [], [], 'd');
r = reach(ls, ellipsoid(m, M), Ell, [k0 k1]);
tube = projection(r, [e1';e2']');
ea = get_ea(cut(r, k0));
ea_pr = projection(ea, [e1';e2']');

% plot_ea(tube);

  time1 = [time1 toc]
end
% figure();
% % colormap copper
% % shading interp
% hold on;
% for k = k0:k1
%   k_sh = k-k0+1;
%   len = size(Solv_set{k_sh}, 2);
% 
%   mesh(linspace(k-1, k, len)'*ones(1, len), ...
%        ones(1, len)'*Solv_set{k_sh}(1, :), ones(1, len)'*Solv_set{k_sh}(2, :));
%   grid on
% end
% xlabel('$t$','interpreter', 'latex');
% ylabel('$e_1$','interpreter', 'latex');
% zlabel('$e_2$','interpreter', 'latex');



[N_dim; time1]