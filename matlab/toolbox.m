function [tube, centers, ellipses] = toolbox(l_1, l_2, k_0, k_1, N_directions)

%A_tool = A(0);
%B_tool = B(0);
%Q_tool = Q(0);
%q_tool = Q_c(0);
%m_tool = X_0_c();
%M_tool = X_0();
syms k;
%A_tool(k)=sym(A);
%B_tool(k)=sym(B);
%Q_tool(k)=sym(Q);
%q_tool(k)=sym(Q_c);
Q_str = struct('center','','shape','');
Q_str.center = Q_c
Q_str.shape = Q

% A_tool = {'sin(k)' '0' '0'; '1/2*cos(k)' '0' '1/2*cos(k)'; '0' 'sin(k)' '0'}

dim = size(A, 1)

h = waitbar(0, 'Calculating...');

e_1 = l_1 / norm(l_1);
e_2 = l_2 - e_1 * (e_1' * l_2);
e_2 = e_2 / norm(e_2);
Ell = cell(1, N_directions);
for i = 1:N_directions
  Ell{i} = rand(1, dim);
  Ell{i} = Ell{i} / norm(Ell{i});
end
Ell = cell2mat(Ell');
waitbar(0.2,h,'Creating linsys...');
ls = linsys(A, B, Q_str, [], [], [], [], 'd');
waitbar(0.4,h,'Reaching...');
r = reach(ls, ellipsoid(X_0_c, X_0), Ell', [k_0 k_1]);
waitbar(0.7,h,'Projecting...');
tube = projection(r,[e_1';e_2']');

ea_pr = (get_ia(r))';
t_count = k_1 - k_0 + 1;
directions_count = N_directions;
ellipses = cell(t_count, directions_count);
centers = cell(t_count);
for (i=1:t_count)
    for (j=1:directions_count)
        [centers{i},ellipses{i,j}] = double(ea_pr(i,j));
    end
end


close(h);