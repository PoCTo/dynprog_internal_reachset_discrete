function tube = toolbox(l_1, l_2, k_0, k_1, N_directions)

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

e_1 = l_1 / norm(l_1);
e_2 = l_2 - e_1 * (e_1' * l_2);
e_2 = e_2 / norm(e_2);
Ell = cell(1, N_directions);
for i = 1:N_directions
  Ell{i} = rand(1, dim);
  Ell{i} = Ell{i} / norm(Ell{i});
end
Ell = cell2mat(Ell');
ls = linsys(A, B, Q_str, [], [], [], [], 'd');
r = reach(ls, ellipsoid(X_0_c, X_0), Ell', [k_0 k_1]);
tube = projection(r,[e_1';e_2']');

end