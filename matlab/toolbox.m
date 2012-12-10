function tube = toolbox(l1, l2, k0, k1, N_directions)

A_tool = A(0);
B_tool = B(0);
Q_tool = Q(0);
q_tool = Q_c(0);
m_tool = X_0_c();
M_tool = X_0();

% A_tool = {'sin(k)' '0' '0'; '1/2*cos(k)' '0' '1/2*cos(k)'; '0' 'sin(k)' '0'}

dim = size(A_tool, 1);

e1 = l1 / norm(l1);
e2 = l2 - e1 * (e1' * l2);
e2 = e2 / norm(e2);
Ell = cell(1, N_directions);
for i = 1:N_directions
  Ell{i} = rand(1, dim);
  Ell{i} = Ell{i} / norm(Ell{i});
end
Ell = cell2mat(Ell');
ls = linsys(A_tool, B_tool, ellipsoid(q_tool, Q_tool), [], [], [], [], 'd');
r = reach(ls, ellipsoid(X_0_c, X_0), Ell', [k0 k1]);
tube = projection(r,[e1';e2']');

end