function [res] = projection(L,l1,l2)
e1 = l1 / norm(l1);
e2 = l2 - e1 * (e1' * l2);
e2 = e2 / norm(e2);

project = @(ell) [ell' * e1; ell' * e2];
%res = cellfun(project, L, 'UniformOutput', 0);
res = project(L);