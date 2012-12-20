function Border = ellipsoidalProjection(q, Q, l1, l2, N)
%(q, Q, l1, l2, count)
%   e1 = l1 / norm(l1);
%   e2 = l2 - e1 * (e1' * l2);
%   e2 = e2 / norm(e2);  
%   angles = linspace(0,2*pi,count);  
%   points = e1*cos(angles) + e2*sin(angles);
%   points = mat2cell(points, length(l1), ones(1,count));
% 
%   find_support_point = @(ell) q + Q * ell / sqrt(ell' * Q * ell);
%   project = @(ell) [ell' * e1; ell' * e2];
%   coords = cellfun(find_support_point, points', 'UniformOutput', 0);
%   coords = cellfun(project, coords', 'UniformOutput', 0);
%   coords = cell2mat(coords);
% end
  e1 = l1 / norm(l1);
  e2 = l2 - e1 * (e1' * l2);
  e2 = e2 / norm(e2);
  
  Alpha = linspace(0, 2*pi, N);
  
  Ell = e1 * cos(Alpha) + e2 * sin(Alpha);
  Ell = mat2cell(Ell, length(l1), ones(1, N));

  find_support_point = @(ell) q + Q * ell / sqrt(ell' * Q * ell);
  project = @(ell) [ell' * e1; ell' * e2];
  Border = cellfun(find_support_point, Ell', 'UniformOutput', 0);
  Border = cellfun(project, Border', 'UniformOutput', 0);
  Border = cell2mat(Border);
