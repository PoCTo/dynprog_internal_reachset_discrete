function coords = ellipsoidalProjection(q,Q,l1,l2,count)
  e1 = l1 / norm(l1);
  e2 = l2 - e1 * (e1' * l2);
  e2 = e2 / norm(e2);  
  angles = linspace(-pi,pi,count);
  points = mat2cell( e1 * cos(angles) + e2 * sin(angles), length(l1), ones(1, count));

  find_support_point = @(points) q + Q * points / ...
                       sqrt(points' * Q * points);
  coords = (cell2mat((cellfun(find_support_point, points', 'UniformOutput', 0))'));
end