function regQMat = regularize_symm(qMat,absTol)
%
% REGULARIZE - regularization of singular symmetric matrix.
%
  [~, n] = size(qMat);
  r      = rank(qMat);
  if r < n
    [U, ~, ~] = svd(qMat);
    E       = absTol * eye(n - r);
    regQMat       = qMat + (U * [zeros(r, r) zeros(r, (n-r)); zeros((n-r), r) E] * U');
    regQMat       = 0.5*(regQMat + regQMat');
  else
    regQMat = qMat;
  end
