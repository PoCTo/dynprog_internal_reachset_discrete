function R = regularize_nonsymm(Q, delta)
%
% ELL_REGULARIZE - regularization of singular matrix.
%

  [m, n] = size(Q);
  if m ~= n
    R = Q;
    return;
  end

  r = rank(Q);

  if r < n
    if min(min(Q == Q')) > 0
      R = Q  +  delta * eye(n);
    else
      [U S V] = svd(Q);
      R       = Q + (delta * U * V');
    end
  else
    R       = Q;
  end

  return;
