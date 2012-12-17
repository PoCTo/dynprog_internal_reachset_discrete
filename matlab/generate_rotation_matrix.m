function S = generate_rotation_matrix(from,to)
[U1,~,V1]=svd(from);
[U2,~,V2]=svd(to);
S=U1*V1*V2'*U2';