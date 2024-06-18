function a_R = alpha_rank(model)
% alpha_rank function definition

S = full(model.S);

% Cantidad de vertices
vertices = zeros(size(S,1),1);
vval = sum(S~=0,2);
vertices(vval>2) = 1;
disp(['El modelo tiene ',num2str(sum(vertices)),' vertices'])

% Cantidad de caminos
vert_index = find(vertices~=0);
C = new_vertex_paths(S,vert_index);

% Rango de alpha
vval(vval<3) = [] ;
edges = unique(vval);
counts = zeros(length(edges),1);
pos = 1;
for i = edges'
    counts(pos) = length(find(vval==i));
    pos = pos + 1;
end
a_R = sum(counts .* (edges -1)) - C;
