function C = new_vertex_paths(S,vertices)
% S: matriz S del modelo
% vertices: posicion (index) de los vertices
ban_v = [];
RB = [0];
NCV = 0;
H = [];
vertices = vertices';
for i=vertices
    disp(['Revisando nodo:              ', num2str(i)])
    pos_v = vertices;
    pos_v(pos_v==i) = []; % quitamos el vertice a analizar de pos_v
    [rxns_baneadas, ncv] = recorrer(S,pos_v,ban_v,i,i,RB,NCV,H);
    ban_v(end+1) = i; % se agrega el vertice utilizado
    RB = rxns_baneadas;
    NCV = ncv;
    disp(['Numero de caminos validos:   ', num2str(NCV)])
end

C = NCV;


function [rxns_baneadas, ncv] = recorrer(S,pos_v,ban_v,nodo,VI,RB,NCV,historial)
%% Definiciones
% S: matriz del modelo
% pos_v: posicion (index) de los vertices MENOS POR EL QUE SE COMIENZA
% ban_v: lista de vertices ya analizados
% nodo: Metabolito actual en el que estamos parados
% VI: Vertice inicial
% RB: reacciones baneadas (evita sobre-contar las rxns-multiples)
% NCV: numero de caminos validos
% historial: camino recorrido hasta el momento
%% Fin del camino
% disp('Estamos en el nodo/met:')
% disp(nodo)
if isequal(nodo,VI) && ~ isempty(historial) % Loop simple sin sentido
    rxns_baneadas = RB;
    ncv = NCV;
%     disp('Se encontro un loop simple')
    return
end
if any(ban_v == nodo) % si llegamos a un vertice ya recorrido
    rxns_baneadas = RB;
    ncv = NCV;
%     disp('llegamos a un vertice ya recorrido')
    return
end
if ~ isempty(find(pos_v==nodo)) % si es vertice (diff de VI)
    rxns_baneadas = RB;
    ncv = NCV + 1;
%     disp('encontramos un vertice')
    return
end
caminos = find(S(nodo,:)~=0); % rxns en las q participa el met
caminos_nuevos = setdiff(caminos,RB); % se quitan las rxns por las que ya pasamos -> todas
if isempty(caminos_nuevos) % Llegamos a un met que se le revisaron todos los bordes
    rxns_baneadas = RB;
    ncv = NCV;
%     disp('Llegamos a un exchange')
    return
end
% disp('los caminos nuevos disponibles son:')
% disp(caminos_nuevos)
%% Recorrer los caminos
for i = caminos_nuevos
%     disp(i)
    if isempty(find(RB==i)) % si i no está en las rxns_baneadas
        mets = find(S(:,i)); % mets de la reaccion
        mets(mets==nodo) = []; % quitamos el met original
%         if size(mets,1) > 1 % se agrega rxn baneada
%             RB(end+1) = i;
%         end
        RB(end+1) = i; % Baneamos todas las rxns pq no podemos pasar 2 veces x ninguna
        historial_nuevo = historial;
        historial_nuevo(end+1) = i;
        for j = mets'
            [new_rxns_baneadas, new_ncv] = recorrer(S,pos_v,ban_v,j,VI,RB,NCV,historial_nuevo);
%             disp('salimos de la funcion con RB:')
%             disp(new_rxns_baneadas)
%             disp('y con NCV:')
%             disp(new_ncv)
            RB = new_rxns_baneadas;
            NCV = new_ncv;
        end

    end
end
rxns_baneadas = RB;
ncv = NCV;
return



