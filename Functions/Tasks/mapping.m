function [xf,yf,map] = mapping(xpositions,ypositions,parameter,xvert,yvert,imagesize)
imagesize = fliplr(imagesize);

%Buscar rectangulo mínimo que encierra al poligono adaptando sistema de coordenadas:
minx = round(min(xvert));
maxx = round(max(xvert));
maxy = round(imagesize(2) - min(yvert));
miny = round(imagesize(2) - max(yvert));

width = maxx - minx;
height = maxy - miny;

xpositions = xpositions - minx;
ypositions = imagesize(2) - ypositions;
ypositions = ypositions - miny;

% Agregar valores de bordes para evitar oscilaciones espurias:
xposbordertop = round(linspace(1,width,10));
yposbordertop = maxy*ones(1,length(xposbordertop));
xposborderbottom = xposbordertop;
yposborderbottom = miny*ones(1,length(xposborderbottom));
yposborderleft = round(linspace(1,height,10));
yposborderright = yposborderleft;
xposborderleft = minx*ones(1,length(yposborderleft));
xposborderright = maxx*ones(1,length(yposborderleft));
borderx = horzcat(xposbordertop,xposborderbottom,xposborderleft,xposborderright);
bordery = horzcat(yposbordertop,yposborderbottom,yposborderleft,yposborderright);

%Se busca para cada posicion del borde la posicion más cercana y se le
%asigna ese valor:
for i=1:length(borderx)
    for k=1:length(xpositions)
        distance(k) = sqrt((xpositions(k)-borderx(i))^2+(ypositions(k)-bordery(i))^2);
    end
    [~,nearestneighbour] = min(distance);
    borderparam(i) = parameter(nearestneighbour);
end

%Se concatenan las posiciones ingresadas del usuario con las del borde:

xpositions = vertcat(xpositions,borderx');
ypositions = vertcat(ypositions,bordery');
parameter = horzcat(parameter,borderparam);

% Primera interpolación chota para equiespaciar datos:
Surf1 = scatteredInterpolant(xpositions,ypositions,parameter');
[xq,yq] = ndgrid(0:width/10:width,0:height/10:height);
vq = Surf1(xq,yq);

%Interpolacion fina para chetizar mapa:
[xf,yf] = ndgrid(1:1:width,1:1:height);
G = griddedInterpolant(xq,yq,vq,'spline');

%Agregado de zonas blancas para rellenar tamaño de plano
mapzone = G(xf,yf);
leftspace = zeros(minx,height);
rightspace = zeros(imagesize(1)-maxx,height);
topspace = zeros(imagesize(1),imagesize(2)-maxy);
bottomspace = zeros(imagesize(1),miny);
centerzone = vertcat(leftspace,mapzone,rightspace);
map = horzcat(bottomspace,centerzone,topspace);

end