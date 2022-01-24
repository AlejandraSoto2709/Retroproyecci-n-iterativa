%DESCRIPCION: Genera una imagen reconstruida mediante el metodo iterativo
%de retroproyeccion. Se presenta en una figura la imagen original y la
%reconstruida. 
%INPUTS:
% NombreImagen -> Nombre de la imagen (string)
% NoIteracion  -> Cantidad de iteraciones
%EJEMPLO:
% Retroproyeccion('Picture0.jpg',10)
%ALEJANDRA SOTO 

function []=Retroproyeccion(NombreImagen,NoIteracion)

    %% LEER IMAGENES
    imOrig = imread(NombreImagen);
    imInfo = imfinfo(NombreImagen);
    %[image, map] = rgb2ind(imOrig,imInfo.NumberOfSamples);
    imRows = imInfo.Height;
    imCols = imInfo.Width;

    %% CONVERSION DE IMAGEN A ESCALA DE GRISES
    imGray = rgb2gray(imOrig);

    %% VALORES DE PROYECCION IMAGEN ORGINAL 
    pDeg0 = zeros(1,imRows);
    pDeg90 = zeros(1,imCols);
    PDeg45 = zeros(1,imRows + imCols - 1);
    pDeg135 = zeros(1,imRows + imCols -1);

    %Proyeccion a 0°
    for i = 1:imRows
        pDeg0(i) = sum(imGray(i,:));
    end    

    %Proyeccion a 90°
    for i = 1:imCols 
        pDeg90(i) = sum(imGray(:,i));
    end

    %Proyeccion a 45°
    pDeg45 =fliplr(sum(rot90(spdiags(imGray))'));

    %Proyeccion a 135°
    pDeg135 = sum(rot90(spdiags(fliplr(imGray)))'); 


    %% METODO DE RECONSTRUCCION
    %Inicializacion para metodo de reconstruccion
    iteracion = ones(imRows,imCols);
    error = ones(imRows,imCols);
    pDeg0_it = zeros(1,imRows);
    pDeg90_it = zeros(1,imCols);
    pDeg45_it = zeros(1,imRows + imCols - 1);
    pDeg135_it = zeros(1,imRows + imCols -1);

    %Figura para animacion
    % figure(1)
    % title('Reconstruccion')
    % set(gcf,'Color','w')

    %Reconstruccion
    
    for k = 1:NoIteracion % K = Numero de iteraciones
        %Estimacion de proyeccion a 0°
        for i = 1:imRows
            pDeg0_it(i) = sum(iteracion(i,:));
        end    

        %Estimacion de proyeccion a 90°
        for i = 1:imCols 
            pDeg90_it(i) = sum(iteracion(:,i));
        end

        %Estimacion de proyeccion a 45°
        pDeg45_it =fliplr(sum(rot90(spdiags(iteracion))'));

        %Proyeccion a 135°
        pDeg135_it = sum(rot90(spdiags(fliplr(iteracion)))'); 

        %Calculo de matriz error - suma de 0° y 90°
        for i = 1:imRows
            for j = 1:imCols
               error(i,j) = pDeg90(j)/pDeg90_it(j) + pDeg0(i)/pDeg0_it(i); 
            end    
        end

        errorDiags = spdiags(error);
        [r,c] = size(errorDiags);
        ind =  fliplr(-c+1:1:c-1);

        %Suma de las proyecciones en 45°
        for i= 1:c %Recorre columnas
            for j = 1:r %Recorre Renglones
                if errorDiags(j,i) ~= 0
                    errorDiags(j,i) =  errorDiags(j,i) + pDeg45(i)/pDeg45_it(i);
                end    
            end    
        end

        %Reacomodo de matriz cuadrada
        if imRows == imCols
            for i = 1:imRows
                error(:,i) = nonzeros(errorDiags(i,:))'; 
            end    
            error = flipud(error);
            errorDiags = spdiags(rot90((error)));
        end

        %Reacomodo matriz con mayor numero de columnas
        if imRows < imCols
            for i = 1:imRows
                error(i,:) = nonzeros(errorDiags(i,:)); 
            end    
            errorDiags = spdiags(rot90((error)));
        end    

        %Reacomodo matriz con mayor numero de renglones
        if imRows > imCols
            for i = 1:imCols
                error(:,i) = nonzeros(errorDiags(i,:))'; 
            end
            error = flipud(error);
            errorDiags = spdiags(rot90((error)));
        end    

        %Suma de las proyecciones en 135°
        for i= 1:c %Recorre columnas
            for j = 1:r %Recorre Renglones
                if errorDiags(j,i) ~= 0
                    errorDiags(j,i) =  errorDiags(j,i) + pDeg135(i)/pDeg135_it(i);
                end    
            end    
        end

        %Reacomodo de matriz cuadrada
        if imRows == imCols
            for i = 1:imRows
                error(:,i) = nonzeros(errorDiags(i,:))'; 
            end    
            error = error'./4;
        end

        %Reacomodo matriz con mayor numero de columnas
        if imRows < imCols
            for i = 1:imRows
                error(i,:) = nonzeros(errorDiags(i,:)); 
            end    
            error = error./4 ; 
        end    

        %Reacomodo matriz con mayor numero de renglones
        if imRows > imCols
            for i = 1:imCols
                error(:,i) = nonzeros(errorDiags(i,:))'; 
            end
            error = fliplr(error)./4; 
        end 

        %Imagen reconstruida iteracion k 
        iteracion = iteracion .* error;

        %Animacion de la reconstruccion
    %     imshow(mat2gray(iteracion))
    %     pause(0.01)

    end

    figure(2)
    subplot(1,2,1)
    imshow(imOrig)
    title('Imagen Original')

    subplot(1,2,2)
    imshow(mat2gray(iteracion))
    sub = "Iteracion " + NoIteracion + " ";
    title(sprintf('Imagen Reconstruida - ' + sub)) 

    set(gcf,'Color','w')

end