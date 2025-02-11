function [B,Bx,By,Bz,H,Hx,Hy,Hz,TotalWireLength,MagneticFlux,Inductance,Resistance]=CoilCalculator(CoilLength,CoilNumber,CoilLayers,CoilDiameter,WireDiameter,CoilPrecision,CalculationZone,RelativePermeability,Current,MaterialResistivity)
tic
mu0=4*pi*1e-7;
mu=mu0*RelativePermeability;
alpha = 2*pi/CoilPrecision;

% Deklaracja macierzy współrzędnych

x(1:CoilPrecision*CoilNumber*CoilLayers)=NaN;
y(1:CoilPrecision*CoilNumber*CoilLayers)=NaN;
z(1:CoilPrecision*CoilNumber*CoilLayers)=NaN;
n=0;

% Tworzenie geometrii cewki 

    if WireDiameter < CoilLength/CoilNumber   
        for CurrentLayer=1:CoilLayers
            for j=1:CoilPrecision*CoilNumber
                  n=n+1;
                  x(n)=((CoilDiameter/2)+(WireDiameter*CurrentLayer))*cos(alpha*j);
                  y(n)=((CoilDiameter/2)+(WireDiameter*CurrentLayer))*sin(alpha*j);
                  z(n)=(-CoilLength/2)+(CoilLength*j)/(CoilPrecision*CoilNumber);
            end
        end

        % Wykres wizualizacji cewki

        figure
        plot3(x,y,z)
        title(['Solenoid generated, total points: ', num2str(n)])
        grid
    else
        return
    end

% Obliczenie długości przewodu i rezystancji cewki

PointDistance(1:length(x)-1)=zeros;
        for i=1:length(x)-1
            if mod(i,CoilNumber*CoilPrecision)~=0
            A=[x(i),y(i),z(i)];
            B=[x(i+1),y(i+1),z(i+1)];
            PointDistance(i)=norm(B-A);
            end
        end
TotalWireLength=sum(PointDistance);
Resistance=MaterialResistivity*(TotalWireLength/(WireDiameter^2*pi));

% Tworzenie macierzy strefy obliczeniowej

if CalculationZone(4)~=0 && CalculationZone(4)<CalculationZone(3)/2 ...
    && CalculationZone(4)<CalculationZone(2)/2 && CalculationZone(4)<CalculationZone(1)/2
    StepsX=abs(CalculationZone(1)/CalculationZone(4));
    StepsY=abs(CalculationZone(2)/CalculationZone(4));
    StepsZ=abs(CalculationZone(3)/CalculationZone(4));
    spacex=linspace(-CalculationZone(1)/2,CalculationZone(1)/2,StepsX);
    spacey=linspace(-CalculationZone(2)/2,CalculationZone(2)/2,StepsY);
    spacez=linspace(-CalculationZone(3)/2,CalculationZone(3)/2,StepsZ);
    [GridX,GridY,GridZ]=meshgrid(spacex,spacey,spacez);
else
    GridX=CalculationZone(1);
    GridY=CalculationZone(2);
    GridZ=CalculationZone(3);
end

% Deklaracja składowych obliczeniowych natężenia pola magnetycznego

lx=size(GridX,1);
ly=size(GridX,2);
lz=size(GridX,3);
H(lx,ly,lz)=zeros;
Hx(lx,ly,lz)=zeros;
Hy(lx,ly,lz)=zeros;
Hz(lx,ly,lz)=zeros;
B(lx,ly,lz)=zeros;

% Pętla obliczania pola magnetycznego

    for i=1:lx
        for j=1:ly
            for k=1:lz
                for m=1:length(x)-1
                    if mod(m,CoilNumber*CoilPrecision)~=0

                        % Przygotowywanie wymaganych danych

                        Lx=x(m+1)-x(m);
                        Ly=y(m+1)-y(m);
                        Lz=z(m+1)-z(m);
                        L=[Lx Ly Lz];
                        Rx=GridX(i,j,k)-x(m);
                        Ry=GridY(i,j,k)-y(m);
                        Rz=GridZ(i,j,k)-z(m);
                        R=[Rx Ry Rz];
                        LcrossR=cross(L,R);
                        RLength=norm(R);
                        Lfi=LcrossR./RLength;

                        % Sprawdzenie czy punkt obliczeniowy znajduje się
                        % wewnątrz przewodu

                        Theta = acos(dot(R,L))/(norm(R)*norm(L));
                        RCheck = RLength*abs(sin(Theta));

                        if RCheck>(WireDiameter/2) || RLength>sqrt(norm(L)^2+(WireDiameter/2)^2)

                            % Jeżeli poza przewodem - prawo Biota-Savarda

                            HVector=(Current/(4*pi)).*(LcrossR./RLength^3);
                        else

                           HVector=((Current/norm(L))*RLength)/(2*pi*(WireDiameter/2)).*Lfi;

                        end

                        % Zsumowanie obliczonych wartości do natężenia pola
                        % w danym punkcie
                            dHx=HVector(1);
                            dHy=HVector(2);
                            dHz=HVector(3);

                    Hx(i,j,k)=Hx(i,j,k)+dHx;
                    Hy(i,j,k)=Hy(i,j,k)+dHy;
                    Hz(i,j,k)=Hz(i,j,k)+dHz;
                    end
                end
            end
        end
    end

% Obliczenie indukcji pola magnetycznego

Bx=Hx.*mu;
By=Hy.*mu;
Bz=Hz.*mu;
    for i=1:lx
        for j=1:ly
            for k=1:lz

                % Połączenie składowych x, y, z w jedną macierz

                dBV=[Bx(i,j,k),By(i,j,k),Bz(i,j,k)];
                dB=norm(dBV);
                B(i,j,k)=dB;
                dHV=[Hx(i,j,k),Hy(i,j,k),Hz(i,j,k)];
                dH=norm(dHV);
                H(i,j,k)=dH;
            end
        end
    end

% Obliczenie indukcyjności

 MagneticFlux=0;
 CrosssectionB=0;
 Inductance=0;
 if CalculationZone(4)~=0
    for i=1:lx
         for j=1:ly
             if abs(GridX(i,j,round(lz/2)))<CoilDiameter/2 && abs(GridY(i,j,round(lz/2)))<CoilDiameter/2
                    CrosssectionB=CrosssectionB+Bz(i,j,round(lz/2));
             end
         end
    end
    Area=(WireDiameter/2*CoilLayers+CoilDiameter^2)*pi/2;
    MagneticFlux=Area*CrosssectionB*CalculationZone(4);
    Inductance=(CoilNumber*CoilLayers)*MagneticFlux/Current;
 end

% Utworzenie wykresów

if CalculationZone(4)~=0
    hold on
    
    % Wykres strzałkowy pola magnetycznego
    
    quiver3(GridX,GridY,GridZ,Hx,Hy,Hz)
    legend("Solenoid coil","Magnetic field intensity vectors")
    hold off
    figure
    tiledlayout(1,2)
    nexttile
    MapX=permute(GridX(1,:,1),[2 1 3]);
    MapY=permute(GridZ(1,1,:),[3 2 1]);
    MapData=permute(H(:,round(ly/2),:),[3 1 2]);
    MapData2=permute(H(:,:,round(lz/2)),[1 2 3]);

    % Mapy przekrojów cewki

    imagesc(-MapX,-MapY,MapData);
    colorbar;
    title("Magnetic field norm [A/m] in XZ plane")
    colormap turbo
    nexttile
    imagesc(-MapX,-MapX,MapData2);
    colorbar;
    title("Magnetic field norm [A/m] in XY plane")
    colormap turbo
    figure
    tiledlayout(1,2)
    nexttile
    MapData3=permute(B(:,round(ly/2),:),[3 1 2]);
    MapData4=permute(B(:,:,round(lz/2)),[1 2 3]);
    imagesc(-MapX,-MapY,MapData3);
    colorbar;
    title("Magnetic flux density norm [T] in XZ plane")
    colormap turbo
    nexttile
    imagesc(-MapX,-MapX,MapData4);
    colorbar;
    title("Magnetic flux density norm [T] in XY plane")
    colormap turbo
    toc
end
end
