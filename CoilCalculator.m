function [B,Bx,By,Bz,H,Hx,Hy,Hz,TotalWireLength,MagneticFlux,Inductance,Resistance]=CoilCalculator(CoilLength,CoilNumber,CoilLayers,CoilDiameter,WireDiameter,CoilPrecision,CalculationZone,CalculationZoneSTEP,RelativePermeability,Current,MaterialResistivity)
tic
mu0=4*pi*10e-7;
mu=mu0*RelativePermeability;
alpha = 2*pi/CoilPrecision;
x(1:CoilPrecision*CoilNumber*CoilLayers)=NaN;
y(1:CoilPrecision*CoilNumber*CoilLayers)=NaN;
z(1:CoilPrecision*CoilNumber*CoilLayers)=NaN;
n=0;
    if WireDiameter < CoilLength/CoilNumber   
        for CurrentLayer=1:CoilLayers
            for j=1:CoilPrecision*CoilNumber
                  n=n+1;
                  x(n)=((CoilDiameter/2)+(WireDiameter*CurrentLayer))*cos(alpha*j);
                  y(n)=((CoilDiameter/2)+(WireDiameter*CurrentLayer))*sin(alpha*j);
                  z(n)=(-CoilLength/2)+(CoilLength*j)/(CoilPrecision*CoilNumber);
            end
        end
        figure
        plot3(x,y,z)
        title(['Solenoid generated, total points: ', num2str(n)])
        grid
    else
        return
    end
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
if CalculationZoneSTEP~='-' && CalculationZone(4)<CalculationZone(3)/2 && CalculationZone(4)<CalculationZone(2)/2 && CalculationZone(4)<CalculationZone(1)/2
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
lx=size(GridX,1);
ly=size(GridX,2);
lz=size(GridX,3);
H(lx,ly,lz)=zeros;
Hx(lx,ly,lz)=zeros;
Hy(lx,ly,lz)=zeros;
Hz(lx,ly,lz)=zeros;
B(lx,ly,lz)=zeros;
    for i=1:lx
        for j=1:ly
            for k=1:lz
                for m=1:length(x)-1
                    if mod(m,CoilNumber*CoilPrecision)~=0
                        Lx=x(m+1)-x(m);
                        Ly=y(m+1)-y(m);
                        Lz=z(m+1)-z(m);
                        L=[Lx Ly Lz];
                        Rx=GridX(i,j,k)-x(m);
                        Ry=GridY(i,j,k)-y(m);
                        Rz=GridZ(i,j,k)-z(m);
                        R=[Rx Ry Rz];
                        LcrossR=cross(L,R);
                        LLength=norm(L);
                        Lfi=LcrossR./LLength;
                        RLength=norm(R);
                        if RLength>(WireDiameter*2)
                            HVector=(Current/(4*pi)).*(LcrossR./RLength^3);
                        else
                            HVector=(Current*RLength)/(2*pi*(WireDiameter/2)).*Lfi;
                        end
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
Bx=Hx.*mu;
By=Hy.*mu;
Bz=Hz.*mu;
    for i=1:lx
        for j=1:ly
            for k=1:lz
                dBV=[Bx(i,j,k),By(i,j,k),Bz(i,j,k)];
                dB=norm(dBV);
                B(i,j,k)=dB;
                dHV=[Hx(i,j,k),Hy(i,j,k),Hz(i,j,k)];
                dH=norm(dHV);
                H(i,j,k)=dH;
            end
        end
    end
 MagneticFlux=0;
 CrosssectionB=0;
 Inductance=0;
 if CalculationZoneSTEP~='-'
    for i=1:lx
         for j=1:ly
             if abs(GridX(i,j,round(lz/2)))<CoilDiameter/2 && abs(GridY(i,j,round(lz/2)))<CoilDiameter/2
                    CrosssectionB=CrosssectionB+B(i,j,round(lz/2));
             end
         end
    end
    Area=(WireDiameter/2*CoilLayers+CoilDiameter^2)*pi/2;
    MagneticFlux=Area*CrosssectionB*CalculationZone(4);
    Inductance=(CoilNumber*CoilLayers)*MagneticFlux/Current;
 end
if CalculationZoneSTEP~='-'
hold on
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
