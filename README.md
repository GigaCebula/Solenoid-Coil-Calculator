# Solenoid-Coil-Calculator
A short MATLAB function calculating magnetic field parameters of an air-core solenoid coil using Biot-Savart law as a basis

![fig1](https://github.com/user-attachments/assets/6af4477c-427c-4134-b00d-e4dc92d716ca)


```
[B,Bx,By,Bz,H,Hx,Hy,Hz,TotalWireLength,MagneticFlux,Inductance,Resistance]=CoilCalculator(CoilLength,CoilNumber,CoilLayers,CoilDiameter,WireDiameter,CoilPrecision,CalculationZone,RelativePermeability,Current,MaterialResistivity)
```
- CoilLenght - total length of the solenoid
- CoilNumber - amount of individual coils in a single layer of the solenoid
- CoilLayers - number of wire coil layers
- CoilDiameter - diameter of first layer
- WireDiameter - diameter of the coil wire
- CoilPrecision - number of corners 1 coil has, higher number gives higher calculation accuracy (I recommend to set it at least 25)
- CalculationZone - 4 dimentional vector that sets up calculation space: \
    `[X, Y, Z STEP]`
- CalculationZoneSTEP - if set as `-` calculation space is set to calculate magnetic field properties only in the point specified by coordinates X, Y, Z. If STEP is set to any number smaller than the smallest of other 3 values, program will create calculation mesh where X, Y and Z are dimensions of a cuboid with center in 0,0,0 and STEP as the distance between calculation points.
- RelativePermeability - relative magnetic permeability of the space around the solenoid
- Current - current flowing through the wire
- MaterialResistivity - resistivity of the wire

- B - scalar magnetic flux density
- Bx, By, Bz - magnetic flux density vector components
- H - scalar magnetic field intensity
- Hx, Hy, Hz - magnetic field intensity vector components
- TotalWireLength - total length of solenoid wire
- MagneticFlux - total magnetic flux going through the coil
- Inductance, Resistance - self explanatory

NOTE: MagneticFlux and Inductance are calculated only if user sets up the function to calculate magnetic field around the coil and won't be initiated otherwise

# Usage
Download the .m file and run it as any other function in either MATLAB or GNU Octave
