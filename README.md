# Solenoid-Coil-Calculator
A short MATLAB function calculating magnetic field parameters of an air-core solenoid coil using Biot-Savart law as a basis
```
[B,Bx,By,Bz,H,Hx,Hy,Hz,TotalWireLength,Inductance,Resistance]=CoilCalculator(CoilLength,CoilNumber,CoilLayers,CoilDiameter,CoilPrecision,CalculationZone,RelativePermeability,Current,MaterialResistivity)
```
- CoilLenght - total length of the solenoid
- CoilNumber - amount of individual coils in a single layer of the solenoid
- CoilLayers - number of wire coil layers
- CoilDiameter - diameter of first layer
- CoilPrecision - number of corners 1 coil has, higher number gives higher calculation accuracy (I recommend to set it at least 25)
- CalculationZone - 4 dimentional vector that sets up calculation space: \
    `[X, Y, Z, STEP]`\
STEP - if set as `-` calculation space is set to calculate magnetic field properties only in the point specified by coordinates X, Y, Z. If STEP is set to any number smaller than the smallest of other 3 values, program will create calculation mesh where X, Y and Z are sizes of a cuboid with center in 0,0,0 and STEP as the step between calculation points.
- RelativePermeability - relative magnetic permeability of the space around the solenoid
- Current - current flowing through the wire
- MaterialResistivity - resistivity of the wire

- B - scalar magnetic flux density
- Bx, By, Bz - magnetic flux density vector components
- H - scalar magnetic field intensity
- Hx, Hy, Hz - magnetic field intensity vector components
- TotalWireLength - total length of solenoid wire
- Inductance, Resistance - self explanatory


# Usage
Download the .m or .p file and run it as a function in either MATLAB or GNU Octave
