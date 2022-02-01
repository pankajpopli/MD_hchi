# LD non-affinity.

## Checking traingle <-> rectangle transformation using hchi.
For this only the files mentioned in the makefile are needed.
This code can be heavily optimised, for instance,
- Deformation matrix can be calculated once and then used in the calculation of X and h_chi_forces.
- Y matrix can be calculated at the beigining of code as it is a geometric quantity and remains constant.
- The particle structure can be shorten by removing unwanted class, eg. ri0 which was originally used in MC.
- Can avoid repeated structure calls. Call it once and store its value in local variable to be used inside the function.
- various other small things.