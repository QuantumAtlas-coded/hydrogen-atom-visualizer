# Hydrogen Atom Wavefunction Visualizer

**Author:** Anubhav  
**Physics:** Quantum Mechanics — Atomic Orbitals  
**Language:** Python (NumPy, SciPy, Matplotlib)

---

## What This Does
Visualizes hydrogen atom electron probability clouds using 
the exact analytical solution to the Schrodinger equation.

## The Wavefunction
psi(n,l,m) = R_nl(r) x Y_lm(theta, phi)

| Symbol | Name | Meaning |
|--------|------|---------|
| n | Principal | Energy level (1,2,3...) |
| l | Azimuthal | Shape: s=sphere, p=dumbbell, d=clover |
| m | Magnetic | Orientation in space |

## Plots Generated
- Orbital probability densities (1s, 2s, 2p, 3s, 3p, 3d)
- Radial probability distributions P(r) = r2 x R(r)2
- Energy level diagram: En = -13.6 eV / n2
- ## Plots Generated
- Orbital probability densities (1s, 2s, 2p, 3s, 3p, 3d)
- Radial probability distributions P(r)
- Energy level diagram

### Orbital Shapes
![Orbitals](hydrogen_orbitals.png)

### Radial Probability
![Radial](hydrogen_radial.png)

### Energy Levels
![Energy](hydrogen_energy_levels.png)

## How to Run
```bash
pip install numpy scipy matplotlib
python hydrogen_atom.py
```

## References
- Griffiths - Introduction to Quantum Mechanics
- Sakurai - Modern Quantum Mechanics
- NIST Atomic Spectra Database
