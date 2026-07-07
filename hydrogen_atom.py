"""
Hydrogen Atom Wavefunction Visualizer
Author: Anubhav
Physics: Quantum Mechanics — Atomic Orbitals
Method: Analytical wavefunctions using SciPy special functions

Plots:
  1. Probability density cross-sections for 6 orbitals
  2. Radial probability distribution R(r) for multiple states
  3. 3D-style orbital shapes (1s, 2p, 3d)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from scipy.special import sph_harm_y, genlaguerre, factorial
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────
# PHYSICAL CONSTANTS (atomic units: a₀ = 1)
# ─────────────────────────────────────────────

a0 = 1.0   # Bohr radius in atomic units

# ─────────────────────────────────────────────
# WAVEFUNCTION FUNCTIONS
# ─────────────────────────────────────────────

def radial_wavefunction(n, l, r):
    """
    Radial wavefunction R_nl(r) for hydrogen atom.
    Uses associated Laguerre polynomials.
    r in units of Bohr radius a0.
    """
    # Normalization constant
    norm = np.sqrt(
        (2.0 / (n * a0))**3 *
        factorial(n - l - 1) /
        (2 * n * factorial(n + l)**3)
    )

    rho = 2.0 * r / (n * a0)

    # Associated Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
    L = genlaguerre(n - l - 1, 2 * l + 1)

    R = norm * np.exp(-rho / 2) * rho**l * L(rho)
    return R


def wavefunction_3d(n, l, m, x, y, z):
    """
    Full 3D wavefunction psi_nlm(x,y,z).
    Returns probability density |psi|^2.
    """
    r     = np.sqrt(x**2 + y**2 + z**2)
    r     = np.where(r == 0, 1e-10, r)  # avoid division by zero
    theta = np.arccos(np.clip(z / r, -1, 1))
    phi   = np.arctan2(y, x)

    R   = radial_wavefunction(n, l, r)
    Y   = sph_harm_y(l, m, theta, phi)

    psi = R * Y
    return np.abs(psi)**2


def wavefunction_2d_slice(n, l, m, grid_size=400, box=20):
    """
    2D cross-section of |psi|^2 in the xz-plane (y=0).
    Returns x, z grids and probability density.
    """
    x = np.linspace(-box, box, grid_size)
    z = np.linspace(-box, box, grid_size)
    X, Z = np.meshgrid(x, z)
    Y    = np.zeros_like(X)

    prob = wavefunction_3d(n, l, m, X, Y, Z)
    return X, Z, prob


# ─────────────────────────────────────────────
# PLOT 1: ORBITAL CROSS-SECTIONS (6 orbitals)
# ─────────────────────────────────────────────

def plot_orbital_crosssections():
    orbitals = [
        (1, 0, 0,  '1s',   15),
        (2, 0, 0,  '2s',   30),
        (2, 1, 0,  '2p₀',  30),
        (3, 0, 0,  '3s',   50),
        (3, 1, 0,  '3p₀',  50),
        (3, 2, 0,  '3d₀',  50),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(
        "Hydrogen Atom — Orbital Probability Densities |ψ|²\n"
        "Cross-section in xz-plane (y = 0)",
        fontsize=14, fontweight='bold', y=0.98
    )

    # Custom colormap: black → purple → white (like electron cloud)
    colors_list = ['#000000', '#1a0050', '#4a0080', '#8800cc',
                   '#cc44ff', '#ffaaff', '#ffffff']
    cmap = LinearSegmentedColormap.from_list('orbital', colors_list)

    for ax, (n, l, m, label, box) in zip(axes.flat, orbitals):
        X, Z, prob = wavefunction_2d_slice(n, l, m, grid_size=300, box=box)

        # Normalize for better visualization
        prob_normalized = prob / prob.max()

        im = ax.contourf(X, Z, prob_normalized,
                         levels=100, cmap=cmap)
        ax.set_aspect('equal')
        ax.set_title(f'n={n}, l={l}, m={m}  |  {label}',
                     fontsize=11, fontweight='bold', color='white')
        ax.set_xlabel('x (a₀)', fontsize=9, color='white')
        ax.set_ylabel('z (a₀)', fontsize=9, color='white')
        ax.tick_params(colors='white')
        ax.set_facecolor('black')
        for spine in ax.spines.values():
            spine.set_edgecolor('white')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('|ψ|² (normalized)', fontsize=8, color='white')
        cbar.ax.yaxis.set_tick_params(color='white')
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

    fig.patch.set_facecolor('#0a0a0a')
    plt.tight_layout()
    plt.savefig('/home/claude/hydrogen_orbitals.png',
                dpi=150, bbox_inches='tight', facecolor='#0a0a0a')
    print("Saved: hydrogen_orbitals.png")


# ─────────────────────────────────────────────
# PLOT 2: RADIAL PROBABILITY DISTRIBUTION
# ─────────────────────────────────────────────

def plot_radial_distributions():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        "Hydrogen Atom — Radial Probability Distribution P(r) = r²|R_nl(r)|²",
        fontsize=13, fontweight='bold'
    )

    r = np.linspace(0.01, 60, 2000)

    # Panel 1: n=1,2,3 s-orbitals (l=0)
    ax = axes[0]
    configs = [(1,0,'#E53935'),(2,0,'#1E88E5'),(3,0,'#43A047')]
    for n, l, color in configs:
        R   = radial_wavefunction(n, l, r)
        P   = r**2 * R**2
        ax.plot(r, P, color=color, linewidth=2, label=f'n={n}, l=0 (s)')
        # Mark most probable radius
        r_max = r[np.argmax(P)]
        ax.axvline(x=r_max, color=color, linestyle='--', alpha=0.4)
    ax.set_title('s orbitals (l=0)', fontsize=11, fontweight='bold')
    ax.set_xlabel('r (a₀)', fontsize=10)
    ax.set_ylabel('P(r) = r²|R|²', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 60)

    # Panel 2: p-orbitals (l=1)
    ax = axes[1]
    configs = [(2,1,'#E53935'),(3,1,'#1E88E5'),(4,1,'#43A047')]
    for n, l, color in configs:
        R   = radial_wavefunction(n, l, r)
        P   = r**2 * R**2
        ax.plot(r, P, color=color, linewidth=2, label=f'n={n}, l=1 (p)')
    ax.set_title('p orbitals (l=1)', fontsize=11, fontweight='bold')
    ax.set_xlabel('r (a₀)', fontsize=10)
    ax.set_ylabel('P(r) = r²|R|²', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 60)

    # Panel 3: d-orbitals (l=2)
    ax = axes[2]
    configs = [(3,2,'#E53935'),(4,2,'#1E88E5'),(5,2,'#43A047')]
    for n, l, color in configs:
        R   = radial_wavefunction(n, l, r)
        P   = r**2 * R**2
        ax.plot(r, P, color=color, linewidth=2, label=f'n={n}, l=2 (d)')
    ax.set_title('d orbitals (l=2)', fontsize=11, fontweight='bold')
    ax.set_xlabel('r (a₀)', fontsize=10)
    ax.set_ylabel('P(r) = r²|R|²', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 60)

    plt.tight_layout()
    plt.savefig('/home/claude/hydrogen_radial.png',
                dpi=150, bbox_inches='tight')
    print("Saved: hydrogen_radial.png")


# ─────────────────────────────────────────────
# PLOT 3: ENERGY LEVEL DIAGRAM
# ─────────────────────────────────────────────

def plot_energy_levels():
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor('#0d1117')
    ax.set_facecolor('#0d1117')

    ax.set_title(
        "Hydrogen Atom — Energy Level Diagram\nEₙ = -13.6 eV / n²",
        fontsize=13, fontweight='bold', color='white', pad=15
    )

    colors = {0:'#64B5F6', 1:'#81C784', 2:'#FFB74D', 3:'#F06292'}
    labels = {0:'s', 1:'p', 2:'d', 3:'f'}

    for n in range(1, 6):
        E_n = -13.6 / n**2
        for l in range(min(n, 4)):
            x_pos = l * 1.5 + 1
            color = colors.get(l, '#FFFFFF')

            # Energy level line
            ax.hlines(E_n, x_pos - 0.4, x_pos + 0.4,
                      color=color, linewidth=2.5)

            # Label
            ax.text(x_pos, E_n + 0.15,
                    f'n={n},{labels[l]}',
                    ha='center', va='bottom', fontsize=8,
                    color=color, fontweight='bold')

        # Energy label on right
        ax.text(8, E_n, f'E_{n} = {E_n:.2f} eV',
                va='center', fontsize=9, color='white')

    ax.set_xlim(0, 9)
    ax.set_ylim(-14.5, 1)
    ax.set_ylabel('Energy (eV)', fontsize=11, color='white')
    ax.set_xticks([1, 2.5, 4, 5.5])
    ax.set_xticklabels(['s (l=0)', 'p (l=1)', 'd (l=2)', 'f (l=3)'],
                       fontsize=10, color='white')
    ax.tick_params(colors='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(y=0, color='white', linestyle='--',
               alpha=0.4, label='Ionization limit (0 eV)')
    ax.grid(True, alpha=0.15, axis='y')

    # Legend for colors
    for l, label in labels.items():
        if l < 4:
            ax.plot([], [], color=colors[l], linewidth=3,
                    label=f'l={l} ({label})')
    ax.legend(fontsize=9, loc='upper right',
              facecolor='#1a1a2e', labelcolor='white')

    plt.tight_layout()
    plt.savefig('/home/claude/hydrogen_energy_levels.png',
                dpi=150, bbox_inches='tight', facecolor='#0d1117')
    print("Saved: hydrogen_energy_levels.png")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 55)
    print("  Hydrogen Atom Wavefunction Visualizer")
    print("  Anubhav | BSc Physics | Utkal University")
    print("=" * 55)

    print("\n[1/3] Generating orbital cross-sections...")
    plot_orbital_crosssections()

    print("[2/3] Generating radial distributions...")
    plot_radial_distributions()

    print("[3/3] Generating energy level diagram...")
    plot_energy_levels()

    print("\n✓ All plots generated successfully!")
    print("\nFiles created:")
    print("  → hydrogen_orbitals.png      (6 orbital shapes)")
    print("  → hydrogen_radial.png        (radial probability)")
    print("  → hydrogen_energy_levels.png (energy diagram)")
