from ase.io import read, write
from ase.constraints import FixAtoms
import numpy as np

# ================== STEP 1: Read the structure ==================
atoms = read("mapbi3_slab.cif")

# ================== STEP 2: Analyze layers properly ==================
print("="*50)
print("STRUCTURE ANALYSIS (MAPbI3 slab - Proper layer detection):")
print("="*50)

z_positions = atoms.positions[:, 2]
# Sort and group atoms into distinct layers with better tolerance
z_sorted = np.sort(z_positions)
z_diff = np.diff(z_sorted)
# Find big gaps between layers (gap > 1 Å usually between PbI planes)
layer_boundaries = [0]
for i in range(len(z_diff)):
    if z_diff[i] > 1.0:  # adjust this if your layer spacing is different
        layer_boundaries.append(i+1)

print(f"Detected {len(layer_boundaries)} distinct layers along Z")

for i in range(len(layer_boundaries)):
    start = layer_boundaries[i]
    end = layer_boundaries[i+1] if i+1 < len(layer_boundaries) else len(z_sorted)
    z_mean = np.mean(z_sorted[start:end])
    n_atoms = end - start
    symbols = atoms[[j for j in range(len(atoms)) if z_positions[j] in z_sorted[start:end]]].get_chemical_symbols()
    print(f"Layer {i+1:2d}: Z ≈ {z_mean:6.2f} Å  |  {n_atoms:3d} atoms  |  {set(symbols)}")

print("="*50)

# ================== STEP 3: Ask how many BOTTOM layers to fix ==================
while True:
    try:
        n_layers_to_fix = int(input("Evlo bottom layers fix pannanum? (1 or 2): "))
        if n_layers_to_fix in [1, 2]:
            break
        print("1 or 2 மட்டுமே enter pannu bro!")
    except ValueError:
        print("Number dhaan venum (1 or 2)!")

# ================== STEP 4: Fix the bottom N layers properly ==================
# Take the first N layers from bottom
fixed_indices = []
for i in range(n_layers_to_fix):
    start = layer_boundaries[i]
    end = layer_boundaries[i+1] if i+1 < len(layer_boundaries) else len(atoms)
    layer_indices = [j for j in range(len(atoms)) if z_positions[j] in z_sorted[start:end]]
    fixed_indices.extend(layer_indices)

constraint = FixAtoms(indices=fixed_indices)
atoms.set_constraint(constraint)

print("\n" + "="*50)
print(f"✓ Fixed {len(fixed_indices)} atoms in bottom {n_layers_to_fix} layer(s)")
print(f"✓ Total atoms: {len(atoms)} | Fixed: {len(fixed_indices)} | Free: {len(atoms)-len(fixed_indices)}")
print("="*50)

# ================== STEP 5: Save files ==================
write("mapbi3_slab_fixed.xyz", atoms)
write("POSCAR", atoms, format="vasp", vasp5=True, direct=True, sort=True)

# Basic QE input with correct constraints
qe_input = f"""&control
    calculation = 'relax'
    prefix = 'mapbi3_slab'
    pseudo_dir = './pseudo/'
    outdir = './out/'
/
&system
    ibrav = 0
    nat = {len(atoms)}
    ntyp = 4
    ecutwfc = 60.0
    ecutrho = 480.0
/
&electrons
/
&ions
/
CELL_PARAMETERS angstrom
{atoms.get_cell().tolist()}
ATOMIC_POSITIONS angstrom
"""

for i, atom in enumerate(atoms):
    x, y, z = atom.position
    if i in fixed_indices:
        qe_input += f"   {atom.symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}  0 0 0\n"
    else:
        qe_input += f"   {atom.symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}  1 1 1\n"

with open("mapbi3_slab_fixed.in", "w") as f:
    f.write(qe_input)

print("Saved:")
print("  - mapbi3_slab_fixed.xyz")
print("  - POSCAR (VASP)")
print("  - mapbi3_slab_fixed.in (QE) ✅")
print("\nIppo bottom layers சரியா fix ஆகிருச்சு! Relax பண்ணு bro! 🎉")
