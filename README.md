## Radial Scope plot
### The matplotlib/rdkit version

## How does it work

1. **Prepare a SMILES string**  
   Get a SMILES string from ChemDraw or any other software able to output SMILES strings (such as [PubChem](https://pubchem.ncbi.nlm.nih.gov/edit3/index.html)).  
   For every \(R_x\) group put a methyl group.

2. **Set up the global settings**  
   In the notebook `radial_scope_plot.ipynb`, define the `settings` dictionary (e.g. `SMILESSTRING`, `white_cutoff`, `scalefactor`, etc.).

3. **Set up per‑substituent scope dictionaries**  
   For each \(R_x\) you want to display (e.g. `radial_scope_setup_R1`, `radial_scope_setup_R2`, …), create a dictionary as shown in the notebook with:
   - **`rest_label`**: the label that appears in the central circle (e.g. `R$_1$`),
   - **`no_wedges` / `coverangle_wedges` / `startangle`**: geometry of the wedges,
   - **`value_inner_circle` / `value_outer_circle`**: numeric (or string) values to plot,
   - **`value_groups`**: SMILES strings or `~index` placeholders for each substituent,
   - **`attach_atom_id`** or **`attach_bond_id`**: where on the core the scope plot is attached,
   - optional **`substituent_labels`**, **`label_startangle`**, **`label_coverangle`** to control text labels.

4. **Run the notebook**  
   Execute all cells in `radial_scope_plot.ipynb` (Run All, or stepwise with Ctrl+Enter).  
   The main call
   ```python
   scope_plot = rs(settings, radial_scope_setup_R1, radial_scope_setup_R2, ...)
   ```
   will generate `substrate_scope.svg` and `substrate_scope_replaced.svg` in the project directory.

## Installation

- **Use a Conda-based Python**  
  Use this within a Python Anaconda distribution, e.g. [Miniconda](https://docs.conda.io/en/latest/miniconda.html).  
  Anaconda is available on Windows, macOS, and Linux, and this project should work on all platforms.

- **Clone this repository**

  ```bash
  git clone https://github.com/Aman-Aryan11/RadialScopePlot.git
  cd RadialScopePlot
  ```

- **Create an environment containing RDKit, Matplotlib, and svgutils**

  ```bash
  conda env create -f environment.yml
  conda activate substrate-scope-plot
  ```

- **Start Jupyter Notebook**

  ```bash
  jupyter-notebook
  ```

  Then open `radial_scope_plot.ipynb` and follow the instructions in the cells.

## Customising colour maps and labels

- **Reversing / un‑reversing the colormap in `radialscope.py`**

  The `RadialScope` class applies a *reversed* colour mapping so that **lower values appear darker and higher values appear lighter**.  
  If you want to *restore the normal (non‑reversed) behaviour* (i.e. higher values → darker colours), remove the explicit reversal logic in `radialscope.py` and use the direct normalised values instead.

  Concretely, in `radialscope.py` inside `plot_figure_and_colorbar`, locate the block where `inner_numeric` / `outer_numeric` are normalised (around the `norm_inner` / `norm_outer` lines). Remove any call that reverses the normalisation (for example a `.reverse()` or equivalent custom reversal) and pass the normalised arrays directly into the colormaps:

  ```python
  norm_outer = plt.Normalize(min_cbar_value[1], max_cbar_value[1])
  outer_colors = cmap_outer(norm_outer(outer_numeric))

  norm_inner = plt.Normalize(min_cbar_value[0], max_cbar_value[0])
  inner_colors = cmap_inner(norm_inner(inner_numeric))
  ```

  If your local version still uses a `.reverse()` operation on the colormap or on the normalised values, **remove that `.reverse()` call** so that the colormap is no longer flipped.

- **Fixing label alignment with wedges**

  If the substituent labels drawn around the circle do not visually align with the wedges, you can disable the automatic label‑placement block in `radialscope.py`.  
  In `plot_figure_and_colorbar`, comment out the label‑rotation section between lines 281 and 304:

  ```python
  # for i, text in enumerate(outer_circle[1]):
  #     ...
  # for i, text in enumerate(inner_circle[1]):
  #     ...
  # for i, text in enumerate(labels_circle[1]):
  #     ...
  ```

  After commenting this block, labels will remain in their default Matplotlib orientation, which is often easier to control manually via `label_startangle`, `label_coverangle`, and the positions of `substituent_labels`.

## Additional tips

- **Passing SMILES directly in `value_groups`**  
  You can now pass SMILES strings directly in `value_groups`; they will automatically be detected and rendered as molecular structures in the final SVG. Existing `~index` placeholders are still supported and are filled via `replace_label_with_smiles`.

- **Using atom vs bond attachment**  
  For radial scopes attached to bonds instead of atoms, use `attach_bond_id` (preferred when substituents conceptually sit on a bond) instead of `attach_atom_id`. The code will compute the bond midpoint and place the scope more symmetrically.

- **Colormap options**  
  Call:
  ```python
  from radialscope import RadialScope as rs
  rs.print_colormap_options()
  ```
  to list all supported colormap families and pick schemes that best match your data (e.g. blues, greens, diverging maps, etc.).
