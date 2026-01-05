from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib import rcParams
import os
from svgutils import transform,compose
import numbers
import math

def is_valid_smiles(smiles_string):
    """
    Check if a string is a valid SMILES string by trying to create a molecule from it.
    
    Parameters
    ----------
    smiles_string : str
        String to check if it's a valid SMILES
        
    Returns
    -------
    bool
        True if the string is a valid SMILES, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        return mol is not None
    except:
        return False

class strSVG(compose.Element):
    """SVG from string.

    Parameters
    ----------
    fname : str
       full path to the file
    """

    def __init__(self, svg):
        obj = transform.fromstring(svg)
        self.root = obj.getroot().root
        
from IPython.display import SVG # /!\ note the 'SVG' function also in svgutils.compose
import numpy as np
from rdkit import Geometry


def draw_with_indeces(settings):
    """
    Drawing function that displays the input smiles string with all atom indeces
    """
    m = Chem.MolFromSmiles(settings['SMILESSTRING'])
    dm = Draw.PrepareMolForDrawing(m)
    d2d = Draw.MolDraw2DSVG(350,350)
    opts = d2d.drawOptions()
    for i in range(m.GetNumAtoms()):
        opts.atomLabels[i] = m.GetAtomWithIdx(i).GetSymbol()+str(i)
    d2d.DrawMolecule(dm)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def compute_label_pos(center, radius, theta_start, theta_end, vertical_offset=0.3):
    # Compute label position above substituent wedge
    # theta is in degrees, convert to radians
    theta = math.radians((theta_start + theta_end) / 2)
    x = center[0] + (radius + vertical_offset) * math.cos(theta)
    y = center[1] + (radius + vertical_offset) * math.sin(theta)
    return x, y


class RadialScope(object):
    """
        Radial Scope Plot using Matplotlib and RDKIT

        Written by Simon Duerr with code from Greg Landrum for the automatic positioning (@rdkit)

        License: MIT 

        This class handles all the heavy work and constructs a colorbar svg and a pie chart svg for each passed radial scope dictionary and then assembles the svg. 

        The size of the matplotlib plot is always the same, which is why the hard coded center used for positioning on the atoms is likely fine.
        
        Enhanced with multiple colormap support and reversed color mapping (lower values = darker colors).
    """
    
    # Available colormaps for different color schemes
    AVAILABLE_COLORMAPS = {
        'blues': ['Blues', 'Blues_r', 'GnBu', 'PuBu'],
        'greens': ['Greens', 'Greens_r', 'BuGn', 'GnBu'],
        'reds': ['Reds', 'Reds_r', 'OrRd', 'PuRd'],
        'purples': ['Purples', 'Purples_r', 'PuBu', 'BuPu'],
        'oranges': ['Oranges', 'Oranges_r', 'OrRd', 'RdOr'],
        'grays': ['Greys', 'Greys_r', 'gray', 'gist_gray'],
        'rainbow': ['viridis', 'plasma', 'inferno', 'magma'],
        'diverging': ['RdBu', 'RdBu_r', 'RdYlBu', 'RdYlBu_r', 'Spectral', 'coolwarm'],
        'sequential': ['Blues', 'Greens', 'Reds', 'Purples', 'Oranges', 'Greys']
    }

    @classmethod
    def get_available_colormaps(cls, category=None):
        """
        Get available colormaps by category or all colormaps.
        
        Parameters
        ----------
        category : str, optional
            Category name ('blues', 'greens', 'reds', etc.) or None for all
            
        Returns
        -------
        list or dict
            List of colormap names if category specified, dict of all categories if None
        """
        if category is None:
            return cls.AVAILABLE_COLORMAPS
        return cls.AVAILABLE_COLORMAPS.get(category, [])
    
    @classmethod
    def print_colormap_options(cls):
        """Print all available colormap options organized by category."""
        print("Available Colormap Categories:")
        print("=" * 40)
        for category, colormaps in cls.AVAILABLE_COLORMAPS.items():
            print(f"\n{category.upper()}:")
            for cmap in colormaps:
                print(f"  - {cmap}")
    
    def __init__(self, settings_dict,*args):
        """
        Calls the main constructur. 


        Parameters
        ----------
        settings_dicte : dict
               main settings
           args: dict
               dictionaries containing the information for the scope plots
        """
        self.settings=settings_dict
        self.plots=[]
        self.cbar_index=0 # if you want to use multiple colorbars you need to uncomment two lines below in the plot_figure_and_colorbar function at the bottom of the functuon
        for plot in args:
            self.plots.append(plot)
        if len(args)<1:
            raise Exception('minimum 1 radial scope dictionary must be given')
        self.main()


    def draw_smiles(self):
        """
            draws the smiles string for the plotting, needs to return the drawing objects also because the atom indeces are extracted.
        """
        m = Chem.MolFromSmiles(self.settings['SMILESSTRING'])
        dm = Draw.PrepareMolForDrawing(m)
        d2d = Draw.MolDraw2DSVG(350,350)
        if self.settings['use_bw_atom_theme']:
            d2d.drawOptions().useBWAtomPalette()
        d2d.DrawMolecule(dm)
        d2d.FinishDrawing()
        return d2d.GetDrawingText(), d2d, dm

    def replace_label_with_smiles(self,svg_file='', smiles='C=C', search_index='~0'):
        """
        draws a small organic rest, looks for the ~index comment in the complete svg, finds the glyphs and replaces them with the organic subsituent. 
        Note that this is hacky and the position of the organic subsituent likely needs to be fixed in a vector software such as Inkscape. 

        Parameters
        ----------
        svg_file : str
           the text of the svg file which contains the comments that are replaced
        search_index: str
            the name of the comment which is indicated by a prepended tilde by the user
        smiles: str
            the smiles string used for the replacement
        """
        m = Chem.MolFromSmiles(smiles)
        dm = Draw.PrepareMolForDrawing(m)
        d2d = Draw.MolDraw2DSVG(100,100)
        d2d.drawOptions().padding=0
        d2d.drawOptions().clearBackground=False
        d2d.drawOptions().useBWAtomPalette()
        d2d.DrawMolecule(dm)
        d2d.FinishDrawing()
        # scale smiles molecule and remove clutter
        group1 = d2d.GetDrawingText()
        replace_str=group1[group1.find('<!-- END OF HEADER -->')+len("<!-- END OF HEADER -->")+1:-8]
        replace_str='<g transform="translate(-300,300)scale(6,-6)">'+replace_str+"</g>"
        # find the index in the pie chart that needs to be replaced, we will geplace the two glyphs with the svg text from rdkit
        index_of_comment=svg_file.find(str(search_index))
        index_of_defsend=svg_file[index_of_comment:].find('</defs>')
        start=svg_file[index_of_comment+index_of_defsend:].find('<use')
        end=svg_file[index_of_comment+index_of_defsend:].find('</g>')
        item_to_replace=svg_file[index_of_comment+index_of_defsend+start:index_of_comment+index_of_defsend+end]
        return svg_file.replace(item_to_replace, replace_str)
    
        
    def plot_figure_and_colorbar(self,radial_scope_setup, vals):
        """
        plots one pie chart and the corresponding colorbar

        Parameters
        ----------
        radial_scope_setup : dict
           contains information about this pie chart like labels, colormap etc.
        vals :list
            contains 5 list with prepended input for the empty wedge in the pie chart and the processed labels.

        Returns
        -------
        fig_svg: str
            the svg for the pie plot
        colorbars
            the svg for the colorbar of the inner and outer circle

        """
        
        fig, ax = plt.subplots(1,figsize=(10,10))

        size = 0.5 # size of inner plot
        alpha = 0
        which_wedge = 0 # first wedge is always transparent
        _=ax.set(aspect="equal")
        circle1 = plt.Circle((0, 0), 0.055, color='w', ls='-', ec='k', lw=1.4, zorder=99, gid='circle_anchor')
        label = ax.annotate(radial_scope_setup['rest_label'], xy=(0, 0), fontsize=13, ha="center", va='center', zorder=100, gid='circle_content')
        
        if len(radial_scope_setup['min_max_value'])==2:
            min_cbar_value=[radial_scope_setup['min_max_value'][0][0],radial_scope_setup['min_max_value'][0][1]]
            max_cbar_value=[radial_scope_setup['min_max_value'][1][0],radial_scope_setup['min_max_value'][1][1]]
        else:
            def _numeric_bounds(values):
                numeric = [v for v in values if isinstance(v, numbers.Number)]
                if len(numeric) == 0:
                    return 0.0, 1.0
                return float(np.min(numeric)), float(np.max(numeric))
            inner_min, inner_max = _numeric_bounds(vals[1])
            outer_min, outer_max = _numeric_bounds(vals[2])
            min_cbar_value=[inner_min, outer_min]
            max_cbar_value=[inner_max, outer_max]

        cmap_inner = plt.get_cmap(radial_scope_setup['CMAPINNER'])
        cmap_outer = plt.get_cmap(radial_scope_setup['CMAPOUTER'])
        # Build numeric arrays for colormaps, mapping string entries to a light intensity
        def _to_numeric_for_colormap(values, vmin, vmax):
            light_fraction = 0.15
            numeric_values = []
            for v in values:
                if isinstance(v, numbers.Number):
                    numeric_values.append(v)
                else:
                    numeric_values.append(vmin + light_fraction * (vmax - vmin))
            return np.array(numeric_values, dtype=float)

        inner_numeric = _to_numeric_for_colormap(vals[1], min_cbar_value[0], max_cbar_value[0])
        outer_numeric = _to_numeric_for_colormap(vals[2], min_cbar_value[1], max_cbar_value[1])

        # Reverse the colormap so lower values get darker colors and higher values get lighter colors
        # We do this by reversing the normalization: higher values get lower normalized values
        norm_outer = plt.Normalize(min_cbar_value[1], max_cbar_value[1])
        # Reverse the mapping: max_value -> 0, min_value -> 1
        outer_norm_reversed = norm_outer(outer_numeric)
        outer_colors = cmap_outer(outer_norm_reversed)
        
        norm_inner = plt.Normalize(min_cbar_value[0], max_cbar_value[0])
        # Reverse the mapping: max_value -> 0, min_value -> 1
        inner_norm_reversed = norm_inner(inner_numeric)
        inner_colors = cmap_inner(inner_norm_reversed)

        labels_circle=ax.pie(vals[0], startangle=radial_scope_setup['startangle'], radius=0.5, colors=['w']*len(vals[0]), labels=vals[5], labeldistance=1.85,
               wedgeprops=dict(width=size, edgecolor='k',linewidth= 1.4), textprops=dict(fontsize='large', weight="semibold",va='center', ha='center') )


        outer_circle=ax.pie(vals[0], radius=0.7,startangle=radial_scope_setup['startangle'], colors=outer_colors, labels=vals[3], labeldistance=0.55,
               wedgeprops=dict(width=size, edgecolor='k',linewidth= 1.4),textprops=dict(fontsize='large',weight="semibold",va='center', ha='center'))

        inner_circle=ax.pie(vals[0], startangle=radial_scope_setup['startangle'], radius=1-size, colors=inner_colors, labels=vals[4], labeldistance=1.2,
               wedgeprops=dict(width=size, edgecolor='k',linewidth= 1.4,), textprops=dict(fontsize='large', weight="semibold",va='center', ha='center') )

        # # Rotate labels to align with wedges (parallel to wedge edges)
        # num_wedges = len(vals[0]) - 1
        # if num_wedges > 0:
        #     label_startangle = radial_scope_setup.get('startangle', 0)
        #     label_coverangle = radial_scope_setup['coverangle_wedges']
        #     angle_step = label_coverangle / num_wedges
            
        #     # Rotate outer circle labels
        #     for i, text in enumerate(outer_circle[1]):
        #         if i > 0 and i <= num_wedges:
        #             angle = label_startangle + (i - 0.5) * angle_step + 90 + 150
        #             text.set_rotation(angle)
            
        #     # Rotate inner circle labels
        #     for i, text in enumerate(inner_circle[1]):
        #         if i > 0 and i <= num_wedges:
        #             angle = label_startangle + (i - 0.5) * angle_step + 90 + 150
        #             text.set_rotation(angle)
            
        #     # Rotate labels circle labels
        #     for i, text in enumerate(labels_circle[1]):
        #         if i > 0 and i <= num_wedges:
        #             angle = label_startangle + (i - 0.5) * angle_step + 90 + 150
        #             text.set_rotation(angle)

        # Annotate substituent labels above each wedge if provided
        if 'substituent_labels' in radial_scope_setup and radial_scope_setup['substituent_labels']:
            num_wedges = len(vals[0]) - 1  # Exclude the first transparent wedge
            # Use label-specific angles if provided, otherwise use plot angles
            label_startangle = radial_scope_setup.get('label_startangle', radial_scope_setup.get('startangle', 0))
            label_coverangle = radial_scope_setup.get('label_coverangle', radial_scope_setup['coverangle_wedges'])
            angle_step = label_coverangle / num_wedges
            center = (0, 0)
            radius = 0.75  # Position labels at the outer edge where structures are placed
            for i, label in enumerate(radial_scope_setup['substituent_labels']):
                if i < num_wedges:  # Only label actual wedges, not the transparent one
                    theta_start = label_startangle + i * angle_step
                    theta_end = theta_start + angle_step
                    lx, ly = compute_label_pos(center, radius, theta_start, theta_end, vertical_offset=0.42)
                    ax.text(lx, ly, label, ha='center', va='bottom', fontsize=13, weight='bold', zorder=50)

        inner_circle[0][which_wedge].set_alpha(alpha)
        outer_circle[0][which_wedge].set_alpha(alpha)
        labels_circle[0][which_wedge].set_alpha(alpha)

        for i in range(len(inner_circle[1])):
            if outer_numeric[i]<self.settings['white_cutoff']:
                _=inner_circle[1][i].set_color('white')
            if inner_numeric[i]<self.settings['white_cutoff']:
                _=outer_circle[1][i].set_color('white')
        _=ax.add_artist(circle1)

        from io import StringIO
        
        sio = StringIO()
        _=fig.savefig(sio, transparent=True, format='SVG')
        fig_svg = sio.getvalue()
        _=plt.close(fig)
        # close first figure as we do not want to display it and start assembling the colorbar
        
        # currently only the first colorbar is used. If you need to print multiple colorbars uncomment line below
        fig, axs = plt.subplots(2, figsize=(3,2))
        
        # Create reversed colormaps for colorbar display to match the reversed mapping
        cmap_inner_reversed = plt.cm.get_cmap(radial_scope_setup['CMAPINNER'])
        cmap_outer_reversed = plt.cm.get_cmap(radial_scope_setup['CMAPOUTER'])

        sm = ScalarMappable(cmap=cmap_inner_reversed, norm=plt.Normalize(min_cbar_value[0],max_cbar_value[0]))
        _=sm.set_array([])
        cbar = plt.colorbar(sm,cax=axs[0], orientation="horizontal",ticks=[min_cbar_value[0], max_cbar_value[0]/2, max_cbar_value[0]])
        _=cbar.ax.set_xticklabels([str(min_cbar_value[0]), str(int(max_cbar_value[0]/2)), str(max_cbar_value[0])]) 
        _=cbar.set_label(radial_scope_setup['INNERLABEL'],weight='bold', fontsize=12)
        
        sm1 = ScalarMappable(cmap=cmap_outer_reversed, norm=plt.Normalize(min_cbar_value[1],max_cbar_value[1]))
        _=sm1.set_array([])
        cbar1 = plt.colorbar(sm1,cax=axs[1], orientation="horizontal",ticks=[min_cbar_value[1], max_cbar_value[1]/2, max_cbar_value[1]])
        _=cbar1.ax.set_xticklabels([str(min_cbar_value[1]), str(int(max_cbar_value[1]/2)), str(max_cbar_value[1])]) 
        _=cbar1.set_label(radial_scope_setup['OUTERLABEL'],weight='bold',fontsize=12)
        figure2=fig.tight_layout()
        sio2 = StringIO()
        # _=fig.savefig('colorbar'+self.cbar_index+'.svg', transparent=True, format='SVG')
        # self.cbar_index+=1
        _=fig.savefig(sio2, transparent=True, format='SVG')
        colorbars = sio2.getvalue()
        _=plt.close(fig)
        return fig_svg, colorbars

    def main(self):
        settings=self.settings
        SMILESSTRING=settings['SMILESSTRING']
        resulting_plots=[]
        pRList=[]
        mol_svg, d2d, dm=self.draw_smiles()
        replace_index=[]
        for scope_plot in self.plots:
            # for each scope plot, make a vals list containing empty first items for the wedge with alpha=0
            if type(scope_plot)!=dict:
                continue
            
            sizes=[360-scope_plot['coverangle_wedges']]+[scope_plot['coverangle_wedges']/scope_plot['no_wedges']]*scope_plot['no_wedges']
            
            label_inner_circle, label_outer_circle=['']+['']*scope_plot['no_wedges'],['']+['']*scope_plot['no_wedges']
            if (len(scope_plot['value_inner_circle'])!=scope_plot['no_wedges'] or len(scope_plot['value_outer_circle'])!=scope_plot['no_wedges']):
                print('not equal')
            value_inner_circle, value_outer_circle=scope_plot['value_inner_circle'], scope_plot['value_outer_circle']
            # Handle rounding boundaries - support both single value and separate inner/outer boundaries
            if isinstance(scope_plot['rounding_boundary'], (list, tuple)) and len(scope_plot['rounding_boundary']) == 2:
                # Separate boundaries for inner and outer circles
                rounding_boundary_inner = scope_plot['rounding_boundary'][0]
                rounding_boundary_outer = scope_plot['rounding_boundary'][1]
            else:
                # Single boundary for both circles (backward compatibility)
                rounding_boundary_inner = scope_plot['rounding_boundary']
                rounding_boundary_outer = scope_plot['rounding_boundary']
            
            value_groups=scope_plot['value_groups']
            
            for i in range(scope_plot['no_wedges']):
                inner_val = value_inner_circle[i]
                outer_val = value_outer_circle[i]
                # If value is numeric, support rounding; if string, use as-is
                if isinstance(inner_val, numbers.Number):
                    # if scope_plot['rounding'] and inner_val >= rounding_boundary_inner:
                    #     label_inner_circle[i+1] = ">"+str(inner_val)
                    if scope_plot['rounding'] and inner_val == rounding_boundary_inner:
                        label_inner_circle[i+1] = "<"+str(inner_val)
                    else:
                        label_inner_circle[i+1] = str(inner_val)
                else:
                    label_inner_circle[i+1] = str(inner_val)

                if isinstance(outer_val, numbers.Number):
                    # if scope_plot['rounding'] and outer_val >= rounding_boundary_outer:
                    #     label_outer_circle[i+1] = ">"+str(outer_val)
                    if scope_plot['rounding'] and outer_val == rounding_boundary_outer:
                        label_outer_circle[i+1] = "<"+str(outer_val)
                    else:
                        label_outer_circle[i+1] = str(outer_val)
                else:
                    label_outer_circle[i+1] = str(outer_val)
                    
            # Add an offset for unique indices based on plot order
            offset = len(replace_index)  # Start new indices after existing ones
            for i, item in enumerate(value_groups):
                if item[0] == '~':
                    # Handle existing ~ prefixed SMILES
                    replace_index.append(('~' + str(offset), item[1:]))
                    value_groups[i] = '~' + str(offset)
                    offset += 1
                elif is_valid_smiles(item):
                    # Handle SMILES strings without ~ prefix
                    replace_index.append(('~' + str(offset), item))
                    value_groups[i] = '~' + str(offset)
                    offset += 1

            vals = [sizes, # size of the wedges, the first wedge is transparent and will not be shown 
                    [0]+value_inner_circle, # colormap values for the inner circle, maximum value determines intensity, first is for the transparent wedge and should stay 0
                    [0]+value_outer_circle, # colormap values for the outer circle, maximum value determines intensity, first is for the transparent wedge and should stay 0
                    label_inner_circle, #labels for the inner circle
                    label_outer_circle, #labels for the outer circle    
                    [""]+value_groups, #groups  
                   ]
            resulting_plots.append(self.plot_figure_and_colorbar(scope_plot,vals))

            # get the position from the settings - support both atom_id and bond_id
            if 'attach_bond_id' in scope_plot and scope_plot['attach_bond_id'] is not None:
                # Position on bond - get midpoint between two atoms
                bond_id = scope_plot['attach_bond_id']
                bond = dm.GetBondWithIdx(bond_id)
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                pos1 = dm.GetConformer().GetAtomPosition(atom1_idx)
                pos2 = dm.GetConformer().GetAtomPosition(atom2_idx)
                # Calculate midpoint
                mid_pos = Geometry.Point2D((pos1.x + pos2.x) / 2, (pos1.y + pos2.y) / 2)
                pRList.append(d2d.GetDrawCoords(mid_pos))
            else:
                # Position on atom (original behavior)
                rIdx = scope_plot['attach_atom_id']
                pRList.append(d2d.GetDrawCoords(Geometry.Point2D(dm.GetConformer().GetAtomPosition(rIdx))))
        
        # take colorbar from first plot  #ToDo extension to multiple colorbars
        colorbar=compose.Panel(strSVG(resulting_plots[0][1]).scale(0.8).move(-350,400))
        panels=[compose.Panel(strSVG('<svg></svg>'))]*len(resulting_plots)
        for i,plot in enumerate(resulting_plots):
            panels[i]=strSVG(resulting_plots[i][0]).move(-369,-358).scale(1).move(pRList[i].x,pRList[i].y)
        compose.Figure("720", "720", 
            compose.Panel(strSVG(mol_svg).scale(1).move(0,0)),
            colorbar,
            *panels
            ).move(350,350).scale(self.settings['scalefactor']).save("substrate_scope.svg")
        new_svg=SVG('substrate_scope.svg')._data
        for item in replace_index:
            new_svg=self.replace_label_with_smiles(svg_file=new_svg, smiles=item[1],search_index= item[0] )

        if settings['use_bold_font']:
            new_svg.replace('font-weight:normal', 'font-weight:bold')
        f = open("substrate_scope_replaced.svg","w") 
        f.write(new_svg)
        f.close()
        print('File written to:',os.getcwd()+'/substrate_scope_replaced.svg')
        
        
        

        
        
