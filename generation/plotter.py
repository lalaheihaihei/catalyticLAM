from matplotlib import pyplot as plt
from ase import Atoms
from ase.visualize.plot import plot_atoms
from pymatgen.io.ase import AseAtomsAdaptor

class SlabPlotter:
    def __init__(self, rows=8, cols=4):
        """
        Initialize the plotter with a specific grid size.
        """
        self.rows = rows
        self.cols = cols
        self.fig = None
        self.axs = None
        self.current_ax_index = 0

    def get_ase_atoms(self, pmg_structure):
        """Convert a pymatgen Structure or Slab to an ASE Atoms object."""
        adaptor = AseAtomsAdaptor()
        return adaptor.get_atoms(pmg_structure)

    def plot_slabs_with_side_view_ase(self, ads_structs, material_id, miller_index):
        """
        Plot the structures on the current figure, creating a new one if necessary,
        and add Miller index above each subplot.

        :param ads_structs: A list of pymatgen Structure objects for plotting.
        :param material_id: The material ID, used for saving the figure.
        :param miller_index: The Miller index of the slab to be displayed on the plot.
        """
        if self.fig is None:
            self.fig, self.axs = plt.subplots(self.rows, self.cols, figsize=(self.cols * 3, self.rows * 3))
            self.axs = self.axs.flatten()  # Flatten the axes array for easier access
        
        for ads_struct in ads_structs:
            ase_atoms = self.get_ase_atoms(ads_struct)
            if self.current_ax_index < len(self.axs):
                ax = self.axs[self.current_ax_index]
                plot_atoms(ase_atoms, ax, rotation='90x')
                ax.title.set_text(f"Miller index: {miller_index}")  # Add Miller index
                self.current_ax_index += 1
            if self.current_ax_index < len(self.axs):
                ax = self.axs[self.current_ax_index]
                plot_atoms(ase_atoms, ax, rotation='0x,0y,0z')
                ax.title.set_text(f"Miller index: {miller_index}")  # Add Miller index
                self.current_ax_index += 1
    def plot_slabs_with_side_view_ase_uudd(self, ads_structs, material_id, miller_index):
        """
        Plot the structures on the current figure, creating a new one if necessary,
        and add Miller index above each subplot.

        :param ads_structs: A list of pymatgen Structure objects for plotting.
        :param material_id: The material ID, used for saving the figure.
        :param miller_index: The Miller index of the slab to be displayed on the plot.
        """
        if self.fig is None:
            self.fig, self.axs = plt.subplots(self.rows, self.cols, figsize=(self.cols * 3, self.rows * 3))
            self.axs = self.axs.flatten()  # Flatten the axes array for easier access

        for ads_struct in ads_structs:
            ase_atoms = self.get_ase_atoms(ads_struct)
            if self.current_ax_index < len(self.axs):
                ax = self.axs[self.current_ax_index]
                plot_atoms(ase_atoms, ax, rotation='90x')
                ax.title.set_text(f"Miller index: {miller_index}")  # Add Miller index
                self.current_ax_index += 1
            if self.current_ax_index < len(self.axs):
                ax = self.axs[self.current_ax_index]
                plot_atoms(ase_atoms, ax, rotation='0x,0y,0z')
                ax.title.set_text(f"Miller index: {miller_index}")  # Add Miller index
                self.current_ax_index += 1

    def save_figure(self, material_id,molecule_type):
        """Save the current figure to a file and reset the plotter for the next use."""
        if self.fig:
            plt.tight_layout()
            filename = f"slabs-{material_id['name']}-{material_id['mp_id']}-{molecule_type}.png"
            self.fig.savefig(filename, dpi=300)
            plt.close(self.fig)
            self.fig = None
            self.axs = None
            self.current_ax_index = 0

class BulkPlotter:
    def __init__(self, rows=7, cols=4):
        """
        Initialize the BulkPlotter with a grid of subplots.

        :param rows: Number of rows in the subplot grid.
        :param cols: Number of columns in the subplot grid.
        :param figsize: Tuple indicating figure size.
        """
        self.fig, self.axs = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3))
        self.current_ax_index = 0  # Track the current axis for plotting
        if rows * cols == 1:
            self.axs = [self.axs]  # Make axs iterable for a single subplot
        else:
            self.axs = self.axs.flatten()  # Flatten in case of a grid

    def plot_bulk(self, bulk_struct, title="Bulk"):
        """
        Plot a bulk structure 

        :param bulk_struct: A pymatgen Structure object representing the bulk material.
        :param title: Title of the plot.
        """
        if self.current_ax_index < len(self.axs):
            ase_atoms = self.get_ase_atoms(bulk_struct)
            ax = self.axs[self.current_ax_index]
            plot_atoms(ase_atoms, ax, rotation='30x,-30y,30z')  # Adjust the rotation as needed
            # Include the chemical formula in the title
            formula = bulk_struct.composition.reduced_formula
            ax.set_title(f"{title}: {formula}") 
            self.current_ax_index += 1
        else:
            print("Reached the maximum number of plots. Increase rows/cols or reset the index.")
    
    def get_ase_atoms(self, pmg_structure):
        """
        Convert a pymatgen Structure to an ASE Atoms object.
        """
        adaptor = AseAtomsAdaptor()
        return adaptor.get_atoms(pmg_structure)

    def save_figure(self, filename="bulk_structures.png"):
        """
        Save the current figure containing all plotted bulk structures to a file.
        
        :param filename: Name of the file to save the figure.
        """
        self.fig.savefig(filename)
        plt.close(self.fig)  # Close the figure to free up memory
