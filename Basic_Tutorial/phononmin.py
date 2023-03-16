import os, json
import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import read_crystal_structure
from ase import Atoms
from ase.io import write, read
from ase.calculators.mixing import LinearCombinationCalculator
from ccs_fit.ase_calculator.ccs_ase_calculator import CCS

def setup_model(CCS_params, **kwargs):
    CCS_calc = CCS(CCS_params=CCS_params)
    calcs = [CCS_calc]
    weights = [1]
    calc = LinearCombinationCalculator(calcs, weights)
    return calc

def run_gap(calc, phonon):
    supercells = phonon.get_supercells_with_displacements()
    print('Number of Cells', len(supercells))
    set_of_forces = []
    for scell in supercells:
        cell = Atoms(symbols=scell.get_chemical_symbols(),
                     scaled_positions=scell.get_scaled_positions(),
                     cell=scell.get_cell(),
                     pbc=True)
        cell.set_calculator(calc)
        forces = cell.get_forces()
        drift_force = forces.sum(axis=0)
        for force in forces:
            force -= drift_force / forces.shape[0]
        set_of_forces.append(forces)
    return set_of_forces

def obtain_phonopy_mesh(phonon, mesh):
    phonon.run_mesh(mesh)
    mesh_dict = phonon.get_mesh_dict()
    return mesh_dict['qpoints'], mesh_dict['weights'], mesh_dict['frequencies']

def print_freq_at_G(phonon):
    print('')
    print("[Phonopy] Phonon frequencies at Gamma:")
    for i, freq in enumerate(phonon.get_frequencies((0, 0, 0))):
        print("[Phonopy] =: .5f THz" %  (i + 1, freq)) # THz

def obtain_phonon_dispersion_bands(phonon,bands_ranges,band_resolution=50,band_connection=False):
    bands = []
    for q_start, q_end in bands_ranges:
        band = []
        for i in range(band_resolution + 1):
            band.append(np.array(q_start) + (np.array(q_end) - np.array(q_start)) / band_resolution * i)
        bands.append(band)
    phonon.run_band_structure(bands, is_band_connection=band_connection, with_eigenvectors=True)

    bands_dict = phonon.get_band_structure_dict()
    return (bands_dict['qpoints'],
            bands_dict['distances'],
            bands_dict['frequencies'],
            bands_dict['eigenvectors'])

def run_dos(phonon, mesh):
    phonon.run_mesh(mesh)
    phonon.run_total_dos(use_tetrahedron_method=True)
    dos_dict = phonon.get_total_dos_dict()
    dos = np.array([dos_dict['frequency_points'], dos_dict['total_dos']])
    np.savetxt("TDOS", dos)
    print("[Phonopy] Phonon DOS:")
    phonon.plot_total_dos().savefig("dos.png")

def run_pdos(phonon, mesh):
    phonon.run_mesh(mesh, with_eigenvectors=True, is_mesh_symmetry=False)
    phonon.run_projected_dos()
    pdos_dict = phonon.get_projected_dos_dict()
    for i in range(len(pdos_dict['projected_dos'])):
        pdos = np.array([pdos_dict['frequency_points'], pdos_dict['projected_dos'][i]])
        np.savetxt(str(i) + ".pdos", pdos)
    phonon.plot_projected_dos().show()

def run_thermal(phonon,mesh):
    phonon.run_mesh(mesh)
    phonon.run_thermal_properties(t_step=10, t_max=1700, t_min=0)
    phonon.write_yaml_thermal_properties()
    tp_dict = phonon.get_thermal_properties_dict()
    temperatures = tp_dict['temperatures']
    free_energy = tp_dict['free_energy']
    entropy = tp_dict['entropy']
    heat_capacity = tp_dict['heat_capacity']
    #for t, F, S, cv in zip(temperatures, free_energy, entropy, heat_capacity):
    #    print(("%12.3f " + "%15.7f" * 3) % ( t, F, S, cv ))
    phonon.plot_thermal_properties().show()  

def main():
    with open("CCS_params.json", "r") as f:
        CCS_params = json.load(f)
    calc = setup_model(CCS_params)
    unitcell = read('LJ.db', do_not_split_by_at_sign=True) #vasp-format can only store 1 Atoms object.
    #write('POSCAR', unitcell)
    cell, _ = read_crystal_structure('POSCAR', interface_mode='vasp')
    N = 3
    smat = [[N, 0, 0], [0, N, 0], [0, 0, N]]
    phonon = Phonopy(cell, smat)
    phonon.generate_displacements(distance=0.02)
    set_of_forces = run_gap(calc, phonon)
    phonon.produce_force_constants(forces=set_of_forces)
    mesh = (20, 20, 20)
    #run_dos(phonon, mesh)
    #run_pdos(phonon, mesh)
    #run_thermal(phonon, mesh)
    #os.chdir(base_dir+"/RUN_SIM/")obtain_phonon_dispersion_bands(phonon, )

if __name__ == "__main__":
    main()