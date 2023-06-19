Visualisation.ipynb: 
Jupyter Notebook to visualise different phonon DOS using plotly.go

getdos.conf, mode_follow.conf:
conf-files for phonopy each for different purposes. 
Example: "phonopy -p getdos.conf" to get the DOS

opt1.traj:
is the trajectory file from geometry optimisation BFGS. It is on the 1400 structures STDEV1E-2 with a fmax of 0.1

mode_follow.bash:
Bash-script to create a desired amount of modulated structures.


#### Making MPOSCAR-files into a db with named keys
```
from ase.io.vasp import read_vasp
from ase.db import connect

db = connect('nymposcar.db')
path = '/home/armin/Downloads/Exjobb/Basic_Tutorial/MPOSCAR/'
for i in range(58, 73):
    for j in np.arange(-1.5, 1.6, 0.5):
        try:
            atoms = read_vasp(path + 'MPOSCAR' + str(i) + str(round(j, 1)))
            db.write(atoms, key=f'Mod {i}, Amp {round(j, 1)}')
        except:
            pass
```

#### Checking largest distance between two atoms
###### Taking account for periodicity 
```
path = '/home/armin/Downloads/Exjobb/Basic_Tutorial/MPOSCAR/'
distances = []
for i in range(58, 73):
    for j in np.arange(-1.5, 1.6, 0.5):
        try:
            atoms1 = read(path + 'MPOSCARorig' + str(i) + str(round(j, 1)))
            atoms2 = read(path + 'MPOSCAR' + str(i) + str(round(j, 1)))

            if len(atoms1) == len(atoms2):
                cell_lengths = atoms1.cell.cellpar()[:3]
                max_distance = min(cell_lengths) / 2.0

                for k in range(len(atoms1)):
                    atom1 = atoms1[k]
                    atom2 = atoms2[k]
                    atoms = Atoms([atom1, atom2])
                    distance = atoms.get_distance(0, 1)
                    if distance > max_distance:
                        atoms.wrap(pbc=True)
                        atoms.minimum_image()
                        distance = atoms.get_distance(0, 1)
                    distances.append(distance)
            else:
                print(':(')
        except:
            pass
print(max(distances))
```

### Choosing x structures from database.
```
from ase.db import connect
import random 

db = connect('tillPeter/1-FP_PBE_U4_stdev_1e-2.db')
struct_id = [row.id for row in db.select()]
selected = random.sample(struct_id, 20)
new_db = connect('1-TwentyFP2.db')

for struct_id in selected:
    structures = db.get(id=struct_id)
    new_db.write(structures)
```
