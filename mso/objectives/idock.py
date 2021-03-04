import subprocess as sp
import numpy as np

#from rdkit import Chem
#from region_mutate import coordinatesPDB
#mport builder
#from builder import MGLTools
import argparse
import os, sys
import time
import uuid
import shutil
import http.client
import mimetypes
import uuid


"""
Docking Routine

1. perform docking with vina
2. output parse
3. docking energy recording

"""

# MGL tools in DOCKER
MGL_RECEPTOR_PATH = "/root/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
MGL_LIGAND_PATH = "/root/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"


def tfold_cli_submit(cmd, task_type, task_tag, cpu, memory, user, 
                     image="docker.oa.com/seven-developer/tencent_threading:latest", 
                     log_path="/dev/null", ceph="realbig", fast_mode=True, 
                     tfold_ms_host="100.109.129.154"):
    fast_mode_str=""
    if fast_mode:
        fast_mode_str="&mode=fast"
    url = "/enqueue"
    conn = http.client.HTTPConnection("100.109.129.154", 8080)
    payload = 'kind=tfoldtask&event={}&task-type={}&key={}&cpu={}&memory={}&user={}&log={}&image={}&ceph={}{}'.format(cmd, task_type, task_tag, cpu, memory, user, log_path, image, ceph, fast_mode_str)
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    while True:
        conn.request("PUT", url, payload, headers)
        response = conn.getresponse()
        if response.status == 200:
            print(response.read().decode('utf-8'))
            break
        else:
            print(response.read().decode('utf-8'))
       	    time.sleep(2)


class MGLTools(object):
    """
    Prepare the ligand and receptor structures for AutoDock Vina with official MGL-tools.

    References
    ----------
    MGLTools http://mgltools.scripps.edu/

    Methods
    -------
    prepare_ligand4vina
    prepare_receptor4vina

    Parameters
    ----------

    Examples
    --------
    >>> # using MGL tools for ligand format converting
    >>> import deepunion
    >>> mgl = deepunion.builder.MGLTools()
    >>> mgl.prepare_ligand4vina("mol1.pdb", "mol1.pdbqt")

    """

    def __init__(self):
        pass

    def prepare_ligand4vina(self, ligand_fn, out_fn,
        script="/Library/MGLTools/1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"):

        cmd = "%s -l %s -o %s -A bonds_hydrogens -g -s -U nphs_lps" %(script, ligand_fn, out_fn)
        try:
            job = sp.Popen(cmd, shell=True)
            job.communicate()
        except:
            print("Prepare molecule error: ", ligand_fn)

        return self

    def prepare_receptor4vina(self, receptor_fn, out_fn,
        script="/Library/MGLTools/1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"):

        cmd = "%s -r %s -o %s -A bonds_hydrogens -e False" %(script, receptor_fn, out_fn)
        try:
            job = sp.Popen(cmd, shell=True)
            job.communicate()
        except:
            print("Prepare molecule error: ", receptor_fn)

        return self


def babel_converter(input, output, babelexe="obabel", mode="general"):
    """
    Convert molecules with openbabel tool.

    Parameters
    ----------
    input : str,
        The input molecule file name.
    output : str,
        The output molecule file name.
    babelexe : str,
        The openbabel exe file.
    mode : str, default is general.
        Method to convert molecules.
        Options:
        general -- generate a molecule based on the previous provided input
        AddPolarH -- generate a new molecule by adding polar hydrogen atoms.
            Polar hydrogens are important for AutoDock Vina docking.

    Returns
    -------

    """

    cmd = ""
    if mode == "general":
        cmd = "%s %s -O %s" % (babelexe, input, output)
        job = sp.Popen(cmd, shell=True)
        job.communicate()
    elif mode == "AddPolarH":
        cmd = "%s %s -O %s -d" % (babelexe, input, output+"_temp_noH.pdbqt")
        job = sp.Popen(cmd, shell=True)
        job.communicate()
        cmd = "%s %s -O %s --AddPolarH" % (babelexe, output+"_temp_noH.pdbqt", output)
        job = sp.Popen(cmd, shell=True)
        job.communicate()
        os.remove(output+"_temp_noH.pdbqt")
    else:
        pass

    return None



def _get_xyz_mol2(filename):
    coords = []
    with open(filename) as lines:
        condition = False
        for s in lines:
            if "@<TRIPOS>ATOM" in s:
                condition = True
            elif "@<TRIPOS>BOND" in s:
                condition = False
            elif condition:
                _xyz = [float(x) for x in s.split()[2:5]]
                coords.append(_xyz)

    return coords


def _get_xyz_pdb(filename):
    coords = []
    with open(filename) as lines:
        lines = [x for x in lines if len(x) > 5 and x.split()[0] in ['ATOM', 'HETATM']]

        for s in lines:
            try:
                _x = float(s[30:38].strip())
                _y = float(s[38:46].strip())
                _z = float(s[46:54].strip())

                coords.append([_x, _y, _z])
            except:
                print("INFO: Parse pdb file error: ", s)

    return coords


class VinaDocking(object):

    def __init__(self, vina_exe="vina"):
        self.vina_exe = vina_exe

        self.config = None
        self.output = None

    def vina_config(self, receptor, ligand, outname,
                    n_cpus, exhaustiveness, center,
                    boxsize=[30, 30, 30], logfile="log.log",
                    n_modes=1, config="vina.config"):
        self.config = config

        config_dir = os.path.dirname(self.config)
        os.makedirs(config_dir, exist_ok=True)

        with open(self.config, "w") as tofile:
            tofile.write("receptor = %s \n" % os.path.abspath(receptor))
            tofile.write("ligand = %s \n" % os.path.abspath(ligand))
            tofile.write("out = %s \n" % os.path.abspath(outname))

            # center of x y z
            tofile.write("center_x = %.3f \n" % center[0])
            tofile.write("center_y = %.3f \n" % center[1])
            tofile.write("center_z = %.3f \n" % center[2])
            # box size of x y z
            tofile.write("size_x = %.2f \n" % boxsize[0])
            tofile.write("size_y = %.2f \n" % boxsize[1])
            tofile.write("size_z = %.2f \n" % boxsize[2])

            if "idock" not in self.vina_exe:
                tofile.write("cpu = %d \n" % n_cpus)
                tofile.write("exhaustiveness = %d \n" % exhaustiveness)
                tofile.write("num_modes = %d \n" % n_modes)
                tofile.write("log = %s \n" % logfile)

        self.output = outname

        return self

    def run_docking(self, environment_str="export CA=10 &&", mode="general"):
        if mode == "general":

            if self.config is not None and os.path.exists(self.config):

                job = sp.Popen("%s %s --config %s " % (environment_str,
                                                       self.vina_exe,
                                                       self.config),
                            shell=True)
                job.communicate()

                job.terminate()
            else:
                print("Please setup config first")

            return ""
        elif mode == "tfoldms":
            cmd = ""
            if self.config is not None and os.path.exists(self.config):
                cmd = "%s %s --config %s" % (environment_str, self.vina_exe, os.path.abspath(self.config)),
            else:
                print("Please setup config first")

            return cmd
        else:
            return ""


class LigandPrepare(object):

    def __init__(self, ligand):
        self.ligand = ligand

    def prepare_ligand(self, out_name):
        mgl = MGLTools()
        mgl.prepare_ligand4vina(self.ligand, out_name, MGL_LIGAND_PATH)

        return self


class ReceptorPrepare(object):

    def __init__(self, receptor):

        self.receptor = receptor

    def pocket_center(self, LIG="", residues=[]):
        if len(LIG):
            print("LIGAND FILE ", LIG)
            with open(LIG) as lines:
                lig_lines = [x for x in lines if len(x.split()) and x.split()[0] in ["ATOM", "HETATM"]]

            # read coordinates
            coord = coordinatesPDB().getAtomCrdFromLines(lig_lines)
            coord = np.array(coord)
            return np.mean(coord, axis=0)

        else:
            # get center of several residues
            coord = self._get_pocket_center(residues)
            return coord

    def _get_pocket_center(self, residue_ids):

        p = mt.load_pdb(self.receptor)
        indices = []

        for resid in residue_ids:
            indices += list(p.topology.select("residue %d" % resid))

        return p.xyz[0][np.array(indices)].mean(axis=0)

    def receptor_addH(self, out_pdb="temp.pdb", method='obabel'):

        if method == "obabel":
            mol = Chem.MolFromPDBFile(self.receptor)

            Chem.AddHs(mol)

            Chem.MolToPDBFile(mol, out_pdb)
        elif method == "mgltools":

            mgl = MGLTools()
            mgl.prepare_receptor4vina(self.receptor, out_pdb, MGL_RECEPTOR_PATH)

        return self


'''def rmsd(mol1, mol2):
    # try:
    #    m1 = mt.load_pdb(mol1).xyz[0]
    #    m2 = mt.load_pdb(mol2).xyz[0]
    # except RuntimeError:
    cpdb = coordinatesPDB()
    with open(mol1) as lines:
        m1 = cpdb.getAtomCrdFromLines([x for x in lines if ("ATOM" in x or "HETATM" in x)])
    with open(mol2) as lines:
        m2 = cpdb.getAtomCrdFromLines([x for x in lines if ("ATOM" in x or "HETATM" in x)])

    rmsd = np.sum((m1 - m2).ravel() ** 2 / m1.shape[0])

    return np.sqrt(rmsd)'''


def pdb2pdbqt(inp, out, keep_polarH=True):
    """Format convert with open babel"""

    if keep_polarH:
        mode = "AddPolarH"
    else:
        mode = "general"
    print(mode)
    babel_converter(inp, out, "obabel")

    return None


def docking_idock(rec, lig, out, ref_lig="ref.mol2",
                  verbose=True, exe='idock', cpu=12,
                  env_variables="export xxx=1 &&", 
                  config="vina.config", pocket_size=[15., 15., 15.]):
    # if prepare ligand

    #if not os.path.exists(lig + ".pdb"):
    #    builder.babel_converter(lig, lig + ".pdb")
    temp_dpath = "/opt/ml/disk/temp_docking"
    os.makedirs(temp_dpath, exist_ok=True)

    ligprep = LigandPrepare(lig)
    if not os.path.exists(lig + "_mgltools.pdbqt"):
        ligand_file_path = os.path.join(temp_dpath, os.path.basename(lig + "_mgltools.pdbqt"))
        #print(lig + "_mgltools.pdbqt")
        if not os.path.exists(ligand_file_path):
            ligprep.prepare_ligand(ligand_file_path)
    else:
        ligand_file_path = lig + "_mgltools.pdbqt"

    # define the pocket center
    if os.path.exists(ref_lig):
        print("INFO: ligand ref format", ref_lig[-4:])
        if ref_lig[-4:] == "mol2":
            #print(ref_lig[-4:])
            _xyz = _get_xyz_mol2(ref_lig)
        else:
           _xyz = _get_xyz_pdb(ref_lig)
    else:
        _xyz = _get_xyz_pdb(rec)
    xyz_c = np.mean(np.array(_xyz), axis=0)

    recprep = ReceptorPrepare(rec)
    # xyz_c = recprep.pocket_center(args.l)
    #if verbose:
    print("Binding Pocket: ", xyz_c)

    if not os.path.exists(rec + "_mgltools.pdbqt"):
        try:
            recprep.receptor_addH(rec + "_mgltools.pdbqt", method='mgltools')
        except:
            print("INFO: processing %s error ..." % rec)

    if not os.path.exists(rec + "_mgltools.pdbqt"):
        pdb2pdbqt(rec, rec + "_obabel.pdbqt", keep_polarH=True)
        shutil.copy(rec + "_obabel.pdbqt", rec + "_mgltools.pdbqt")

    if not os.path.exists(rec + "_mgltools.pdbqt"):
        print("INFO: processing %s error ..." % rec)
        sys.exit(0)

    # define the vina docking executable
    if os.path.exists(exe.split()[0]):
        _exe = exe
    else:
        _exe = os.path.join(os.path.dirname(__file__), "idock_backup")
        print("INFO: warning, cound not find {}, fall back to idock".format(exe.split()[0]))
    # configure file
    _configure = config #"vina_" + str(uuid.uuid4().hex)[:8] + ".config"

    docking = VinaDocking(_exe)
    docking.vina_config(rec + "_mgltools.pdbqt", ligand_file_path, out,
                        n_cpus=cpu, exhaustiveness=32,
                        center=xyz_c, boxsize=pocket_size,
                        logfile="log_vina.log", n_modes=20, config=_configure)
    docking.run_docking(env_variables)
    print("INFO: complete docking")


def docking_idock_tfoldms(rec, lig, out, ref_lig="ref.mol2",
                          verbose=False, exe='idock', cpu=12,
                          pocket_size=[15., 15., 15.],
                          temp_dpath="temp",
                          env_variables="MKL_THREADING_LAYER=GNU ", 
                          config="vina.config", 
                          user="zhengliangzhen", 
                          ceph="realbig", 
                          log_dpath="/opt/ml/disk/logs", 
                          image="docker.oa.com/seven-developer/tencent_threading:latest"):

    temp_dpath = os.path.join("/opt/ml/disk/docking_" + temp_dpath, os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(out)))))
    os.makedirs(temp_dpath, exist_ok=True)

    if os.path.exists(os.path.join(out, "log.csv")):
        print("INFO: find output file {}, skip ......".format(os.path.join(out, "log.csv")))
        return None

    ligprep = LigandPrepare(lig)
    if not os.path.exists(lig + "_mgltools.pdbqt"):
        #print(lig + "_mgltools.pdbqt")
        ligand_file_path = os.path.join(temp_dpath, os.path.basename(lig) + "_mgltools.pdbqt")
        #print(lig + "_mgltools.pdbqt")
        if not os.path.exists(ligand_file_path):
            ligprep.prepare_ligand(ligand_file_path)
    else:
        ligand_file_path = lig + "_mgltools.pdbqt"

    # define the pocket center
    if os.path.exists(ref_lig):
        if ref_lig[-4:] == "mol2":
            _xyz = _get_xyz_mol2(ref_lig)
        else:
            _xyz = _get_xyz_pdb(ref_lig)
    else:
        _xyz = _get_xyz_pdb(rec)
    xyz_c = np.mean(np.array(_xyz), axis=0)

    recprep = ReceptorPrepare(rec)
    # xyz_c = recprep.pocket_center(args.l)
    if verbose:
        print("Binding Pocket: ", xyz_c)

    if not os.path.exists(rec + "_mgltools.pdbqt"):
        try:
            recprep.receptor_addH(rec + "_mgltools.pdbqt", method='mgltools')
        except:
            print("INFO: processing %s error ..." % rec)

    if not os.path.exists(rec + "_mgltools.pdbqt"):
        pdb2pdbqt(rec, rec + "_obabel.pdbqt", keep_polarH=True)
        shutil.copy(rec + "_obabel.pdbqt", rec + "_mgltools.pdbqt")

    if not os.path.exists(rec + "_mgltools.pdbqt"):
        print("INFO: processing %s error ..." % rec)
        sys.exit(0)

    # define the vina docking executable
    if os.path.exists(exe.split()[0]):
        _exe = exe
    else:
        _exe = os.path.join(os.path.dirname(__file__), "idock_backup")
        print("INFO: warning, cound not find {}, fall back to idock".format(exe.split()[0]))

    # configure file
    _configure = config #"vina_" + str(uuid.uuid4().hex)[:8] + ".config"

    docking = VinaDocking(_exe)
    docking.vina_config(rec + "_mgltools.pdbqt", ligand_file_path, out,
                        n_cpus=cpu, exhaustiveness=32,
                        center=xyz_c, boxsize=pocket_size,
                        logfile="log_vina.log", n_modes=20, config=_configure)

    cmd = docking.run_docking(env_variables, mode='tfoldms')
    if len(cmd) and not os.path.exists(os.path.join(out, "log.csv")):
        if "idock" in _exe:
            tagid="idock" #os.path.basename(rec)[:4]
        elif "pyvina" in _exe:
            tagid = "pyvina"
        else:
            tagid = "docking"

        print("Submiting task: \n========================")
        print(cmd[0])
        print("======================")
     
        tfold_cli_submit(cmd[0], "idock", tagid, cpu=cpu, memory=(cpu*2), user=user, 
                         ceph=ceph, image=image, log_path=log_dpath)
    else:
        print("INFO: find output file {}, skip ......".format(os.path.join(out, "log.csv")))

    #print("INFO: complete docking")


def arguments():
    d = """
    Perform molecular docking using AutoDock Vina.

    Input protein and ligand structures with pdb, pdbqt and mol2 file formats.

    Vina only accepts .pdbqt files, thus the input coordination files would be 
    converted into pdbqt format. Only polar hydrogen atoms would be kept. 


    Examples:

    python vina_app.py -rec receptor.pdb -lig ligand.mol2 -out output_vina.pdbqt


    """

    parser = argparse.ArgumentParser(description=d)

    parser.add_argument("-p", type=str, default="receptor.pdb",
                        help="Input. Default is receptor.pdb. \n"
                             "The input receptor conformation.")
    parser.add_argument("-l", type=str, default="ligand.mol2",
                        help="Input. Default is ligand.mol2 . "
                             "The input ligand conformation.")
    parser.add_argument("-o", type=str, default="output_",
                        help="Output. Optional. Default is output_ \n"
                             "The prefix of the output")
    parser.add_argument("-c", type=str, default="ligand.mol2",
                        help="Get the pocket center of coordinates in ligands file.")
    parser.add_argument("-v", default=1, type=int,
                        help="Verbose.")
    parser.add_argument("-exe", default="vina", type=str,
                        help="Autodock Vina (semia/idock) executable path. ")
    parser.add_argument("-cpu", default=4, type=int,
                        help="Number of CPUs per task.")
    parser.add_argument("-bs", type=float, default=15,
                        help="Box size. Default 15 Angstrom.")
    parser.add_argument("-e", default=32, type=int,
                        help="Exhaustiveness value. Default 32.")
    parser.add_argument("-n", default=20, type=int,
                        help="Number of modes in output.")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = arguments()
    docking_idock(args.p, args.l, args.o, args.c, exe=args.exe, cpu=args.cpu)
