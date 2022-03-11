def run(step_number, steps, dt, temp, press=None, gen_vel=False):
    cont = "no" if gen_vel else "yes"
    mdp = [
        "integrator              = md        ; leap-frog integrator",
        f"nsteps                  = {steps}   ; {steps*dt/1000} ps",
        f"dt                      = {dt/1000}      ; {dt} fs",
        "; Output control",
        "nstxout-compressed      = 10000     ; save compressed coordinates",
        "nstxout                 = 0         ; save coordinates",
        "nstvout                 = 0         ; save velocities",
        "nstenergy               = 1000      ; save energies",
        "nstlog                  = 1000      ; update log file",
        "; Bond parameters",
        f"continuation            = {cont}    ; set velocities before first step",
        "constraint_algorithm    = lincs     ; holonomic constraints",
        "constraints             = h-bonds   ; bonds involving H are constrained",
        "lincs_iter              = 1         ; accuracy of LINCS",
        "lincs_order             = 4         ; also related to accuracy",
        "; Nonbonded settings",
        "cutoff-scheme           = Verlet    ; Buffered neighbor searching",
        "ns_type                 = grid      ; search neighboring grid cells",
        "nstlist                 = 10        ; 10 fs, largely irrelevant with Verlet",
        "rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)",
        "rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)",
        "DispCorr                = EnerPres  ; account for cut-off vdW scheme",
        "; Electrostatics",
        "coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics",
        "pme_order               = 4         ; cubic interpolation",
        "fourierspacing          = 0.16      ; grid spacing for FFT",
        "; Temperature coupling is on",
        "tcoupl                  = V-rescale ; modified Berendsen thermostat",
        "tc-grps                 = System",
        "tau-t                   = 0.1       ; time constant, in ps",
        f"ref-t                   = {temp}    ; reference temperature, one for each group, in K",
        "; Periodic boundary conditions",
        "pbc                     = xyz       ; 3-D PBC",
    ]

    if press is not None:
        mdp += [
            "; Pressure coupling",
            "pcoupl                = Berendsen",
            "pcoupltype            = isotropic",
            "tau-p                 = 1.0",
            f"ref-p                 = {press}",
            "compressibility       = 4.5e-5",
        ]
    if gen_vel:
        mdp += [
            "; Velocity generation",
            "gen_vel                 = yes       ; assign velocities from Maxwell distribution",
            f"gen_temp                = {temp}       ; temperature for Maxwell distribution",
            "gen_seed                = -1        ; generate a random seed",
        ]
    else:
        mdp += [
            "; Velocity generation",
            "gen_vel                 = no       ; use velocities from previous run",
        ]

    with open(f"{step_number}.mdp", "w") as f:
        for line in mdp:
            f.write(f"{line}\n")


run(1, 25000, 2, 1000, gen_vel=True)
run(2, 25000, 2, 300)
run(3, 25000, 2, 300, press=1000)
run(4, 25000, 2, 1000)
run(5, 50000, 2, 300)
run(6, 25000, 2, 300, press=30000)
run(7, 25000, 2, 1000)
run(8, 50000, 2, 300)
run(9, 25000, 2, 300, press=50000)
run(10, 25000, 2, 1000)
run(11, 50000, 2, 300)
run(12, 2500, 2, 300, press=25000)
run(13, 2500, 2, 1000)
run(14, 5000, 2, 300)
run(15, 2500, 2, 300, press=5000)
run(16, 2500, 2, 1000)
run(17, 5000, 2, 300)
run(18, 2500, 2, 300, press=500)
run(19, 2500, 2, 1000)
run(20, 5000, 2, 300)
run(21, 400000, 2, 300, press=1)


slurm_cmds = [
    "#!/bin/bash",
    "#SBATCH -J 21-step",
    "#SBATCH -e run.%j.err",
    "#SBATCH -o run.%j.out",
    "#SBATCH -N 1",
    "#SBATCH -n 48",
    "#SBATCH -p skx-dev",
    "#SBATCH -t 2:00:00",
    "# #SBATCH --mail-user=thomas.mason1@monash.edu",
    "# #SBATCH --mail-type=all",
    "",
    "module load gromacs/2019.6",
]


# maxwarn 1 for charge
def generate_script(
    steps=21,
    ref_coord="../npt.gro",
    top="../topol.top",
    gmxexe="gmx",
    mdrun="ibrun mdrun_mpi",
    scheduler_cmds=slurm_cmds,
    maxwarn=1,
):
    with open("run.sh", "w") as f:
        if scheduler_cmds is not None:
            for line in scheduler_cmds:
                f.write(f"{line}\n")
        # assigning velocities here, so using coords only is sufficient
        # - no need to use restart
        f.write("# Step 1\n")
        f.write(
            f"{gmxexe} grompp -f 1.mdp -c {ref_coord} -p {top} -o 1.tpr -maxwarn {maxwarn}\n"
        )
        f.write(f"{mdrun} -s 1.tpr -deffnm 1 -g 1.log\n")
        for n in range(2, steps + 1):
            f.write(f"# Step {n}\n")
            # read restart from previous step
            f.write(
                f"{gmxexe} grompp -f {n}.mdp -c {n - 1}.gro -t {n - 1}.cpt -p {top} -o {n}.tpr -maxwarn {maxwarn}\n"
            )
            f.write(f"{mdrun} -s {n}.tpr -deffnm {n} -g {n}.log\n")


generate_script(ref_coord="../npt.gro", top="../topol.top", maxwarn=2)
