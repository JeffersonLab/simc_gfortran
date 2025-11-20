Objective:
There's slight difference how simc and hcana computes physics variables. To ensure faithful comparison between Data and Simc, recon_hcana computes simc generated variables in hcana style. 

Acknowledgement: This is based on prior work of Carlos Yero and Richard Trotta.

Procedure:
recon_hcana program grabs simc generated <filename_stem>.root file, copies branches/variables that don't require any reconstruction, reconstruct branches/variables the way how hcana computes, put all branches to tree, and outputs a new .root file with name under worksim/ as "recon_hcana_<filename_stem>.root".
You can either 
1. run simc usual way and then run recon_hcana.C; OR
2. run the shell script run_simc_recon.sh to run simc and recon_hcana.C all together. 
For case (2), you must run the code from your simc_gfortran directory, bcz, SIMC_DIR="$(pwd)"
If you find run_simc_recon.sh non-executable, make it executable by running
chmod +x run_simc_recon.sh

run commands:
1. running recon_hcana.C separately:

root -l -q -b 'recon_hcana.C+("filename_stem","reaction", "hadron_type", Earm_HMS)'

example:
root -l -b -q 'recon_hcana.C+("heep_rsidis_4pass_elastic1","heep", "mk", kFALSE)'
root -l -b -q 'recon_hcana.C+("coin_7p87deg_3p632gev_hyd_rsidis","sidis", "mpi", kTRUE)'
root -l -b -q 'recon_hcana.C+("rho_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2","rho")'
root -l -b -q 'recon_hcana.C+("delta_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2","delta")'
root -l -b -q 'recon_hcana.C+("exclusive_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2","exclusive")'


2. running shell script:
./run_simc_recon.sh <stem> <reaction> <hadron_type> <Earm_HMS>
where:
<stem>        = base name of your SIMC input/output (without .inp)
<reaction>    = "heep", "sidis", "rho", "delta", "exclusive".
<hadron_type> = "mpi" (default) or "mk". For heep, proton mass is hard-coded
<Earm_HMS>    = 1 (default: electron in HMS) or 0 (electron in SHMS)

example:
./run_simc_recon.sh heep_12p47deg_7p07gev_hyd_rsidis heep mk 0
./run_simc_recon.sh coin_7p87deg_3p632gev_hyd_rsidis sidis mpi 1
./run_simc_recon.sh rho_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2 rho mpi 1
./run_simc_recon.sh delta_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2 delta mpi 1
./run_simc_recon.sh exclusive_kin37_q24p4_x0p44_z0p67_thpi12p9_pip_lh2 sidis mpi 1

Note:
1. recon_hcana.h declares all necessary class, functions, variables; recon_hcana.C computes.
2. Need to have <filname_stem>.{hist,root} under {outfiles/,worksim/} respectively.
3. Output tree will have both simc and hcana styled computation for some physics variables. Example: Q2, Q2_recon, where Q2 is computed in simc, Q2_recon computed in hcana style.
4. Output tree will contain same set of branches for all reaction_type except "heep". If some branches are not found in original simc generated file, that branch will be valued ErrorVal = -999.0
5. Don't bother about progress bar not completing 100%. It's computing (i%1000) where i is nentries with nentries = tree->GetEntries();
6. reaction_types = {"heep", "sidis", "rho", "delta", "exclusive"}
7. argument hadron_type selects final hadrons mass, can be mass of {kaon,pion}. If not passed this argument, mpi (pion mass) will be seleted.
8. Argument Earm_HMS is kTRUE if electron arm is HMS. For heep, kFALSE. If not passed, default kTRUE.
9. Added Ngen, normfac, fWeight (= Weight*normfac/Ngen) in output tree.
10. normfac is generally a big number. new TBrowser may show empty plot. You need to use specific range to see normfac. 

Contact: if you see any problem/bug, please contact rparvez@jlab.org
