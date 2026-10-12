STIMTEC: this is the single fracture model approach (rectangle fracture)

revison 1

if you are the STIMTEC simulation folders
$cd target-directory
> BH10_20180716_40.6
> BH10_20180716_51.6

1) create a directory results: mkdir results
the subdirectory should be:
> results
> simulation
> workflow

2) go into directory simulation and run the ogs simulation, e.g.: ogs stimtec.prj

3) go into directory workflow and the python-function: python compare_sim_meas_en.py
it will show all available results (pvd files) in folder results

4) may also select a specific result with: python compare_sim_meas_en.py ../results/stimtec_bc_right_kS_t.pvd



claude instructions

11.10.2026: STIMTEC
(OK1)
I need a need Python program to compare the actual simulation (in subfolder results) with the measured pressure data (in bh10_40_gemessen_digitalisiert.csv).

(OK2)
Please show the results as well (e.g. plt.show()).

