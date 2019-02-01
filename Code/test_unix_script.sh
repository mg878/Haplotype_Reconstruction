#s=1
#t=0
#cd /rds/user/mg878/hpc-work/
#sbatch submit1.sandybridge /rds/user/mg878/hpc-work/Transmission_Project_050318/Transmission_ProjectNeutral_8genes_5loci_C1000000_simple_maxParams5_110118/Nt_50/Seed_$s/SimulatedData_Mahan_Gene_$t.dat
#cd ../ for C=200, Nt+50 folder is switched with Nt+1 folder

#for s in `seq 1 100`; do
#   for t in `seq 0 7`; do
#        cd /rds/user/mg878/hpc-work/Simulations/SimulatedData_step1_N_100_C_1000000/seed_$s/seeds/test_$t
#         cd /rds/user/mg878/hpc-work/HIV/length_50/
#         sbatch submit1.sandybridge /rds/user/mg878/hpc-work/HIV/length_50/Multi_locus_trajectories.out outcome_100.txt
#		 sbatch submit1.sandybridge /rds/user/mg878/hpc-work/Transmission_Project_050318/Transmission_ProjectNeutral_8genes_5loci_C1000000_simple_maxParams5_110118/Nt_100/Seed_$s/SimulatedData_Mahan_Gene_$t.dat
#		cd ../
#	done
#done 

#n=100
for s in `seq 1 43`; do
#   t=7
#   for t in `seq 0 7`; do
		cd /rds/user/mg878/hpc-work/Transmission_data_highN_noends_step1/Transmission$s/seeds/test_0
#        ./run_known_haps /rds/user/mg878/hpc-work/Transmission_Project_050318/Transmission_ProjectNeutral_8genes_5loci_C1000000_simple_maxParams5_110118/Nt_$n/Seed_$s/SimulatedData_Mahan_Gene_$t.dat /rds/user/mg878/hpc-work/Transmission_Project_050318/Transmission_ProjectNeutral_8genes_5loci_C1000000_simple_maxParams5_110118/Nt_$n/Seed_$s/SimulatedData_Gene_$t.dat /rds/user/mg878/hpc-work/Simulations/SimulatedData_step1_N_5_C_1000000/seed_$s/seeds/test_$t outcome_1.txt
#        echo $s $t
#		./run_known_haps /rds/user/mg878/hpc-work/Transmission_Project_050318/Transmission_ProjectNeutral_8genes_5loci_C200_simple_maxParams5_110118/Nt_$n/Seed_$s/SimulatedData_Gene_$t.dat /rds/user/mg878/hpc-work/Simulations/SimulatedData_step1_N_100_C_200/seed_$s/seeds/test_$t outcome_1.txt
#        ./run_good_bad /rds/user/mg878/hpc-work/Simulations/SimulatedData_step1_N_10_C_1000000/seed_$s/seeds/test_$t/ outcome_1_1.txt outcome_1.txt
		 sbatch submit1.sandybridge /rds/user/mg878/hpc-work/Transmission_data_highN_noends_step1/Transmission$s/HA/Multi_locus_trajectories.out
		cd ../
	done
#done 
#6_1 8_1 (2 no problem) and 17_6 from 8 to 16 
#for p in `seq 58 100`; do
#  echo $p
# for s in `seq 1001 1100`; do
#   echo $s
#   for t in `seq 0 7`; do
#   t=2
#   ./run_code /rds/user/mg878/hpc-work/Michigan_data_step2/Transmission$p/Seed_$s/SimulatedData_Mahan_Gene_$t.dat /rds/user/mg878/hpc-work/Michigan_data_step1/Transmission$p/seeds/test_$t/outcome_1.txt /rds/user/mg878/hpc-work/Michigan_data_step2/Transmission$p/Seed_$s test_$t.txt
#  ./run_test_2 /rds/user/mg878/hpc-work/Simulations/SimulatedData_step2_N_50_C_1000000/Seed_$p/Seed_$s/SimulatedData_Mahan_Gene_$t.dat /rds/user/mg878/hpc-work/Simulations/SimulatedData_step1_N_50_C_1000000/seed_$p/seeds/test_$t/outcome_1.txt /rds/user/mg878/hpc-work/Simulations/SimulatedData_step2_N_50_C_1000000/Seed_$p/Seed_$s test_$t.txt
#   done
# done 
#done
#for s in `seq 600 800`; do
#  scancel 3778$s
#done 