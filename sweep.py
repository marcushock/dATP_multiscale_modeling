import os
import subprocess
import time
import numpy as np 

# Make sure to also change baseline ATP params in force_pCa_curve.cu (if you want to change them for ATP)!
# Based on the above line, it appears that these changes down below will influence the dATP parameters when we change it. 
# When the percent_cycle is set to 1, then it usees the new parameters completely. 
# When it is set to something less, than it uses that fraction for dATP. 
# So if you just want to change paramters, not necessarily comapring them to ATP, you can leave percent_cycle at 1
# Then you can adjust as you see fit. 

def rename_undated(results_directory= 'MCMC_simulation_results/', extra_append = ""):
	file_names = os.listdir(results_directory)
	time_str = time.strftime("%y%m%d-%H%M_")+extra_append
	for name in file_names:
		if "States_" == name[0:7] or "Force_" == name[0:6]:
			new_filename = time_str+fix_trailing_zeros(name)
			print("Renaming '{}' \nto \n'{}'\n".format(name, new_filename))
			original_path = results_directory+'/'+name
			new_path = results_directory + '/' + new_filename
			try:
				os.rename(original_path, new_path)
			except:
				print('Could not rename. ')
	return 

def fix_trailing_zeros(filename_string):
    new_name = ''
    for filename_piece in filename_string.strip('.csv').split(' '):
        try:
            new_value = float(filename_piece)
            new_name += str(new_value)
            new_name += ' '
        except:
            new_name += filename_piece
            new_name += ' '
    if new_name[-1] == ' ':
        new_name = new_name[0:-1]
    new_name += '.csv'
    return new_name


binName= "MCMC_CUDA_10States"
defaults = {
    "gamma_B": 45,
    "gamma_M": 21,
    "mu_M": 2,
    "k2_plus_ref": 0.0025, #kf+
    "k3_plus": 0.05, #kp+
    "k4_plus_ref": 0.135, #kg+
    "kB_plus_ref": 13, 
    "kB_minus_ref": 0.1, 
    "lambda": 0,
    "kCa_plus_ref": 0.09, 
    "kCa_minus_ref": 0.57,
    "percent_dATP": 0,
    "k_force": 0.2, #krecruit
    "k_plus_SR_ref": 16, #km+
    "k_minus_SR_ref": 15, #km-
}

def main():
# Change parameters here for dATP simulation
	k2_cycle = [0.0025] #kf+ [0.0025, ]   # Original [0.0025] #[0.00478] #kf+
	k3_cycle = [0.05] #[0.08] #kp+
	k4_cycle = [0.135] #[0.23] #kg+ 
	percent_cycle = [1] #percent dATP [0.05,0.1,0.15,0.25,0.5,0.75]  # Original [1] #[0.01] #percent dATP 
	k_force_cycle = [0.207] #krecruit default = 0.2  # Original [0.2] #[779] #krecruit 
	k_plus_SR_ref_cycle = [5e6] #km+
	k_minus_SR_ref_cycle = [15] #km-
	for k2_plus_ref in k2_cycle:
		for k3_plus in k3_cycle:
			for k4_plus_ref in k4_cycle:
				for k_force in k_force_cycle:
					for k_plus_SR_ref in k_plus_SR_ref_cycle:
						for k_minus_SR_ref in k_minus_SR_ref_cycle:
							for percent_dATP in percent_cycle:
								# This is skipping the definition of "parameters", by just manually entering a dictionary with 
								# the needed values for the analysis. 
								# Note that not all of the parameters necessarily need to be defined in the input argument 
								# because they will be fille
								callProgram({
									"k2_plus_ref": k2_plus_ref,
									"k3_plus": k3_plus,
									"k4_plus_ref": k4_plus_ref,
									"percent_dATP": percent_dATP,
									"k_force": k_force,
									"k_plus_SR_ref": k_plus_SR_ref,
									"k_minus_SR_ref": k_minus_SR_ref,
									},
									folder="sweep/sweep2"
									)
	print("Done first cycle")
	rename_undated('MCMC_simulation_results', 'check_following_') ### Can change to add to file name appending ###





def callProgram(parameters, folder="sweep"):
	for key, value in defaults.items():
		if key not in parameters:
			parameters[key] = value

	args = [
		"bin/" + str(binName),
		"expData/Force_pCa_Optmz_morePts.csv", #Force_pCa_3pts_2 #Force_pCa_Optmz_morePts
		str(parameters["gamma_B"]),
		str(parameters["gamma_M"]),
		str(parameters["mu_M"]),
		str(parameters["k2_plus_ref"]),
		str(parameters["k3_plus"]),
		str(parameters["k4_plus_ref"]),
		str(parameters["kB_plus_ref"]),
		str(parameters["kB_minus_ref"]),
		str(parameters["lambda"]),
		str(parameters["kCa_plus_ref"]),
		str(parameters["kCa_minus_ref"]),
		str(parameters["percent_dATP"]),
		str(parameters["k_force"]),
		str(parameters["k_plus_SR_ref"]),
		str(parameters["k_minus_SR_ref"]),
		 ]
	try:
		os.makedirs(folder)
	except OSError:
		pass
	fname = folder+"/"
	for key, value in parameters.items():
		fname += str(key) + " " + str(value) + " "
	fname = fname.strip()
	with open(fname+".out", 'w') as out:
		with open(fname+".err", 'w') as err:
			print("Running " + str(args))
			time_start = time.time()
			subprocess.call(args, stdout=out, stderr=err, shell=False)
			print("Time (seconds): ", time.time() - time_start)

if __name__ == "__main__":
	main()
