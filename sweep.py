import os
import subprocess

# Make sure to also change baseline ATP params in force_pCa_curve.cu (if you want to change them for ATP)!

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
	k2_cycle = [0.0025] #[0.00478] #kf+
	k3_cycle = [0.05] #[0.08] #kp+
	k4_cycle = [0.135] #[0.23] #kg+ 
	percent_cycle = [1] #[0.01] #percent dATP 
	k_force_cycle = [0.2] #[779] #krecruit 
	k_plus_SR_ref_cycle = [16] #km+
	k_minus_SR_ref_cycle = [15] #km-
	for k2_plus_ref in k2_cycle:
		for k3_plus in k3_cycle:
			for k4_plus_ref in k4_cycle:
				for k_force in k_force_cycle:
					for k_plus_SR_ref in k_plus_SR_ref_cycle:
						for k_minus_SR_ref in k_minus_SR_ref_cycle:
							for percent_dATP in percent_cycle:
						   		callProgram({
						     		"k2_plus_ref": k2_plus_ref,
						     		"k3_plus": k3_plus,
								"k4_plus_ref": k4_plus_ref,
						      		"percent_dATP": percent_dATP,
						      		"k_force": k_force,
						      		"k_plus_SR_ref": k_plus_SR_ref,
						      		"k_minus_SR_ref": k_minus_SR_ref
						       		},
						     		folder="sweep/sweep2"
						      		)
	print("Done first cycle")


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
			subprocess.call(args, stdout=out, stderr=err, shell=False)

if __name__ == "__main__":
	main()
