#pragma once
#include "MTwisterFunctions.h"
#include <time.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include "ChiSquare.h"


class epistatic_variance_list
{
public:
	double ep_var_0by0_on_0;
	double ep_var_1by1_on_1;
	double ep_var_0by1_on_0;
	double ep_var_0by1_on_1;
	double ep_var_0by0_on_1;
	double ep_var_1by1_on_0;

	double pleio_ep_var_0by0_on_0;
	double pleio_ep_var_1by1_on_1;
	double pleio_ep_var_0by1_on_0;
	double pleio_ep_var_0by1_on_1;
	double pleio_ep_var_0by0_on_1;
	double pleio_ep_var_1by1_on_0;
};

class parameter_value_set
{
public:
	int p_reps, p_init_gens;
	double p_init_w00, p_init_w11, p_init_sel_corr, p_init_opt0, p_init_opt1;
	bool p_init_sex_lim;
	int p_gens, p_carry_cap, p_fecund, p_mate_enc;
	double p_pref_var;
	int p_loci_trt0, p_loci_trt1, p_loci_pleio; // these values indicate the numbers of qtls per chromosome
	double p_env_var0, p_env_var1, p_mut_var0, p_mut_var1, p_mut_corr, p_mut_rate;
	double p_exp_w00, p_exp_w11, p_exp_sel_corr, p_exp_opt0, p_exp_opt1;
	bool p_exp_sex_lim;
	double p_migration_rate;
	int p_number_chromosomes, p_markers_per_chromosome;
	double p_exp_recomb_rate, p_starting_QTL_stdev;
	double p_marker_mut_rate;
	int p_max_alleles_per_marker;
	int p_sample_size;
	double pop2_opt0, pop2_opt1, pop2_w00, pop2_w11, pop2_sel_corr;
	int p_all_freq_calc_interval;
	double p_fst_weighting_variance;
	int p_intervening_gens;

	bool p_permit_epistasis;
	epistatic_variance_list p_epistatic_vars;

	std::string file_name;
};

class epistatic_parameter_collection
{
public:
	std::vector<double> trait0loci_by_trait0loci_effect_on_trait0;
	std::vector<double> trait1loci_by_trait1loci_effect_on_trait1;
	std::vector<double> trait0loci_by_trait1loci_effect_on_trait0;
	std::vector<double> trait0loci_by_trait1loci_effect_on_trait1;
	std::vector<double> trait0loci_by_trait0loci_effect_on_trait1;
	std::vector<double> trait1loci_by_trait1loci_effect_on_trait0;

	std::vector<double> pleio_trait0_by_trait0_effect_on_trait0;
	std::vector<double> pleio_trait1_by_trait1_effect_on_trait1;
	std::vector<double> pleio_trait0_by_trait1_effect_on_trait0;
	std::vector<double> pleio_trait0_by_trait1_effect_on_trait1;
	std::vector<double> pleio_trait0_by_trait0_effect_on_trait1;
	std::vector<double> pleio_trait1_by_trait1_effect_on_trait0;

	void set_epistatic_effects(int n_chrom, int n_trt0_loci_per, int n_trt1_loci_per, int n_pleio_loci_per, epistatic_variance_list &ep_vars)
	{
		int i, number_needed;
		int n_trt0_loci, n_trt1_loci, n_pleio_loci;

		n_trt0_loci = n_chrom * n_trt0_loci_per;
		n_trt1_loci = n_chrom * n_trt1_loci_per;
		n_pleio_loci = n_chrom * n_pleio_loci_per;

		number_needed = (n_trt0_loci)*(n_trt0_loci - 1) / 2;
		for (i = 0; i < number_needed; i++)
			trait0loci_by_trait0loci_effect_on_trait0.push_back(randnorm(0, ep_vars.ep_var_0by0_on_0));

		number_needed = (n_trt1_loci)*(n_trt1_loci - 1) / 2;
		for (i = 0; i < number_needed; i++)
			trait1loci_by_trait1loci_effect_on_trait1.push_back(randnorm(0, ep_vars.ep_var_1by1_on_1));

		number_needed = (n_trt0_loci)*(n_trt1_loci);
		for (i = 0; i < number_needed; i++)
			trait0loci_by_trait1loci_effect_on_trait0.push_back(randnorm(0, ep_vars.ep_var_0by1_on_0));

		number_needed = (n_trt0_loci)*(n_trt1_loci);
		for (i = 0; i < number_needed; i++)
			trait0loci_by_trait1loci_effect_on_trait1.push_back(randnorm(0, ep_vars.ep_var_0by1_on_1));

		number_needed = (n_trt0_loci)*(n_trt0_loci - 1) / 2;
		for (i = 0; i < number_needed; i++)
			trait0loci_by_trait0loci_effect_on_trait1.push_back(randnorm(0, ep_vars.ep_var_0by0_on_1));

		number_needed = (n_trt1_loci)*(n_trt1_loci - 1) / 2;
		for (i = 0; i < number_needed; i++)
			trait1loci_by_trait1loci_effect_on_trait0.push_back(randnorm(0, ep_vars.ep_var_1by1_on_0));

		// set epistatic effects for pleiotropic loci

		number_needed = (n_pleio_loci)*(n_pleio_loci - 1) / 2;
		for (i = 0; i < number_needed; i++)
			pleio_trait0_by_trait0_effect_on_trait0.push_back(randnorm(0, ep_vars.pleio_ep_var_0by0_on_0));

		for (i = 0; i < number_needed; i++)
			pleio_trait1_by_trait1_effect_on_trait1.push_back(randnorm(0, ep_vars.pleio_ep_var_1by1_on_1));

		for (i = 0; i < number_needed; i++)
			pleio_trait0_by_trait1_effect_on_trait0.push_back(randnorm(0, ep_vars.pleio_ep_var_0by1_on_0));
		
		for (i = 0; i < number_needed; i++)
			pleio_trait0_by_trait1_effect_on_trait1.push_back(randnorm(0, ep_vars.pleio_ep_var_0by1_on_1));
		
		for (i = 0; i < number_needed; i++)
			pleio_trait0_by_trait0_effect_on_trait1.push_back(randnorm(0, ep_vars.pleio_ep_var_0by0_on_1));
		
		for (i = 0; i < number_needed; i++)
			pleio_trait1_by_trait1_effect_on_trait0.push_back(randnorm(0, ep_vars.pleio_ep_var_1by1_on_0));

	}
};

int parse_command_line_arguments(int arg_c, char* arg_v[], parameter_value_set &parm_set)
{
	// First set default values
	parm_set.p_reps = 1;

	// Initial Generations Parameters
	parm_set.p_init_gens = 10000;
	parm_set.p_init_w00 = 49;
	parm_set.p_init_w11 = 49;
	parm_set.p_init_sel_corr = 0;
	parm_set.p_init_opt0 = 0;
	parm_set.p_init_opt1 = 0;
	parm_set.p_init_sex_lim = false;
	parm_set.p_starting_QTL_stdev = 0.05;
	parm_set.p_intervening_gens = 2000;

	// Demographic Parameters
	parm_set.p_gens = 2000;
	parm_set.p_carry_cap = 500;
	parm_set.p_fecund = 4;
	parm_set.p_sample_size = 100;

	// Mating Parameters
	parm_set.p_mate_enc = 50;
	parm_set.p_pref_var = 0;

	// Genetic Parameters
	parm_set.p_number_chromosomes = 4;
	parm_set.p_loci_trt0 = 1; // loci per chromosome
	parm_set.p_loci_trt1 = 1;
	parm_set.p_loci_pleio = 0;
	parm_set.p_env_var0 = 1;
	parm_set.p_env_var1 = 1;
	parm_set.p_exp_recomb_rate = 0.5;

	// Mutation Parameters for QTL
	parm_set.p_mut_var0 = 0.2;
	parm_set.p_mut_var1 = 0.2;
	parm_set.p_mut_corr = 0;
	parm_set.p_mut_rate = 0.0002; 

	// Selection Parameters
	parm_set.p_exp_w00 = 49;
	parm_set.p_exp_w11 = 49;
	parm_set.p_exp_sel_corr = 0;
	parm_set.p_exp_opt0 = 4;
	parm_set.p_exp_opt1 = 0;
	parm_set.p_exp_sex_lim = false;

	// Migration Rate
	parm_set.p_migration_rate = 0.016; 

	// Marker Loci
	parm_set.p_markers_per_chromosome = 2000;
	parm_set.p_marker_mut_rate = 0.0002;
	parm_set.p_max_alleles_per_marker = 4;
	parm_set.p_fst_weighting_variance = 500;

	// Population two parameters
	parm_set.pop2_opt0 = -4;
	parm_set.pop2_opt1 = 0;
	parm_set.pop2_w00 = 49;
	parm_set.pop2_w11 = 49;
	parm_set.pop2_sel_corr = 0;

	// Epistatic parameters
	parm_set.p_permit_epistasis = false;
	parm_set.p_epistatic_vars.ep_var_0by0_on_0 = 0;
	parm_set.p_epistatic_vars.ep_var_1by1_on_1 = 0;
	parm_set.p_epistatic_vars.ep_var_0by1_on_0 = 0;
	parm_set.p_epistatic_vars.ep_var_0by1_on_1 = 0;
	parm_set.p_epistatic_vars.ep_var_1by1_on_0 = 0;
	parm_set.p_epistatic_vars.ep_var_0by0_on_1 = 0;
	parm_set.p_epistatic_vars.pleio_ep_var_0by0_on_0 = 0;
	parm_set.p_epistatic_vars.pleio_ep_var_1by1_on_1 = 0;
	parm_set.p_epistatic_vars.pleio_ep_var_0by1_on_0 = 0;
	parm_set.p_epistatic_vars.pleio_ep_var_0by1_on_1 = 0;
	parm_set.p_epistatic_vars.pleio_ep_var_1by1_on_0 = 0;
	parm_set.p_epistatic_vars.pleio_ep_var_0by0_on_1 = 0;

	
	parm_set.file_name = "outfile";
	parm_set.p_all_freq_calc_interval = 2000;

	if (arg_c < 2) // There are no arguments so defaults are in use
		return 0;

	if (arg_v[1] == "-h" || arg_v[1] == "--help") // user wants help
		return 1;

	int i;
	std::string tstr1, tstr2;
	for (i = 1; i < arg_c - 1; i++)
	{
		tstr1 = arg_v[i];
		tstr2 = arg_v[i + 1];
		if (tstr1 == "--reps")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_reps;
		}
		if (tstr1 == "--init_gens")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_init_gens;
		}
		if (tstr1 == "--init_w00")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_init_w00;
		}
		if (tstr1 == "--init_w11")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_init_w11;
		}
		if (tstr1 == "--init_sel_corr")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_init_sel_corr;
		}
		if (tstr1 == "--init_opt0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_init_opt0;
		}
		if (tstr1 == "--init_opt1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_init_opt1;
		}
		if (tstr1 == "--init_sex_lim")
		{
			if (tstr2[0] == 't' || tstr2[0] == 'T')
				parm_set.p_init_sex_lim = true;
		}
		if (tstr1 == "--gens")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_gens;
		}
		if (tstr1 == "--carry_cap")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_carry_cap;
		}
		if (tstr1 == "--fecund")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_fecund;
		}
		if (tstr1 == "--mate_enc")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_mate_enc;
		}
		if (tstr1 == "--pref_var")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_pref_var;
		}
		if (tstr1 == "--loci_trt0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_loci_trt0;
		}
		if (tstr1 == "--loci_trt1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_loci_trt1;
		}
		if (tstr1 == "--loci_pleio")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_loci_pleio;
		}
		if (tstr1 == "--env_var0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_env_var0;
		}
		if (tstr1 == "--env_var1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_env_var1;
		}
		if (tstr1 == "--mut_var0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_mut_var0;
		}
		if (tstr1 == "--mut_var1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_mut_var1;
		}
		if (tstr1 == "--mut_corr")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_mut_corr;
		}
		if (tstr1 == "--mut_rate")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_mut_rate;
		}
		if (tstr1 == "--exp_w00")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_exp_w00;
		}
		if (tstr1 == "--exp_w11")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_exp_w11;
		}
		if (tstr1 == "--exp_sel_corr")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_exp_sel_corr;
		}
		if (tstr1 == "--exp_opt0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_exp_opt0;
		}
		if (tstr1 == "--exp_opt1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_exp_opt1;
		}
		if (tstr1 == "--exp_sex_lim")
		{
			if (tstr2[0] == 't' || tstr2[0] == 'T')
				parm_set.p_exp_sex_lim = true;
		}
		if (tstr1 == "--filename")
		{
			parm_set.file_name = tstr2;
		}
		if (tstr1 == "--migr_r")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_migration_rate;
		}
		if (tstr1 == "--n_markers")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_markers_per_chromosome;
		}
		if (tstr1 == "--n_chrom")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_number_chromosomes;
		}
		if (tstr1 == "--start_qtl_sd")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_starting_QTL_stdev;
		}
		if (tstr1 == "--marker_mut_rate")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_marker_mut_rate;
		}
		if (tstr1 == "--marker_max_alleles")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_max_alleles_per_marker;
		}
		if (tstr1 == "--recomb_rate")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_exp_recomb_rate;
		}
		if (tstr1 == "--sample_size")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_sample_size;
		}
		if (tstr1 == "--pop2_opt0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.pop2_opt0;
		}
		if (tstr1 == "--pop2_opt1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.pop2_opt1;
		}
		if (tstr1 == "--pop2_w00")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.pop2_w00;
		}
		if (tstr1 == "--pop2_w11")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.pop2_w11;
		}
		if (tstr1 == "--pop2_selcorr")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.pop2_sel_corr;
		}
		if (tstr1 == "--all_freq_calc_int")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_all_freq_calc_interval;
		}

		if (tstr1 == "--use_epistasis")
		{
			if (tstr2[0] == 'y' || tstr2[0] == 'Y')
				parm_set.p_permit_epistasis = true;
			else
				parm_set.p_permit_epistasis = false;
		}

		if (tstr1 == "--non_pleio_ep_var_0by0_on_0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.ep_var_0by0_on_0;
		}
		if (tstr1 == "--non_pleio_ep_var_1by1_on_1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.ep_var_1by1_on_1;
		}
		if (tstr1 == "--non_pleio_ep_var_0by1_on_0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.ep_var_0by1_on_0;
		}
		if (tstr1 == "--non_pleio_ep_var_0by1_on_1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.ep_var_0by1_on_1;
		}
		if (tstr1 == "--non_pleio_ep_var_1by1_on_0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.ep_var_1by1_on_0;
		}
		if (tstr1 == "--non_pleio_ep_var_0by0_on_1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.ep_var_0by0_on_1;
		}
		
		if (tstr1 == "--pleio_ep_var_0by0_on_0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.pleio_ep_var_0by0_on_0;
		}
		if (tstr1 == "--pleio_ep_var_1by1_on_1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.pleio_ep_var_1by1_on_1;
		}
		if (tstr1 == "--pleio_ep_var_0by1_on_0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.pleio_ep_var_0by1_on_0;
		}
		if (tstr1 == "--pleio_ep_var_0by1_on_1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.pleio_ep_var_0by1_on_1;
		}
		if (tstr1 == "--pleio_ep_var_1by1_on_0")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.pleio_ep_var_1by1_on_0;
		}
		if (tstr1 == "--pleio_ep_var_0by0_on_1")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_epistatic_vars.pleio_ep_var_0by0_on_1;
		}
		if (tstr1 == "--fst_weight_var")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_fst_weighting_variance;
		}
		if (tstr1 == "--intervening_gens")
		{
			std::stringstream ss(tstr2);
			ss >> parm_set.p_intervening_gens;
		}
		
	}

	return 2; // means parameters have changed
}

void convert_string_to_char(std::string my_string, char *my_char)
{
	size_t i;
	size_t str_length;

	str_length = my_string.size();

	if (str_length > 255)
		str_length = 255;

	for (i = 0; i < str_length; i++)
	{
		my_char[i] = my_string[i];
	}

	my_char[i] = '\0';
}

class marker_allele_frequencies
{
public:
	std::vector<double> freq;
};

class qtl_allele_frequencies
{
public:
	std::vector<double> trt0_effect;
	std::vector<double> trt1_effect;
	std::vector<double> freq;
};

class peak_data
{
public:
	int numberofpeaks;
	std::vector<int> peaklocation;
	std::vector<double> peaksmoothedFst;
	std::vector<double> highestFstonplateau;
	std::vector<int> highestFstIndex;
	std::vector<bool> sigFDR;
	std::vector<bool> sig80CI;
	std::vector<bool> sig90CI;
	std::vector<bool> sig95CI;
	std::vector<bool> sig98CI;
	std::vector<bool> sig99CI;
	std::vector<bool> sig05FSTprime;
	std::vector<bool> sig01FSTprime;
	std::vector<bool> sigFDRFSTprime;
	std::vector<int> nearestQTL_0;
	std::vector<int> distancetoQTL_0;
	std::vector<int> nearestQTL_1;
	std::vector<int> distancetoQTL_1;
	std::vector<int> nearestQTL_P;
	std::vector<int> distancetoQTL_P;
};

class chromosomedata
{
public:
	std::vector<qtl_allele_frequencies> qtl_afs_trt0;
	std::vector<qtl_allele_frequencies> qtl_afs_trt1;
	std::vector<qtl_allele_frequencies> qtl_afs_pleio;
	std::vector<marker_allele_frequencies> marker_afs;
	std::vector<int> qtl_loc_trt0;
	std::vector<int> qtl_loc_trt1;
	std::vector<int> qtl_loc_pleio;

	std::vector<double> qtl_fst_values_trt0;
	std::vector<double> qtl_fst_values_trt1;
	std::vector<double> qtl_fst_values_pleio;
	std::vector<double> marker_fst_values;
	std::vector<double> smoothed_marker_fsts;
	std::vector<double> fst_chi_sqr_p;
	std::vector<bool> is_significant;
	std::vector<bool> smoothedFstpeak;
	std::vector<double> marker_heterozygosity;

	std::vector<double> qtl_var_trt0;
	std::vector<double> qtl_mean_trt0;
	std::vector<double> qtl_var_trt1;
	std::vector<double> qtl_mean_trt1;
	std::vector<double> qtl_var_pleio_0;
	std::vector<double> qtl_mean_pleio_0;
	std::vector<double> qtl_var_pleio_1;
	std::vector<double> qtl_mean_pleio_1;

	std::vector<double> marker_fst_prime_value;
	std::vector<double> fst_prime_p_value;
	std::vector<bool> fst_prime_is_significant_point05;
	std::vector<bool> fst_prime_is_significant_point01;
	std::vector<bool> fst_prime_is_significant_fdrpoint05;

	double SmoothedCritValue99, SmoothedCritValue98, SmoothedCritValue95;
	double SmoothedCritValue90, SmoothedCritValue80, FDRp;


	peak_data peakdata;

	int find_all_local_smoothed_fst_maxima(int halfregressioninterval)
	{
		int jjj, kkk;

		int NumberMarkerLociPerChromosome;
		NumberMarkerLociPerChromosome = static_cast<int>(smoothed_marker_fsts.size());
		if (NumberMarkerLociPerChromosome <= halfregressioninterval)
			return 0;

		double currentregression;
		double previousregression;
		double regressionproduct;
		double meanX, meanY, covarXY, varX;
		int start, stop;
		start = halfregressioninterval;
		stop = NumberMarkerLociPerChromosome - halfregressioninterval;
		double numbersummed;
		double tempdub;
		bool isamax;
		double lowestpeak, highestpeak;
		double highestpointinendsegment;
		int indexofhighestpoint;

		smoothedFstpeak.clear();
		smoothedFstpeak.resize(NumberMarkerLociPerChromosome);
	
		previousregression = -1000000;
		lowestpeak = 2;
		highestpeak = 0;
		for (jjj = start; jjj < stop; jjj++)
		{
			numbersummed = 0;
			meanY = 0;
			meanX = 0;
			for (kkk = jjj - halfregressioninterval; kkk < jjj + halfregressioninterval; kkk++)
			{
				meanX = meanX + kkk;
				meanY = meanY + smoothed_marker_fsts[kkk];
				numbersummed++;
			}
			meanX = meanX / numbersummed;
			meanY = meanY / numbersummed;

			covarXY = 0;
			varX = 0;
			for (kkk = jjj - halfregressioninterval; kkk < jjj + halfregressioninterval; kkk++)
			{
				tempdub = kkk;
				covarXY = covarXY + (tempdub - meanX)*(smoothed_marker_fsts[kkk] - meanY);
				varX = varX + (smoothed_marker_fsts[kkk] - meanY)*(smoothed_marker_fsts[kkk] - meanY);
			}
			covarXY = covarXY / (numbersummed - 1);
			varX = varX / (numbersummed - 1);

			currentregression = covarXY / varX;
			isamax = false;
			if (previousregression > -1000000)
			{
				regressionproduct = currentregression * previousregression;
				if (regressionproduct <= 0 && previousregression > 0)
					isamax = true;
			}

			if (isamax)
			{
				smoothedFstpeak[jjj - 1] = true;
				if (smoothed_marker_fsts[jjj - 1] < lowestpeak)
					lowestpeak = smoothed_marker_fsts[jjj - 1];
				if (smoothed_marker_fsts[jjj - 1] > highestpeak)
					highestpeak = smoothed_marker_fsts[jjj - 1];
			}
			else
			{
				smoothedFstpeak[jjj - 1] = false;
			}

			previousregression = currentregression;
		} // jjj

		// Now we have to deal with the ends
		// See if there's a local maximum at the end, and if it's in the 
		// range of the other peaks, count it.  Otherwise don't.

		highestpointinendsegment = 0;
		for (jjj = 0; jjj < halfregressioninterval * 2; jjj++)
		{
			smoothedFstpeak[jjj] = false;
			if (smoothed_marker_fsts[jjj] > highestpointinendsegment)
			{
				highestpointinendsegment = smoothed_marker_fsts[jjj];
				indexofhighestpoint = jjj;
			}
		}

		if (indexofhighestpoint != halfregressioninterval * 2 - 1)
		{
			if (highestpointinendsegment > lowestpeak)
				smoothedFstpeak[indexofhighestpoint] = true;
		}

		highestpointinendsegment = 0;
		for (jjj = NumberMarkerLociPerChromosome - halfregressioninterval * 2; jjj < NumberMarkerLociPerChromosome; jjj++)
		{
			smoothedFstpeak[jjj] = false;
			if (smoothed_marker_fsts[jjj] > highestpointinendsegment)
			{
				highestpointinendsegment = smoothed_marker_fsts[jjj];
				indexofhighestpoint = jjj;
			}
		}

		if (indexofhighestpoint != NumberMarkerLociPerChromosome - halfregressioninterval * 2)
		{
			if (highestpointinendsegment > lowestpeak)
			{
				smoothedFstpeak[indexofhighestpoint] = true;
			}
		}
		
		// Remove all peaks that are adjacent to one another (because in this case
		// it's really just one peak).  In each case, keep the one with the higher Fst
		for (jjj = 0; jjj < NumberMarkerLociPerChromosome - 1; jjj++)
		{
			if (smoothedFstpeak[jjj] && smoothedFstpeak[jjj + 1])
			{
				if (marker_fst_values[jjj] > marker_fst_values[jjj + 1])
					smoothedFstpeak[jjj + 1] = false;
				else
					smoothedFstpeak[jjj] = false;
			}
		}

		return 1;
	}

	void compile_peak_data(int peakplateauwidth)
	{
		int jjj, kkk;
		double highestFstnearPeak;
		int highestFstnearPeakindex;
		bool nearbysignificance;
		int platstart, platstop;
		int smallestQTLdist;
		int closestQTLloc;
		int NumberMarkerLociPerChromosome;
		NumberMarkerLociPerChromosome = static_cast<int>(marker_fst_values.size());

		size_t n_qtl_0, n_qtl_1, n_qtl_p;
		n_qtl_0 = qtl_loc_trt0.size();
		n_qtl_1 = qtl_loc_trt1.size();
		n_qtl_p = qtl_loc_pleio.size();

		// Clear all previous data
		
		peakdata.highestFstonplateau.clear();
		peakdata.highestFstIndex.clear();
		peakdata.peaklocation.clear();
		peakdata.peaksmoothedFst.clear();
		peakdata.sig80CI.clear();
		peakdata.sig90CI.clear();
		peakdata.sig95CI.clear();
		peakdata.sig98CI.clear();
		peakdata.sig99CI.clear();
		peakdata.sigFDR.clear();
		peakdata.sig01FSTprime.clear();
		peakdata.sig05FSTprime.clear();
		peakdata.sigFDRFSTprime.clear();
		peakdata.nearestQTL_0.clear();
		peakdata.distancetoQTL_0.clear();
		peakdata.nearestQTL_1.clear();
		peakdata.distancetoQTL_1.clear();
		peakdata.nearestQTL_P.clear();
		peakdata.distancetoQTL_P.clear();

		peakdata.numberofpeaks = 0;
		for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
		{
			if (smoothedFstpeak[jjj])  // Is this marker a peak?
			{
				peakdata.peaksmoothedFst.push_back(smoothed_marker_fsts[jjj]);
				peakdata.peaklocation.push_back(jjj);

				platstart = jjj - peakplateauwidth;
				platstop = jjj + peakplateauwidth + 1;

				if (platstart < 0)
					platstart = 0;
				if (platstop > NumberMarkerLociPerChromosome)
					platstop = NumberMarkerLociPerChromosome;

				highestFstnearPeak = 0;
				highestFstnearPeakindex = 0;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (marker_fst_values[kkk] > highestFstnearPeak)
					{
						highestFstnearPeak = marker_fst_values[kkk];
						highestFstnearPeakindex = kkk;
					}
				}
				peakdata.highestFstonplateau.push_back(highestFstnearPeak);
				peakdata.highestFstIndex.push_back(highestFstnearPeakindex);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (is_significant[kkk])
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sigFDR.push_back(true);
				else
					peakdata.sigFDR.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (smoothed_marker_fsts[kkk] >= SmoothedCritValue99)
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sig99CI.push_back(true);
				else
					peakdata.sig99CI.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (smoothed_marker_fsts[kkk] >= SmoothedCritValue98)
						nearbysignificance = true;
				}
					if (nearbysignificance)
						peakdata.sig98CI.push_back(true);
					else
						peakdata.sig98CI.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (smoothed_marker_fsts[kkk] >= SmoothedCritValue95)
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sig95CI.push_back(true);
				else
					peakdata.sig95CI.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (smoothed_marker_fsts[kkk] >= SmoothedCritValue90)
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sig90CI.push_back(true);
				else
					peakdata.sig90CI.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (smoothed_marker_fsts[kkk] >= SmoothedCritValue80)
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sig80CI.push_back(true);
				else
					peakdata.sig80CI.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (fst_prime_p_value[kkk] <= 0.05)
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sig05FSTprime.push_back(true);
				else
					peakdata.sig05FSTprime.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (fst_prime_p_value[kkk] <= 0.01)
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sig01FSTprime.push_back(true);
				else
					peakdata.sig01FSTprime.push_back(false);

				nearbysignificance = false;
				for (kkk = platstart; kkk < platstop; kkk++)
				{
					if (fst_prime_is_significant_fdrpoint05[kkk])
						nearbysignificance = true;
				}
				if (nearbysignificance)
					peakdata.sigFDRFSTprime.push_back(true);
				else
					peakdata.sigFDRFSTprime.push_back(false);

				closestQTLloc = -1;
				smallestQTLdist = NumberMarkerLociPerChromosome;
				for (kkk = 0; kkk < n_qtl_0; kkk++)
				{
					if (abs(jjj - qtl_loc_trt0[kkk]) < smallestQTLdist)
					{
						smallestQTLdist = abs(jjj - qtl_loc_trt0[kkk]);
						closestQTLloc = qtl_loc_trt0[kkk];
					}
				}
				peakdata.distancetoQTL_0.push_back(smallestQTLdist);
				peakdata.nearestQTL_0.push_back(closestQTLloc);
				
				closestQTLloc = -1;
				smallestQTLdist = NumberMarkerLociPerChromosome;
				for (kkk = 0; kkk < n_qtl_1; kkk++)
				{
					if (abs(jjj - qtl_loc_trt1[kkk]) < smallestQTLdist)
					{
						smallestQTLdist = abs(jjj - qtl_loc_trt1[kkk]);
						closestQTLloc = qtl_loc_trt1[kkk];
					}
				}
				peakdata.distancetoQTL_1.push_back(smallestQTLdist);
				peakdata.nearestQTL_1.push_back(closestQTLloc);

				closestQTLloc = -1;
				smallestQTLdist = NumberMarkerLociPerChromosome;
				for (kkk = 0; kkk < n_qtl_p; kkk++)
				{
					if (abs(jjj - qtl_loc_pleio[kkk]) < smallestQTLdist)
					{
						smallestQTLdist = abs(jjj - qtl_loc_pleio[kkk]);
						closestQTLloc = qtl_loc_pleio[kkk];
					}
				}
				peakdata.distancetoQTL_P.push_back(smallestQTLdist);
				peakdata.nearestQTL_P.push_back(closestQTLloc);
				
				peakdata.numberofpeaks++;

			}

		}
		

	}

};

class Chromosome
{
public:
	int* MarkerLoci;
	double* QTLeffect; // loci affecting trait 0
	double* QTLeffect_1; // loci affecting trait 1
	double* QTLeffect_P0; //effect on trait 0 for pleiotropic loci
	double* QTLeffect_P1; //effect on trait 1 for pleiotropic loci
};

class QTLloc
{
public:
	int* QTLlocation;
	int* QTLlocation_1;
	int* QTLlocation_P;
};

class Gamete
{
public:
	Chromosome* GameteChromosome;
};

class individual
{
public:
	double Genotype[2];
	double Phenotype[2];
	bool Female;
	bool Alive;
	double MatingSuccess;
	Chromosome* PaternalChromosome;
	Chromosome* MaternalChromosome;

	void calculate_genotypic_values(int number_chromosomes, int number_qtl_per, int number_qtl_per_1, int number_qtl_per_p, bool using_epistasis, epistatic_parameter_collection &ep_parms)
	{
		int jj, kk;
		Genotype[0] = 0;
		Genotype[1] = 0;
		for (jj = 0; jj < number_chromosomes; jj++)
		{
			for (kk = 0; kk < number_qtl_per; kk++)
			{
				Genotype[0] = Genotype[0] + PaternalChromosome[jj].QTLeffect[kk] + MaternalChromosome[jj].QTLeffect[kk];
			}
			for (kk = 0; kk < number_qtl_per_1; kk++)
			{
				Genotype[1] = Genotype[1] + PaternalChromosome[jj].QTLeffect_1[kk] + MaternalChromosome[jj].QTLeffect_1[kk];
			}
			for (kk = 0; kk < number_qtl_per_p; kk++)
			{
				Genotype[0] = Genotype[0] + PaternalChromosome[jj].QTLeffect_P0[kk] + MaternalChromosome[jj].QTLeffect_P0[kk];
				Genotype[1] = Genotype[1] + PaternalChromosome[jj].QTLeffect_P1[kk] + MaternalChromosome[jj].QTLeffect_P1[kk];
			}
		}

		// Add epistatic effects
		if (using_epistasis)
		{
			int counter;
			int i, j;
			int n_0, n_1, n_p;

			n_0 = number_chromosomes * number_qtl_per;
			n_1 = number_chromosomes * number_qtl_per_1;
			n_p = number_chromosomes * number_qtl_per_p;

			double* qtl_0 = new double[n_0];
			double* qtl_1 = new double[n_1];
			double* qtl_p0 = new double[n_p];
			double* qtl_p1 = new double[n_p];

			// linearize the qtls across chromosomes

			counter = 0;
			for (i = 0; i < number_chromosomes; i++)
			{
				for (j = 0; j < number_qtl_per; j++)
				{
					qtl_0[counter] = PaternalChromosome[i].QTLeffect[j] + MaternalChromosome[i].QTLeffect[j];
					counter++;
				}
			}

			counter = 0;
			for (i = 0; i < number_chromosomes; i++)
			{
				for (j = 0; j < number_qtl_per_1; j++)
				{
					qtl_1[counter] = PaternalChromosome[i].QTLeffect_1[j] + MaternalChromosome[i].QTLeffect_1[j];
					counter++;
				}
			}

			counter = 0;
			for (i = 0; i < number_chromosomes; i++)
			{
				for (j = 0; j < number_qtl_per_p; j++)
				{
					qtl_p0[counter] = PaternalChromosome[i].QTLeffect_P0[j] + MaternalChromosome[i].QTLeffect_P0[j];
					counter++;
				}
			}

			counter = 0;
			for (i = 0; i < number_chromosomes; i++)
			{
				for (j = 0; j < number_qtl_per_p; j++)
				{
					qtl_p1[counter] = PaternalChromosome[i].QTLeffect_P1[j] + MaternalChromosome[i].QTLeffect_P1[j];
					counter++;
				}
			}

			// Non-pleiotropic 0by0 on 0 epistatic effects
			counter = 0;
			for (i = 0; i < n_0; i++)
			{
				for (j = i + 1; j < n_0; j++)
				{
					Genotype[0] = Genotype[0] + qtl_0[i] * qtl_0[j] * ep_parms.trait0loci_by_trait0loci_effect_on_trait0[counter];
					counter++;
				}
			}

			// Non-pleiotropic 1by1 on 1 epistatic effects
			counter = 0;
			for (i = 0; i < n_1; i++)
			{
				for (j = i + 1; j < n_1; j++)
				{
					Genotype[1] = Genotype[1] + qtl_1[i] * qtl_1[j] * ep_parms.trait1loci_by_trait1loci_effect_on_trait1[counter];
					counter++;
				}
			}

			// Non-pleiotropic 0by1 on 0 epistatic effects
			counter = 0;
			for (i = 0; i < n_0; i++)
			{
				for (j = 0; j < n_1; j++)
				{
					Genotype[0] = Genotype[0] + qtl_0[i] * qtl_1[j] * ep_parms.trait0loci_by_trait1loci_effect_on_trait0[counter];
					counter++;
				}
			}

			// Non-pleiotropic 0by1 on 1 epistatic effects
			counter = 0;
			for (i = 0; i < n_0; i++)
			{
				for (j = 0; j < n_1; j++)
				{
					Genotype[1] = Genotype[1] + qtl_0[i] * qtl_1[j] * ep_parms.trait0loci_by_trait1loci_effect_on_trait1[counter];
					counter++;
				}
			}

			// Non-pleiotropic 1by1 on 0 epistatic effects
			counter = 0;
			for (i = 0; i < n_1; i++)
			{
				for (j = i + 1; j < n_1; j++)
				{
					Genotype[0] = Genotype[0] + qtl_1[i] * qtl_1[j] * ep_parms.trait1loci_by_trait1loci_effect_on_trait0[counter];
					counter++;
				}
			}

			// Non-pleiotropic 0by0 on 1 epistatic effects
			counter = 0;
			for (i = 0; i < n_0; i++)
			{
				for (j = i + 1; j < n_0; j++)
				{
					Genotype[1] = Genotype[1] + qtl_0[i] * qtl_0[j] * ep_parms.trait0loci_by_trait0loci_effect_on_trait1[counter];
					counter++;
				}
			}

			// Pleiotropic 0by0 on 0 epistatic effects
			counter = 0;
			for (i = 0; i < n_p; i++)
			{
				for (j = i + 1; j < n_p; j++)
				{
					Genotype[0] = Genotype[0] + qtl_p0[i] * qtl_p0[j] * ep_parms.pleio_trait0_by_trait0_effect_on_trait0[counter];
					counter++;
				}
			}

			// Pleiotropic 1by1 on 1 epistatic effects
			counter = 0;
			for (i = 0; i < n_p; i++)
			{
				for (j = i + 1; j < n_p; j++)
				{
					Genotype[1] = Genotype[1] + qtl_p1[i] * qtl_p1[j] * ep_parms.pleio_trait1_by_trait1_effect_on_trait1[counter];
					counter++;
				}
			}

			// Pleiotropic 0by1 on 0 epistatic effects
			counter = 0;
			for (i = 0; i < n_p; i++)
			{
				for (j = i + 1; j < n_p; j++)
				{
					Genotype[0] = Genotype[0] + qtl_p0[i] * qtl_p1[j] * ep_parms.pleio_trait0_by_trait1_effect_on_trait0[counter];
					Genotype[0] = Genotype[0] + qtl_p1[i] * qtl_p0[j] * ep_parms.pleio_trait0_by_trait1_effect_on_trait0[counter];
					counter++;
				}
			}

			// Pleiotropic 0by1 on 1 epistatic effects
			counter = 0;
			for (i = 0; i < n_p; i++)
			{
				for (j = i + 1; j < n_p; j++)
				{
					Genotype[1] = Genotype[1] + qtl_p0[i] * qtl_p1[j] * ep_parms.pleio_trait0_by_trait1_effect_on_trait1[counter];
					Genotype[1] = Genotype[1] + qtl_p1[i] * qtl_p0[j] * ep_parms.pleio_trait0_by_trait1_effect_on_trait1[counter];
					counter++;
				}
			}

			// Pleiotropic 0by0 on 1 epistatic effects
			counter = 0;
			for (i = 0; i < n_p; i++)
			{
				for (j = i + 1; j < n_p; j++)
				{
					Genotype[1] = Genotype[1] + qtl_p0[i] * qtl_p0[j] * ep_parms.pleio_trait0_by_trait0_effect_on_trait1[counter];
					counter++;
				}
			}

			// Pleiotropic 1by1 on 0 epistatic effects
			counter = 0;
			for (i = 0; i < n_p; i++)
			{
				for (j = i + 1; j < n_p; j++)
				{
					Genotype[0] = Genotype[0] + qtl_p1[i] * qtl_p1[j] * ep_parms.pleio_trait1_by_trait1_effect_on_trait0[counter];
					counter++;
				}
			}



			delete[] qtl_0;
			delete[] qtl_1;
			delete[] qtl_p0;
			delete[] qtl_p1;

		} // if (using_epistasis)

	}

	void calculate_phenotype(double env_st_dev_0, double env_st_dev_1)
	{
		Phenotype[0] = Genotype[0] + randnorm(0, env_st_dev_0);
		Phenotype[1] = Genotype[1] + randnorm(0, env_st_dev_1);
	}

	void set_sex()
	{
		if (genrand() < 0.5)
			Female = true;
		else
			Female = false;
	}

};

class simulation_data
{
public:
	std::vector<double> vN, vZbar0, vZbar1;
	std::vector<double> vP00, vP11, vP01, vRp, vGbar0, vGbar1;
	std::vector<double> vG00, vG11, vG01, vRg, vLambda1, vLambda2;
	std::vector<double> vEvecX, vEvecY, vAngle, vSize, vEccen;
	std::vector<double> vStrt0, vStrt1, vcLmbd1, vcLmbd2, vcAng;
	std::vector<double> vcSize, vcEcc, vcG00, vcG11, vcG01, vcRg;
	std::vector<double> vASR, vIm, vIf, vMdifM, vMdifF;
	std::vector<double> vV0, vV1, vV01;

	double meanN, meanZbar0, meanZbar1;
	double meanP00, meanP11, meanP01, meanRp, meanGbar0, meanGbar1;
	double meanG00, meanG11, meanG01, meanRg, meanLambda1, meanLambda2;
	double meanEvecX, meanEvecY, meanAngle, meanSize, meanEccen;
	double meanStrt0, meanStrt1, meancLmbd1, meancLmbd2, meancAng;
	double meancSize, meancEcc, meancG00, meancG11, meancG01, meancRg;
	double meanASR, meanIm, meanIf, meanMdifM, meanMdifF;
	double meanV0, meanV1, meanV01;

	void calculate_means()
	{
		size_t i;
		double dNgens = 0;

		// Set all means to zero
		meanN = 0;
		meanZbar0 = 0;
		meanZbar1 = 0;
		meanP00 = 0;
		meanP11 = 0;
		meanP01 = 0;
		meanRp = 0;
		meanGbar0 = 0;
		meanGbar1 = 0;
		meanG00 = 0;
		meanG11 = 0;
		meanG01 = 0;
		meanRg = 0;
		meanLambda1 = 0;
		meanLambda2 = 0;
		meanEvecX = 0;
		meanEvecY = 0;
		meanAngle = 0;
		meanSize = 0;
		meanEccen = 0;
		meanStrt0 = 0;
		meanStrt1 = 0;
		meancLmbd1 = 0;
		meancLmbd2 = 0;
		meancAng = 0;
		meancSize = 0;
		meancEcc = 0;
		meancG00 = 0;
		meancG11 = 0;
		meancG01 = 0;
		meancRg = 0;
		meanASR = 0;
		meanIm = 0;
		meanIf = 0;
		meanMdifM = 0;
		meanMdifF = 0;
		meanV0 = 0; 
		meanV1 = 0;
		meanV01 = 0;

		for (i = 0; i < vN.size(); i++)
		{
			meanN = meanN + vN[i];
			meanZbar0 = meanZbar0 + vZbar0[i];
			meanZbar1 = meanZbar1 + vZbar1[i];
			meanP00 = meanP00 + vP00[i];
			meanP11 = meanP11 + vP11[i];
			meanP01 = meanP01 + vP01[i];
			meanRp = meanRp + vRp[i];
			meanGbar0 = meanGbar0 + vGbar0[i];
			meanGbar1 = meanGbar1 + vGbar1[i];
			meanG00 = meanG00 + vG00[i];
			meanG11 = meanG11 + vG11[i];
			meanG01 = meanG01 + vG01[i];
			meanRg = meanRg + vRg[i];
			meanLambda1 = meanLambda1 + vLambda1[i];
			meanLambda2 = meanLambda2 + vLambda2[i];
			meanEvecX = meanEvecX + vEvecX[i];
			meanEvecY = meanEvecY + vEvecY[i];
			meanAngle = meanAngle + vAngle[i];
			meanSize = meanSize + vSize[i];
			meanEccen = meanEccen + vEccen[i];
			meanStrt0 = meanStrt0 + vStrt0[i];
			meanStrt1 = meanStrt1 + vStrt1[i];
			meancLmbd1 = meancLmbd1 + vcLmbd1[i];
			meancLmbd2 = meancLmbd2 + vcLmbd2[i];
			meancAng = meancAng + vcAng[i];
			meancSize = meancSize + vcSize[i];
			meancEcc = meancEcc + vcEcc[i];
			meancG00 = meancG00 + vcG00[i];
			meancG11 = meancG11 + vcG11[i];
			meancG01 = meancG01 + vcG01[i];
			meancRg = meancRg + vcRg[i];
			meanASR = meanASR + vASR[i];
			meanIm = meanIm + vIm[i];
			meanIf = meanIf + vIf[i];
			meanMdifM = meanMdifM + vMdifM[i];
			meanMdifF = meanMdifF + vMdifF[i];
			meanV0 = meanV0 + vV0[i];
			meanV1 = meanV1 + vV1[i];
			meanV01 = meanV01 + vV01[i];

			dNgens++;
		}

		meanN = meanN / dNgens;
		meanZbar0 = meanZbar0 / dNgens;
		meanZbar1 = meanZbar1 / dNgens;
		meanP00 = meanP00 / dNgens;
		meanP11 = meanP11 / dNgens;
		meanP01 = meanP01 / dNgens;
		meanRp = meanRp / dNgens;
		meanGbar0 = meanGbar0 / dNgens;
		meanGbar1 = meanGbar1 / dNgens;
		meanG00 = meanG00 / dNgens;
		meanG11 = meanG11 / dNgens;
		meanG01 = meanG01 / dNgens;
		meanRg = meanRg / dNgens;
		meanLambda1 = meanLambda1 / dNgens;
		meanLambda2 = meanLambda2 / dNgens;
		meanEvecX = meanEvecX / dNgens;
		meanEvecY = meanEvecY / dNgens;
		meanAngle = meanAngle / dNgens;
		meanSize = meanSize / dNgens;
		meanEccen = meanEccen / dNgens;
		meanStrt0 = meanStrt0 / dNgens;
		meanStrt1 = meanStrt1 / dNgens;
		meancLmbd1 = meancLmbd1 / dNgens;
		meancLmbd2 = meancLmbd2 / dNgens;
		meancAng = meancAng / dNgens;
		meancSize = meancSize / dNgens;
		meancEcc = meancEcc / dNgens;
		meancG00 = meancG00 / dNgens;
		meancG11 = meancG11 / dNgens;
		meancG01 = meancG01 / dNgens;
		meancRg = meancRg / dNgens;
		meanASR = meanASR / dNgens;
		meanIm = meanIm / dNgens;
		meanIf = meanIf / dNgens;
		meanMdifM = meanMdifM / dNgens;
		meanMdifF = meanMdifF / dNgens;
		meanV0 = meanV0 / dNgens;
		meanV1 = meanV1 / dNgens;
		meanV01 = meanV01 / dNgens;
	}

};

class mean_recorder
{
public:
	std::vector<double> mean_list;
};

class marker_fst_data
{
public:
	int chromosome;
	int location;
	double fst;
	double smoothed;
};

class qtl_fst_data
{
public:
	int chromosome;
	int location;
	int type; // 0 = trait0, 1 = trait1, 2 = pleiotropic
	double fst;
};

class simulation_engine
{
private:
	individual *adult;
	individual *progeny;
	individual *temp_progeny;
	int NadultMax, NprogMax;
	int NumberOfGenerations;
	int PopulationSize;
	double MutationalCorrelation, SelectionalCorrelation;
	int Fecundity; // Number of offspring each female can produce
	double MutationRatePerLocus; // The per-locus mutation rate
	int CarryingCapacity; // The maximum adult population size
	double MutationalVariance[2];
	double SelectionStrength[2];
	double EnvironmentalVariance[2];
	double EnvironmentalStDev[2];
	int MaxMatingEncounters;
	int Nprogeny;
	double GaussianPreferenceVariance;
	bool PopulationExtinct;
	int NumberOfInitialGenerations;
	double InitialSelectionStrength[2], InitialSelectionalCorrelation;
	double InitialOptimum[2];
	double Optimum[2];
	bool InitialSelectionSexLimited, ExperimentalSelectionSexLimited;
	double fst_weighting_variance;
	int NumberOfInterveningGens;

	// Variables for Chromosomes and Marker Loci
	int NumberChromosomes;
	int NumberMarkerLociPerChromosome;
	int NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P;
	QTLloc* QTLlocus; // An array of QTL locations that spans the chromosomes
	double MutationRatePerMarker;
	double ExpRecombPerChromosome;
	double ExpectedMarkerMutationsPerChromosome;
	int MaxNumAllelesPerMarker;
	double ExpectedQTLMutationsPerChromosome, ExpectedQTLMutationsPerChromosome_1, ExpectedQTLMutationsPerChromosome_P;
	double StartingQTLAllelicEffectStdDev;
	double sqrtEMMPC, sqrtEQTLMPC, sqrtEQTLMPC_1, sqrtEQTLMPC_P;

	// Some additional new parameters
	int SampleSizeAdults, ActualAdultSampleSize;
	int *adult_sample;
	int all_calc_interval;

	// Epistatic parameters
	bool epistasis_allowed;
	epistatic_parameter_collection ep_par_collection;
	epistatic_variance_list ep_var_list;

	// Variables corresponding to population-level summary statistics:
	double phenotypic_mean[2], genotypic_mean[2];
	double phenotypic_variance[2], genotypic_variance[2];
	double phenotypic_covariance, genotypic_covariance;
	double phenotypic_correlation, genotypic_correlation;
	double EigenValue[2];
	double EigenVector1[2];
	double EigenVector2[2];
	double LeadAngle, Sigma, Epsilon;
	double sel_diff_trt_0, sel_diff_trt_1;
	double est_va00, est_va11, est_va01;

	double PrevEval[2], PrevAngle, PrevSigma;
	double PrevEpsilon, PrevG00, PrevG11, PrevG01, PrevRg;
	double cPrevEval[2], cPrevAngle, cPrevSigma;
	double cPrevEpsilon, cPrevG00, cPrevG11, cPrevG01, cPrevRg;

	double SexRatio, Im, If, MdiffMales, MdiffFemales;

	simulation_data sim_data;
	mean_recorder m_rec;
	std::vector<chromosomedata> chr_data;

public:
	individual migrant;

	simulation_engine() // Empty. Set default parameter values in the parse command line function (above)
	{

	}

	void display_parameters()
	{
		std::cout << "Parameter_Values:\n";
		// Demographic Parameters
		std::cout << "Demographic_Parameters:\n";
		std::cout << "No_Generations:   \t" << NumberOfGenerations << "\n";
		std::cout << "Initial_Pop_Size: \t" << PopulationSize << "\n";
		std::cout << "Carrying_Capacity:\t" << CarryingCapacity << "\n";
		std::cout << "Female_Fecundity: \t" << Fecundity << "\n";

		// Mating Parameters
		std::cout << "Mating_Parameters:\n";
		std::cout << "Max_Mating_Enc:   \t" << MaxMatingEncounters << "\n";
		std::cout << "Gaussian_Pref_Var:\t" << GaussianPreferenceVariance << "\n";

		// Quantitative Genetic Parameters
		std::cout << "Quantitative_Genetic_Parameters:\n";
		std::cout << "Env_Variance_Trt0:\t" << EnvironmentalVariance[0] << "\n";
		std::cout << "Env_Variance_Trt1:\t" << EnvironmentalVariance[1] << "\n";

		// Mutational Parameters
		std::cout << "Mutational_Parameters:\n";
		std::cout << "Mut_Var_Trait0:   \t" << MutationalVariance[0] << "\n";
		std::cout << "Mut_Var_Trait1:   \t" << MutationalVariance[1] << "\n";
		std::cout << "Mut_Correlation:  \t" << MutationalCorrelation << "\n";
		std::cout << "Mutation_Rate:    \t" << MutationRatePerLocus << "\n";

		// Selection Parameters
		std::cout << "Selection_Parameters:\n";
		std::cout << "Omega_Trait0:     \t" << SelectionStrength[0] << "\n";
		std::cout << "Omega_Trait1:     \t" << SelectionStrength[1] << "\n";
		std::cout << "Selection_Corr:   \t" << SelectionalCorrelation << "\n";
		std::cout << "Optimum_Trait0    \t" << Optimum[0] << "\n";
		std::cout << "Optimum_Trait1    \t" << Optimum[1] << "\n";
		if (ExperimentalSelectionSexLimited)
			std::cout << "Sex_Limited_Sel:  \ttrue\n";
		else
			std::cout << "Sex_Limited_Sel:  \tfalse\n";

		// Initial Generations
		std::cout << "Initial_Generations_Parameters:\n";
		std::cout << "No_Initial_Gens:  \t" << NumberOfInitialGenerations << "\n";
		std::cout << "Init_Omega_Trait0:\t" << InitialSelectionStrength[0] << "\n";
		std::cout << "Init_Omega_Trait1:\t" << InitialSelectionStrength[1] << "\n";
		std::cout << "Init_Sel_Corr:    \t" << InitialSelectionalCorrelation << "\n";
		std::cout << "Init_Opt_Trt0     \t" << InitialOptimum[0] << "\n";
		std::cout << "Init_Opt_Trt1     \t" << InitialOptimum[1] << "\n";
		if (InitialSelectionSexLimited)
			std::cout << "Init_Sex_Lim_Sel: \ttrue\n";
		else
			std::cout << "Init_Sex_Lim_Sel: \tfalse\n";
	}

	void initialize_population()
	{
		int i, j, k;

		// Initialize previous generation variables to zero
		PrevEval[0] = 0;
		PrevEval[1] = 0;
		PrevAngle = 0;
		PrevSigma = 0;
		PrevEpsilon = 0;
		PrevG00 = 0;
		PrevG11 = 0;
		PrevG01 = 0;
		PrevRg = 0;

		PopulationExtinct = false;

		// Set the maximum number of adults and progeny in the population
		// Allocate memory for the pointers for the adults and progeny
		NadultMax = CarryingCapacity;
		NprogMax = CarryingCapacity * Fecundity;
		adult = new individual[NadultMax];
		progeny = new individual[NprogMax];
		temp_progeny = new individual[10];
		adult_sample = new int[SampleSizeAdults];
		QTLlocus = new QTLloc[NumberChromosomes];
		for (i = 0; i < NumberChromosomes; i++)
		{
			QTLlocus[i].QTLlocation = new int[NumberQTLsPerChromosome];
			QTLlocus[i].QTLlocation_1 = new int[NumberQTLsPerChromosome_1];
			QTLlocus[i].QTLlocation_P = new int[NumberQTLsPerChromosome_P];
		}

		// Allocate memory for the adults
		for (i = 0; i < NadultMax; i++)
		{
			adult[i].MaternalChromosome = new Chromosome[NumberChromosomes];
			adult[i].PaternalChromosome = new Chromosome[NumberChromosomes];
			for (j = 0; j < NumberChromosomes; j++)
			{
				adult[i].MaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				adult[i].PaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				adult[i].MaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
				adult[i].PaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
				adult[i].MaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
				adult[i].PaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
				adult[i].MaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
				adult[i].PaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
				adult[i].MaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
				adult[i].PaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
			}
		}

		// Allocate memory for progeny
		for (i = 0; i < NprogMax; i++)
		{
			progeny[i].MaternalChromosome = new Chromosome[NumberChromosomes];
			progeny[i].PaternalChromosome = new Chromosome[NumberChromosomes];
			for (j = 0; j < NumberChromosomes; j++)
			{
				progeny[i].MaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				progeny[i].PaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				progeny[i].MaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
				progeny[i].PaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
				progeny[i].MaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
				progeny[i].PaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
				progeny[i].MaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
				progeny[i].PaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
				progeny[i].MaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
				progeny[i].PaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
			}
		}

		// Allocate memory for temp_progeny
		for (i = 0; i < 10; i++)
		{
			temp_progeny[i].MaternalChromosome = new Chromosome[NumberChromosomes];
			temp_progeny[i].PaternalChromosome = new Chromosome[NumberChromosomes];
			for (j = 0; j < NumberChromosomes; j++)
			{
				temp_progeny[i].MaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				temp_progeny[i].PaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				temp_progeny[i].MaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
				temp_progeny[i].PaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
				temp_progeny[i].MaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
				temp_progeny[i].PaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
				temp_progeny[i].MaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
				temp_progeny[i].PaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
				temp_progeny[i].MaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
				temp_progeny[i].PaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
			}
		}

		// Allocate memory for migrant
		migrant.MaternalChromosome = new Chromosome[NumberChromosomes];
		migrant.PaternalChromosome = new Chromosome[NumberChromosomes];
		for (j = 0; j < NumberChromosomes; j++)
		{
			migrant.MaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
			migrant.PaternalChromosome[j].MarkerLoci = new int[NumberMarkerLociPerChromosome];
			migrant.MaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
			migrant.PaternalChromosome[j].QTLeffect = new double[NumberQTLsPerChromosome];
			migrant.MaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
			migrant.PaternalChromosome[j].QTLeffect_1 = new double[NumberQTLsPerChromosome_1];
			migrant.MaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
			migrant.PaternalChromosome[j].QTLeffect_P0 = new double[NumberQTLsPerChromosome_P];
			migrant.MaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
			migrant.PaternalChromosome[j].QTLeffect_P1 = new double[NumberQTLsPerChromosome_P];
		}

		// Initialize Adults
		double* tempQTLallele = new double[MaxNumAllelesPerMarker];
		double* tempQTLallele_1 = new double[MaxNumAllelesPerMarker];
		double* tempQTLallele_P0 = new double[MaxNumAllelesPerMarker];
		double* tempQTLallele_P1 = new double[MaxNumAllelesPerMarker];
		
		for (j = 0; j < NumberChromosomes; j++)
		{
			for (i = 0; i < MaxNumAllelesPerMarker; i++)
			{
				tempQTLallele[i] = randnorm(0, StartingQTLAllelicEffectStdDev);
				tempQTLallele_1[i] = randnorm(0, StartingQTLAllelicEffectStdDev);
				randbivnorm(StartingQTLAllelicEffectStdDev, StartingQTLAllelicEffectStdDev, MutationalCorrelation, tempQTLallele_P0[i], tempQTLallele_P1[i]);
			}

			for (i = 0; i < NadultMax; i++)
			{
				for (k = 0; k < NumberMarkerLociPerChromosome; k++)
				{
					adult[i].MaternalChromosome[j].MarkerLoci[k] = i % MaxNumAllelesPerMarker;
					adult[i].PaternalChromosome[j].MarkerLoci[k] = i % MaxNumAllelesPerMarker;
				}

				for (k = 0; k < NumberQTLsPerChromosome; k++)
				{
					adult[i].MaternalChromosome[j].QTLeffect[k] = tempQTLallele[i%MaxNumAllelesPerMarker];
					adult[i].PaternalChromosome[j].QTLeffect[k] = tempQTLallele[i%MaxNumAllelesPerMarker];
				}
				for (k = 0; k < NumberQTLsPerChromosome_1; k++)
				{
					adult[i].MaternalChromosome[j].QTLeffect_1[k] = tempQTLallele_1[i%MaxNumAllelesPerMarker];
					adult[i].PaternalChromosome[j].QTLeffect_1[k] = tempQTLallele_1[i%MaxNumAllelesPerMarker];
				}
				for (k = 0; k < NumberQTLsPerChromosome_P; k++)
				{
					adult[i].MaternalChromosome[j].QTLeffect_P0[k] = tempQTLallele_P0[i%MaxNumAllelesPerMarker];
					adult[i].PaternalChromosome[j].QTLeffect_P0[k] = tempQTLallele_P0[i%MaxNumAllelesPerMarker];
					adult[i].MaternalChromosome[j].QTLeffect_P1[k] = tempQTLallele_P1[i%MaxNumAllelesPerMarker];
					adult[i].PaternalChromosome[j].QTLeffect_P1[k] = tempQTLallele_P1[i%MaxNumAllelesPerMarker];
				}
				
			}
		}

		for (i = 0; i < NadultMax; i++)
		{
			adult[i].calculate_genotypic_values(NumberChromosomes, NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P, epistasis_allowed, ep_par_collection);
			adult[i].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
			adult[i].set_sex();
		} // end of i

		// Initialize progeny with zeros for everything
		for (i = 0; i < NprogMax; i++)
		{
			for (j = 0; j < NumberChromosomes; j++)
			{
				for (k = 0; k < NumberMarkerLociPerChromosome; k++)
				{
					progeny[i].MaternalChromosome[j].MarkerLoci[k] = 0;
					progeny[i].PaternalChromosome[j].MarkerLoci[k] = 0;
				}

				for (k = 0; k < NumberQTLsPerChromosome; k++)
				{
					progeny[i].MaternalChromosome[j].QTLeffect[k] = 0;
					progeny[i].PaternalChromosome[j].QTLeffect[k] = 0;
				}

				for (k = 0; k < NumberQTLsPerChromosome_1; k++)
				{
					progeny[i].MaternalChromosome[j].QTLeffect_1[k] = 0;
					progeny[i].PaternalChromosome[j].QTLeffect_1[k] = 0;
				}

				for (k = 0; k < NumberQTLsPerChromosome_P; k++)
				{
					progeny[i].MaternalChromosome[j].QTLeffect_P0[k] = 0;
					progeny[i].PaternalChromosome[j].QTLeffect_P0[k] = 0;
					progeny[i].MaternalChromosome[j].QTLeffect_P1[k] = 0;
					progeny[i].PaternalChromosome[j].QTLeffect_P1[k] = 0;
				}

			}

			progeny[i].calculate_genotypic_values(NumberChromosomes, NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P, epistasis_allowed, ep_par_collection);
			progeny[i].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
			progeny[i].set_sex();
		} // end of i

		// Set the QTL locations -- for the multi-pop model, these need to be set to be consistent across pops at the metapop level
		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberQTLsPerChromosome; j++)
				QTLlocus[i].QTLlocation[j] = randnum(NumberMarkerLociPerChromosome);
			for (j = 0; j < NumberQTLsPerChromosome_1; j++)
				QTLlocus[i].QTLlocation_1[j] = randnum(NumberMarkerLociPerChromosome);
			for (j = 0; j < NumberQTLsPerChromosome_P; j++)
				QTLlocus[i].QTLlocation_P[j] = randnum(NumberMarkerLociPerChromosome);
		}

		// prepare chromosome data vector for use
		
		chr_data.resize(NumberChromosomes);
		for (i = 0; i < NumberChromosomes; i++)
		{
			chr_data[i].marker_afs.resize(NumberMarkerLociPerChromosome);
			chr_data[i].qtl_afs_trt0.resize(NumberQTLsPerChromosome);
			chr_data[i].qtl_afs_trt1.resize(NumberQTLsPerChromosome_1);
			chr_data[i].qtl_afs_pleio.resize(NumberQTLsPerChromosome_P);
			chr_data[i].qtl_loc_trt0.resize(NumberQTLsPerChromosome);
			chr_data[i].qtl_loc_trt1.resize(NumberQTLsPerChromosome_1);
			chr_data[i].qtl_loc_pleio.resize(NumberQTLsPerChromosome_P);
		}
	
		delete[] tempQTLallele;
		delete[] tempQTLallele_1;
		delete[] tempQTLallele_P0;
		delete[] tempQTLallele_P1;
		
	}

	void deinitialize_population()
	{
		int i, j;

		for (i = 0; i < NadultMax; i++)
		{
			for (j = 0; j < NumberChromosomes; j++)
			{
				delete[] adult[i].MaternalChromosome[j].MarkerLoci;
				delete[] adult[i].PaternalChromosome[j].MarkerLoci;
				delete[] adult[i].MaternalChromosome[j].QTLeffect;
				delete[] adult[i].PaternalChromosome[j].QTLeffect;
				delete[] adult[i].MaternalChromosome[j].QTLeffect_1;
				delete[] adult[i].PaternalChromosome[j].QTLeffect_1;
				delete[] adult[i].MaternalChromosome[j].QTLeffect_P0;
				delete[] adult[i].PaternalChromosome[j].QTLeffect_P0;
				delete[] adult[i].MaternalChromosome[j].QTLeffect_P1;
				delete[] adult[i].PaternalChromosome[j].QTLeffect_P1;
			}
			delete[] adult[i].MaternalChromosome;
			delete[] adult[i].PaternalChromosome;
		}

		for (i = 0; i < NprogMax; i++)
		{
			for (j = 0; j < NumberChromosomes; j++)
			{
				delete[] progeny[i].MaternalChromosome[j].MarkerLoci;
				delete[] progeny[i].PaternalChromosome[j].MarkerLoci;
				delete[] progeny[i].MaternalChromosome[j].QTLeffect;
				delete[] progeny[i].PaternalChromosome[j].QTLeffect;
				delete[] progeny[i].MaternalChromosome[j].QTLeffect_1;
				delete[] progeny[i].PaternalChromosome[j].QTLeffect_1;
				delete[] progeny[i].MaternalChromosome[j].QTLeffect_P0;
				delete[] progeny[i].PaternalChromosome[j].QTLeffect_P0;
				delete[] progeny[i].MaternalChromosome[j].QTLeffect_P1;
				delete[] progeny[i].PaternalChromosome[j].QTLeffect_P1;
			}
			delete[] progeny[i].MaternalChromosome;
			delete[] progeny[i].PaternalChromosome;
		}

		for (i = 0; i < 10; i++)
		{
			for (j = 0; j < NumberChromosomes; j++)
			{
				delete[] temp_progeny[i].MaternalChromosome[j].MarkerLoci;
				delete[] temp_progeny[i].PaternalChromosome[j].MarkerLoci;
				delete[] temp_progeny[i].MaternalChromosome[j].QTLeffect;
				delete[] temp_progeny[i].PaternalChromosome[j].QTLeffect;
				delete[] temp_progeny[i].MaternalChromosome[j].QTLeffect_1;
				delete[] temp_progeny[i].PaternalChromosome[j].QTLeffect_1;
				delete[] temp_progeny[i].MaternalChromosome[j].QTLeffect_P0;
				delete[] temp_progeny[i].PaternalChromosome[j].QTLeffect_P0;
				delete[] temp_progeny[i].MaternalChromosome[j].QTLeffect_P1;
				delete[] temp_progeny[i].PaternalChromosome[j].QTLeffect_P1;
			}
			delete[] temp_progeny[i].MaternalChromosome;
			delete[] temp_progeny[i].PaternalChromosome;
		}

		for (j = 0; j < NumberChromosomes; j++)
		{
			delete[] migrant.MaternalChromosome[j].MarkerLoci;
			delete[] migrant.PaternalChromosome[j].MarkerLoci;
			delete[] migrant.MaternalChromosome[j].QTLeffect;
			delete[] migrant.PaternalChromosome[j].QTLeffect;
			delete[] migrant.MaternalChromosome[j].QTLeffect_1;
			delete[] migrant.PaternalChromosome[j].QTLeffect_1;
			delete[] migrant.MaternalChromosome[j].QTLeffect_P0;
			delete[] migrant.PaternalChromosome[j].QTLeffect_P0;
			delete[] migrant.MaternalChromosome[j].QTLeffect_P1;
			delete[] migrant.PaternalChromosome[j].QTLeffect_P1;
		}
		delete[] migrant.MaternalChromosome;
		delete[] migrant.PaternalChromosome;

		for (i = 0; i < NumberChromosomes; i++)
			delete[] QTLlocus[i].QTLlocation;

		delete[] adult;
		delete[] progeny;
		delete[] temp_progeny;
		delete[] adult_sample;
		delete[] QTLlocus;
	}

	void ProduceRecombinedChromosome(Chromosome &RecombinedChromosome, individual &Parent, int WhichChromosome, double ExpectedRecombinationEvents) // Chromosome Result, Adult Parent, Chromosome Number, Max of 20 events
	{
		int PRi, PRj;
		int NumberRecombinationEvents = 0;
		int SegmentStart[22], SegmentEnd[22];
		int BreakPoint[20];

		if (ExpectedRecombinationEvents < 6)
			NumberRecombinationEvents = poissonrand(ExpectedRecombinationEvents);
		if (ExpectedRecombinationEvents >= 6)
			NumberRecombinationEvents = positiveroundnorm(ExpectedRecombinationEvents, sqrt(ExpectedRecombinationEvents));

		for (PRi = 0; PRi < 20; PRi++)
			BreakPoint[PRi] = NumberMarkerLociPerChromosome + 1; // Make sure the unused breakpoints in the array are larger than the used breakpoints
		bool SegmentMaternal[22];
		int NumberSegments;
		bool StartMaternal;
		if (NumberRecombinationEvents > 20)
			NumberRecombinationEvents = 20;

		if (NumberRecombinationEvents > 0)
		{
			for (PRi = 0; PRi < NumberRecombinationEvents; PRi++)
				BreakPoint[PRi] = randnum(NumberMarkerLociPerChromosome);

			// Have to sort the breakpoints in ascending order
			std::sort(BreakPoint, BreakPoint+20);

			// Is the first segment maternal or paternal?
			if (genrand() < 0.5)
				StartMaternal = true;
			else
				StartMaternal = false;

			NumberSegments = 1;
			SegmentStart[0] = 0;
			SegmentMaternal[0] = StartMaternal;
			for (PRi = 0; PRi < NumberRecombinationEvents; PRi++)
			{
				SegmentEnd[PRi] = BreakPoint[PRi];
				SegmentStart[PRi + 1] = BreakPoint[PRi];
				if (SegmentMaternal[PRi])
					SegmentMaternal[PRi + 1] = false;
				else
					SegmentMaternal[PRi + 1] = true;
				NumberSegments++;
			} // end of PRi
			SegmentEnd[PRi] = NumberMarkerLociPerChromosome;

			// Finally we can pass the allelic information to the recombined chromosome
			for (PRi = 0; PRi < NumberSegments; PRi++)
			{
				if (SegmentMaternal[PRi])
				{
					for (PRj = SegmentStart[PRi]; PRj < SegmentEnd[PRi]; PRj++)
						RecombinedChromosome.MarkerLoci[PRj] = Parent.MaternalChromosome[WhichChromosome].MarkerLoci[PRj];
				}
				else
				{
					for (PRj = SegmentStart[PRi]; PRj < SegmentEnd[PRi]; PRj++)
						RecombinedChromosome.MarkerLoci[PRj] = Parent.PaternalChromosome[WhichChromosome].MarkerLoci[PRj];
				}
			} // end of PRi

			// Next do it for the QTLs for trait 0
			for (PRi = 0; PRi < NumberQTLsPerChromosome; PRi++)
			{
				for (PRj = 0; PRj < NumberSegments; PRj++)
				{
					if (QTLlocus[WhichChromosome].QTLlocation[PRi] >= SegmentStart[PRj] && QTLlocus[WhichChromosome].QTLlocation[PRi] < SegmentEnd[PRj])
					{
						if (SegmentMaternal[PRj])
							RecombinedChromosome.QTLeffect[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect[PRi];
						else
							RecombinedChromosome.QTLeffect[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect[PRi];
					}
				}
			} // end of PRi

			// Next do it for the QTLs for trait 1
			for (PRi = 0; PRi < NumberQTLsPerChromosome_1; PRi++)
			{
				for (PRj = 0; PRj < NumberSegments; PRj++)
				{
					if (QTLlocus[WhichChromosome].QTLlocation_1[PRi] >= SegmentStart[PRj] && QTLlocus[WhichChromosome].QTLlocation_1[PRi] < SegmentEnd[PRj])
					{
						if (SegmentMaternal[PRj])
							RecombinedChromosome.QTLeffect_1[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect_1[PRi];
						else
							RecombinedChromosome.QTLeffect_1[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect_1[PRi];
					}
				}
			} // end of PRi

			// Next do it for the QTLs for pleiotropic loci
			for (PRi = 0; PRi < NumberQTLsPerChromosome_P; PRi++)
			{
				for (PRj = 0; PRj < NumberSegments; PRj++)
				{
					if (QTLlocus[WhichChromosome].QTLlocation_P[PRi] >= SegmentStart[PRj] && QTLlocus[WhichChromosome].QTLlocation_P[PRi] < SegmentEnd[PRj])
					{
						if (SegmentMaternal[PRj])
						{
							RecombinedChromosome.QTLeffect_P0[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect_P0[PRi];
							RecombinedChromosome.QTLeffect_P1[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect_P1[PRi];
						}
						else
						{
							RecombinedChromosome.QTLeffect_P0[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect_P0[PRi];
							RecombinedChromosome.QTLeffect_P1[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect_P1[PRi];
						}
					}
				}
			} // end of PRi

		}  // end of if(NumberRecombinationEvents > 0)
		else
		{
			// No recombination
			if (genrand() < 0.5)
			{
				for (PRi = 0; PRi < NumberMarkerLociPerChromosome; PRi++)
					RecombinedChromosome.MarkerLoci[PRi] = Parent.MaternalChromosome[WhichChromosome].MarkerLoci[PRi];
				for (PRi = 0; PRi < NumberQTLsPerChromosome; PRi++)
					RecombinedChromosome.QTLeffect[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect[PRi];
				for (PRi = 0; PRi < NumberQTLsPerChromosome_1; PRi++)
					RecombinedChromosome.QTLeffect_1[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect_1[PRi];
				for (PRi = 0; PRi < NumberQTLsPerChromosome_P; PRi++)
				{
					RecombinedChromosome.QTLeffect_P0[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect_P0[PRi];
					RecombinedChromosome.QTLeffect_P1[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect_P1[PRi];
				}
			}
			else
			{
				for (PRi = 0; PRi < NumberMarkerLociPerChromosome; PRi++)
					RecombinedChromosome.MarkerLoci[PRi] = Parent.PaternalChromosome[WhichChromosome].MarkerLoci[PRi];
				for (PRi = 0; PRi < NumberQTLsPerChromosome; PRi++)
					RecombinedChromosome.QTLeffect[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect[PRi];
				for (PRi = 0; PRi < NumberQTLsPerChromosome_1; PRi++)
					RecombinedChromosome.QTLeffect_1[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect_1[PRi];
				for (PRi = 0; PRi < NumberQTLsPerChromosome_P; PRi++)
				{
					RecombinedChromosome.QTLeffect_P0[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect_P0[PRi];
					RecombinedChromosome.QTLeffect_P1[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect_P1[PRi];
				}
			}
		} // else
	}  // End of ProduceRecombinedChromosome

	void deep_copy_individuals(individual &to, individual &from)
	{
		int ii, jj;
		to.Alive = from.Alive;
		to.Female = from.Female;
		to.Genotype[0] = from.Genotype[0];
		to.Genotype[1] = from.Genotype[1];
		to.Phenotype[0] = from.Phenotype[0];
		to.Phenotype[1] = from.Genotype[1];
		to.MatingSuccess = from.MatingSuccess;
		for (ii = 0; ii < NumberChromosomes; ii++)
		{
			for (jj = 0; jj < NumberMarkerLociPerChromosome; jj++)
			{
				to.MaternalChromosome[ii].MarkerLoci[jj] = from.MaternalChromosome[ii].MarkerLoci[jj];
				to.PaternalChromosome[ii].MarkerLoci[jj] = from.PaternalChromosome[ii].MarkerLoci[jj];
			}
			for (jj = 0; jj < NumberQTLsPerChromosome; jj++)
			{
				to.MaternalChromosome[ii].QTLeffect[jj] = from.MaternalChromosome[ii].QTLeffect[jj];
				to.PaternalChromosome[ii].QTLeffect[jj] = from.PaternalChromosome[ii].QTLeffect[jj];
			}
			for (jj = 0; jj < NumberQTLsPerChromosome_1; jj++)
			{
				to.MaternalChromosome[ii].QTLeffect_1[jj] = from.MaternalChromosome[ii].QTLeffect_1[jj];
				to.PaternalChromosome[ii].QTLeffect_1[jj] = from.PaternalChromosome[ii].QTLeffect_1[jj];
			}
			for (jj = 0; jj < NumberQTLsPerChromosome_P; jj++)
			{
				to.MaternalChromosome[ii].QTLeffect_P0[jj] = from.MaternalChromosome[ii].QTLeffect_P0[jj];
				to.PaternalChromosome[ii].QTLeffect_P0[jj] = from.PaternalChromosome[ii].QTLeffect_P0[jj];
				to.MaternalChromosome[ii].QTLeffect_P1[jj] = from.MaternalChromosome[ii].QTLeffect_P1[jj];
				to.PaternalChromosome[ii].QTLeffect_P1[jj] = from.PaternalChromosome[ii].QTLeffect_P1[jj];
			}
		}
	}

	void output_adults()
	{
		std::cout << "\n\nID\tSex\tGeno\tPheno\t";
		int i, j, k;

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberMarkerLociPerChromosome; j++)
				std::cout << "C" << i << "M" << j << "\t";
			for (j = 0; j < NumberChromosomes; j++)
				std::cout << "C" << i << "Q" << j << "\t";
		}

		for (i = 0; i < PopulationSize; i++)
		{
				std::cout << "\n" << i << "\t";
				if (adult[i].Female)
					std::cout << "female\t";
				else
					std::cout << "male\t";
				std::cout << std::setprecision(3) << std::fixed << adult[i].Genotype[0] << "\t";
				std::cout << std::setprecision(3) << std::fixed << adult[i].Phenotype[0] << "\t";
				
				for (j = 0; j < NumberChromosomes; j++)
				{
					for (k = 0; k < NumberMarkerLociPerChromosome; k++)
					{
						std::cout << adult[i].MaternalChromosome[j].MarkerLoci[k] << "x" << adult[i].PaternalChromosome[j].MarkerLoci[k] << "\t";
					}
					for (k = 0; k < NumberQTLsPerChromosome; k++)
					{
						std::cout << adult[i].MaternalChromosome[j].QTLeffect[k] << "x" << adult[i].PaternalChromosome[j].QTLeffect[k] << "\t";
					}
				}
		} // end of i loop

	}

	void polygynous_mating()
	{
		// This function implements strict polygyny.
		// Under this mating system, each female mates once, but
		// each male can mate an unlimited number of times.
		// Females choose males at random.

		int i, m, n;
		int iPC;
		bool mate_found;
		int mateID, counter, rnum;

		// Check to make sure at least one male is present in the population
		bool males_present = false;
		for (i = 0; i < PopulationSize; i++)
		{
			if (!adult[i].Female)
				males_present = true;
		}

		iPC = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female && males_present)
			{

				// Find a mate for this female
				// Mating is random, so any male will do

				mate_found = false;
				counter = 0;
				while (counter < MaxMatingEncounters && !mate_found)
				{
					rnum = randnum(PopulationSize);
					if (!adult[rnum].Female)
					{
						mateID = rnum;
						mate_found = true;
					}
					counter++;
				} // end of while

				// If a mate is found, produce progeny
				if (mate_found)
				{
					for (m = 0; m < Fecundity; m++)
					{
						if (iPC >= NprogMax)
							iPC = NprogMax - 1;

						for (n = 0; n < NumberChromosomes; n++)
						{
							ProduceRecombinedChromosome(progeny[iPC].MaternalChromosome[n], adult[i], n, ExpRecombPerChromosome);
							ProduceRecombinedChromosome(progeny[iPC].PaternalChromosome[n], adult[mateID], n, ExpRecombPerChromosome);
						}

						progeny[iPC].set_sex();
						progeny[iPC].calculate_genotypic_values(NumberChromosomes, NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P, epistasis_allowed, ep_par_collection);
						progeny[iPC].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
						progeny[iPC].Alive = true;
						iPC++;

					} // end of m loop
					adult[i].MatingSuccess++;
					adult[mateID].MatingSuccess++;
				} // end of if (mate_found)
			} // end of if (adult[i].Female && males_present)
		} // end of i loop

		Nprogeny = iPC;

	}

	void output_progeny()
	{
		std::cout << "\n\nID\tSex\tGeno\tPheno\t";
		int i, j, k;

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberMarkerLociPerChromosome; j++)
				std::cout << "C" << i << "M" << j << "\t";
			for (j = 0; j < NumberChromosomes; j++)
				std::cout << "C" << i << "Q" << j << "\t";
		}

		for (i = 0; i < Nprogeny; i++)
		{
			std::cout << "\n" << i << "\t";
			if (progeny[i].Female)
				std::cout << "female\t";
			else
				std::cout << "male\t";
			std::cout << std::setprecision(3) << std::fixed << progeny[i].Genotype[0] << "\t";
			std::cout << std::setprecision(3) << std::fixed << progeny[i].Phenotype[0] << "\t";

			for (j = 0; j < NumberChromosomes; j++)
			{
				for (k = 0; k < NumberMarkerLociPerChromosome; k++)
				{
					std::cout << progeny[i].MaternalChromosome[j].MarkerLoci[k] << "x" << progeny[i].PaternalChromosome[j].MarkerLoci[k] << "\t";
				}
				for (k = 0; k < NumberQTLsPerChromosome; k++)
				{
					std::cout << progeny[i].MaternalChromosome[j].QTLeffect[k] << "x" << progeny[i].PaternalChromosome[j].QTLeffect[k] << "\t";
				}
			}

			if (progeny[i].Alive)
				std::cout << "Alive";
			else
				std::cout << "Dead";

		} // end of i loop

	}

	void mutation()
	{
		int i, j, k, number_mutations;
		int whichlocus;
		double drnum0, drnum1;
		double mutational_std_dev_0, mutational_std_dev_1;
		mutational_std_dev_0 = sqrt(MutationalVariance[0]);
		mutational_std_dev_1 = sqrt(MutationalVariance[1]);

		for (i = 0; i < Nprogeny; i++)
		{
			for (j = 0; j < NumberChromosomes; j++)
			{
				// Mutations at the marker loci
				// Maternal Chromosome
				if (ExpectedMarkerMutationsPerChromosome < 6)
					number_mutations = poissonrand(ExpectedMarkerMutationsPerChromosome);
				else
					number_mutations = positiveroundnorm(ExpectedMarkerMutationsPerChromosome,sqrtEMMPC);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberMarkerLociPerChromosome);
					progeny[i].MaternalChromosome[j].MarkerLoci[whichlocus] = randnum(MaxNumAllelesPerMarker);
				}
				// Paternal Chromosome
				if (ExpectedMarkerMutationsPerChromosome < 6)
					number_mutations = poissonrand(ExpectedMarkerMutationsPerChromosome);
				else
					number_mutations = positiveroundnorm(ExpectedMarkerMutationsPerChromosome, sqrtEMMPC);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberMarkerLociPerChromosome);
					progeny[i].PaternalChromosome[j].MarkerLoci[whichlocus] = randnum(MaxNumAllelesPerMarker);
				}

				// Mutations at the QTL for trait 0
				// Maternal Chromosome
				if (ExpectedQTLMutationsPerChromosome < 6)
					number_mutations = poissonrand(ExpectedQTLMutationsPerChromosome);
				else
					number_mutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome, sqrtEQTLMPC);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome);
					progeny[i].MaternalChromosome[j].QTLeffect[whichlocus] = progeny[i].MaternalChromosome[j].QTLeffect[whichlocus] + randnorm(0, mutational_std_dev_0);
				}

				// Paternal Chromosome
				if (ExpectedQTLMutationsPerChromosome < 6)
					number_mutations = poissonrand(ExpectedQTLMutationsPerChromosome);
				else
					number_mutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome, sqrtEQTLMPC);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome);
					progeny[i].PaternalChromosome[j].QTLeffect[whichlocus] = progeny[i].PaternalChromosome[j].QTLeffect[whichlocus] + randnorm(0, mutational_std_dev_0);
				}

				// Mutations at the QTL for trait 1
				// Maternal Chromosome
				if (ExpectedQTLMutationsPerChromosome_1 < 6)
					number_mutations = poissonrand(ExpectedQTLMutationsPerChromosome_1);
				else
					number_mutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome_1, sqrtEQTLMPC_1);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome_1);
					progeny[i].MaternalChromosome[j].QTLeffect_1[whichlocus] = progeny[i].MaternalChromosome[j].QTLeffect_1[whichlocus] + randnorm(0, mutational_std_dev_1);
				}

				// Paternal Chromosome
				if (ExpectedQTLMutationsPerChromosome_1 < 6)
					number_mutations = poissonrand(ExpectedQTLMutationsPerChromosome_1);
				else
					number_mutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome_1, sqrtEQTLMPC_1);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome_1);
					progeny[i].PaternalChromosome[j].QTLeffect_1[whichlocus] = progeny[i].PaternalChromosome[j].QTLeffect_1[whichlocus] + randnorm(0, mutational_std_dev_1);
				}

				// Mutations at the QTL for pleiotropic loci
				// Maternal Chromosome
				if (ExpectedQTLMutationsPerChromosome_P < 6)
					number_mutations = poissonrand(ExpectedQTLMutationsPerChromosome_P);
				else
					number_mutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome_P, sqrtEQTLMPC_P);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome_P);
					randbivnorm(mutational_std_dev_0, mutational_std_dev_1, MutationalCorrelation, drnum0, drnum1);
					progeny[i].MaternalChromosome[j].QTLeffect_P0[whichlocus] = progeny[i].MaternalChromosome[j].QTLeffect_P0[whichlocus] + drnum0;
					progeny[i].MaternalChromosome[j].QTLeffect_P1[whichlocus] = progeny[i].MaternalChromosome[j].QTLeffect_P1[whichlocus] + drnum1;
				}

				// Paternal Chromosome
				if (ExpectedQTLMutationsPerChromosome_P < 6)
					number_mutations = poissonrand(ExpectedQTLMutationsPerChromosome_P);
				else
					number_mutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome_P, sqrtEQTLMPC_P);
				for (k = 0; k < number_mutations; k++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome_P);
					randbivnorm(mutational_std_dev_0, mutational_std_dev_1, MutationalCorrelation, drnum0, drnum1);
					progeny[i].PaternalChromosome[j].QTLeffect_P0[whichlocus] = progeny[i].PaternalChromosome[j].QTLeffect_P0[whichlocus] + drnum0;
					progeny[i].PaternalChromosome[j].QTLeffect_P1[whichlocus] = progeny[i].PaternalChromosome[j].QTLeffect_P1[whichlocus] + drnum1;
				}

			}

			progeny[i].calculate_genotypic_values(NumberChromosomes, NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P, epistasis_allowed, ep_par_collection);
			progeny[i].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
			progeny[i].MatingSuccess = 0;
			progeny[i].Alive = true;

		} // end of i loop

	}

	void natural_selection(bool init_gens)
	{
		int i;
		double survival_prob, dRnum1;
		double dSu, dSv, dSc;
		double SSsqrt[2];
		double local_sel_str[2], local_sel_corr, local_opt[2];

		if (init_gens)
		{
			for (i = 0; i < 2; i++)
			{
				local_sel_str[i] = InitialSelectionStrength[i];
				local_opt[i] = InitialOptimum[i];
			}
			local_sel_corr = InitialSelectionalCorrelation;
		}
		else
		{
			for (i = 0; i < 2; i++)
			{
				local_sel_str[i] = SelectionStrength[i];
				local_opt[i] = Optimum[i];
			}
			local_sel_corr = SelectionalCorrelation;
		}

		dSc = 2 * (1 - local_sel_corr * local_sel_corr);
		SSsqrt[0] = sqrt(local_sel_str[0]);
		SSsqrt[1] = sqrt(local_sel_str[1]);

		for (i = 0; i < Nprogeny; i++)
		{
			//Bivariate selection on traits 0 and 1:
			if (SSsqrt[0] > 0 && SSsqrt[1] > 0)
			{
				dSu = (progeny[i].Phenotype[0] - local_opt[0]) / SSsqrt[0];
				dSv = (progeny[i].Phenotype[1] - local_opt[1]) / SSsqrt[1];
				survival_prob = exp((2 * local_sel_corr*dSu*dSv - dSu * dSu - dSv * dSv) / dSc);
			}

			// Selection on trait 0 only
			if (SSsqrt[0] > 0 && SSsqrt[1] == 0)
			{
				survival_prob = exp(-1.0 *
					(progeny[i].Phenotype[0] - local_opt[0])*
					(progeny[i].Phenotype[0] - local_opt[0])
					/ (2 * local_sel_str[0]));
			}

			// Selection on trait 1 only
			if (SSsqrt[0] == 0 && SSsqrt[1] > 0)
			{
				survival_prob = exp(-1.0 *
					(progeny[i].Phenotype[1] - local_opt[1])*
					(progeny[i].Phenotype[1] - local_opt[1])
					/ (2 * local_sel_str[1]));
			}

			// No selection on either trait
			if (SSsqrt[0] == 0 && SSsqrt[1] == 0)
			{
				survival_prob = 1;
			}

			dRnum1 = genrand();
			if (dRnum1 < survival_prob)
				progeny[i].Alive = true;
			else
				progeny[i].Alive = false;
		} // end of i loop

	}

	void population_regulation()
	{
		int i;
		double carrying_capacity_unfilled, progeny_left, keep_prob;
		int number_adults_chosen;
		double drnd1;

		progeny_left = 0;
		for (i = 0; i < Nprogeny; i++)
			if (progeny[i].Alive)
				progeny_left++;

		carrying_capacity_unfilled = CarryingCapacity;
		number_adults_chosen = 0;
		for (i = 0; i < Nprogeny; i++)
		{
			if (progeny[i].Alive)
			{
				keep_prob = carrying_capacity_unfilled / progeny_left;
				drnd1 = genrand();
				if (drnd1 < keep_prob)
				{
					deep_copy_individuals(adult[number_adults_chosen], progeny[i]);
					carrying_capacity_unfilled = carrying_capacity_unfilled - 1;
					number_adults_chosen++;
				}
				progeny_left = progeny_left - 1;
			} // end of if (progeny[i].Alive)
		} // end of i
		PopulationSize = number_adults_chosen;
	}

	void calculate_values_progeny()
	{
		int i;
		double dNP = Nprogeny;
		for (i = 0; i < 2; i++)
		{
			phenotypic_mean[i] = 0;
			genotypic_mean[i] = 0;
			phenotypic_variance[i] = 0;
			genotypic_variance[i] = 0;
		}
		phenotypic_covariance = 0;
		genotypic_covariance = 0;

		if (dNP > 0)
		{
			// Calculate the Means
			for (i = 0; i < Nprogeny; i++)
			{
				phenotypic_mean[0] = phenotypic_mean[0] + progeny[i].Phenotype[0];
				phenotypic_mean[1] = phenotypic_mean[1] + progeny[i].Phenotype[1];
				genotypic_mean[0] = genotypic_mean[0] + progeny[i].Genotype[0];
				genotypic_mean[1] = genotypic_mean[1] + progeny[i].Genotype[1];
			}
			phenotypic_mean[0] = phenotypic_mean[0] / dNP;
			phenotypic_mean[1] = phenotypic_mean[1] / dNP;
			genotypic_mean[0] = genotypic_mean[0] / dNP;
			genotypic_mean[1] = genotypic_mean[1] / dNP;

			// Calculate Variances and Covariances
			for (i = 0; i < Nprogeny; i++)
			{
				phenotypic_variance[0] = phenotypic_variance[0]
					+ (progeny[i].Phenotype[0] - phenotypic_mean[0])
					* (progeny[i].Phenotype[0] - phenotypic_mean[0]);
				phenotypic_variance[1] = phenotypic_variance[1]
					+ (progeny[i].Phenotype[1] - phenotypic_mean[1])
					* (progeny[i].Phenotype[1] - phenotypic_mean[1]);
				phenotypic_covariance = phenotypic_covariance
					+ (progeny[i].Phenotype[0] - phenotypic_mean[0])
					* (progeny[i].Phenotype[1] - phenotypic_mean[1]);

				genotypic_variance[0] = genotypic_variance[0]
					+ (progeny[i].Genotype[0] - genotypic_mean[0])
					* (progeny[i].Genotype[0] - genotypic_mean[0]);
				genotypic_variance[1] = genotypic_variance[1]
					+ (progeny[i].Genotype[1] - genotypic_mean[1])
					* (progeny[i].Genotype[1] - genotypic_mean[1]);
				genotypic_covariance = genotypic_covariance
					+ (progeny[i].Genotype[0] - genotypic_mean[0])
					* (progeny[i].Genotype[1] - genotypic_mean[1]);
			} // end of i
			phenotypic_variance[0] = phenotypic_variance[0] / dNP;
			phenotypic_variance[1] = phenotypic_variance[1] / dNP;
			genotypic_variance[0] = genotypic_variance[0] / dNP;
			genotypic_variance[1] = genotypic_variance[1] / dNP;
			phenotypic_covariance = phenotypic_covariance / dNP;
			genotypic_covariance = genotypic_covariance / dNP;
		}

		// Calculate the phenotypic and genotypic correlations
		double dtemp;
		dtemp = sqrt(phenotypic_variance[0] * phenotypic_variance[1]);
		if (dtemp > 0)
			phenotypic_correlation = phenotypic_covariance / dtemp;
		else
			phenotypic_correlation = 0;

		dtemp = sqrt(genotypic_variance[0] * genotypic_variance[1]);
		if (dtemp > 0)
			genotypic_correlation = genotypic_covariance / dtemp;
		else
			genotypic_correlation = 0;

		// Check to make sure the variances are both non-zero before
		// attempting to calculate the G-matrix.

		if (genotypic_variance[0] > 0 && genotypic_variance[1] > 0)
		{

			// Calculate the eigenvectors and eigenvalues of the G-matrix

			double a, b, c, d;
			double ev[2];
			a = genotypic_variance[0];
			b = genotypic_covariance;
			c = genotypic_covariance;
			d = genotypic_variance[1];
			const double m_pi = 3.14159265358979323846;

			// Use the formula for a 2x2 matrix to calculate eigenvalues
			if ((a + d)*(a + d) > 4 * (a*d - b * c))
			{
				ev[0] = ((a + d) + sqrt((a + d)*(a + d) - 4 * (a*d - b * c))) / 2;
				ev[1] = ((a + d) - sqrt((a + d)*(a + d) - 4 * (a*d - b * c))) / 2;
			}
			else
			{
				ev[0] = 0;
				ev[1] = 0;
			}

			// Make sure EigenValue[0] is the larger of the two
			if (ev[0] > ev[1])
			{
				EigenValue[0] = ev[0];
				EigenValue[1] = ev[1];
			}
			else
			{
				EigenValue[0] = ev[1];
				EigenValue[1] = ev[0];
			}

			// Calculate eigenvectors using 2x2 matrix formulae
			if (c != 0)
			{
				EigenVector1[0] = (EigenValue[0] - d) / sqrt((EigenValue[0] - d)*(EigenValue[0] - d) + c * c);
				EigenVector1[1] = c / sqrt((EigenValue[0] - d)*(EigenValue[0] - d) + c * c);
				EigenVector2[0] = (EigenValue[1] - d) / sqrt((EigenValue[1] - d)*(EigenValue[1] - d) + c * c);
				EigenVector2[1] = c / sqrt((EigenValue[1] - d)*(EigenValue[1] - d) + c * c);
			}
			if (c == 0)
			{
				if (d > a)
				{
					EigenVector1[0] = 0;
					EigenVector1[1] = 1;
					EigenVector2[0] = 1;
					EigenVector2[1] = 0;
				}
				else
				{
					EigenVector1[0] = 1;
					EigenVector1[1] = 0;
					EigenVector2[0] = 0;
					EigenVector2[1] = 1;
				}
			}
			// Calculate the angle of the leading eigenvector in degrees
			if (EigenVector1[0] != 0)
			{
				LeadAngle = atan(EigenVector1[1] / EigenVector1[0])*(180 / m_pi);
			}
			else
			{
				LeadAngle = 90;
			}

			Sigma = EigenValue[0] + EigenValue[1];
			Epsilon = EigenValue[1] / EigenValue[0];
		}
		else // deal with the case where one or both genetic variances are zero
		{
			// If either variance is zero, then the covariance is also zero.
			// If the covariance is zero, then the eigenvalues are the same as the variances.
			EigenValue[0] = genotypic_variance[0];
			EigenValue[1] = genotypic_variance[1];
			Sigma = EigenValue[0] + EigenValue[1];
			Epsilon = 0; // If either variance is zero, so is epsilon

						 // if only one variance is zero, we will define the angle to 
						 // point along the axis of the trait with the non-zero variance.

			EigenVector1[0] = 0;
			EigenVector1[1] = 0;
			EigenVector2[0] = 0;
			EigenVector2[1] = 0;
			LeadAngle = 0;

			if (EigenValue[0] > 0)
			{
				EigenVector1[0] = 1;
				EigenVector1[1] = 0;
				EigenVector2[0] = 0;
				EigenVector2[1] = 1;
			}
			if (EigenValue[1] > 0)
			{
				EigenVector1[0] = 0;
				EigenVector1[1] = 1;
				EigenVector2[0] = 1;
				EigenVector2[1] = 0;
				LeadAngle = 90;
			}
		}

		// Calculate selection differentials on trait 0 and trait 1
		// These selection differentials represent viability 
		// selection, because they include only progeny survival
		// during the juvenile phase of the lifecycle.

		// The selection differential is the covariance between
		// trait values and fitness. Here, fitness is based on 
		// whether or not the individual survived viability selection.

		// We use some of the variables calculated above:
		// phenotypic_mean[2], dNP (number of progeny).

		double mean_fitness;
		if (dNP > 0)
		{
			// Tally the mean fitness of the offspring
			mean_fitness = 0;
			for (i = 0; i < Nprogeny; i++)
			{
				if (progeny[i].Alive)
				{
					mean_fitness++;
				}
			}
			mean_fitness = mean_fitness / dNP;
		}

		if (dNP > 0 && mean_fitness > 0 && mean_fitness < 1)
		{
			sel_diff_trt_0 = 0;
			sel_diff_trt_1 = 0;
			for (i = 0; i < Nprogeny; i++)
			{
				if (progeny[i].Alive) // Alive and dead have different fitnesses (1 vs. 0)
				{
					sel_diff_trt_0 = sel_diff_trt_0 + (progeny[i].Phenotype[0] 
						- phenotypic_mean[0]) * (1.0 / mean_fitness - 1.0);
					sel_diff_trt_1 = sel_diff_trt_1 + (progeny[i].Phenotype[1] 
						- phenotypic_mean[1]) * (1.0 / mean_fitness - 1.0);
				}
				else
				{
					sel_diff_trt_0 = sel_diff_trt_0 + (progeny[i].Phenotype[0] 
						- phenotypic_mean[0]) * (-1.0);
					sel_diff_trt_1 = sel_diff_trt_1 + (progeny[i].Phenotype[1] 
						- phenotypic_mean[1]) * (-1.0);
				}
			}
			sel_diff_trt_0 = sel_diff_trt_0 / dNP;
			sel_diff_trt_1 = sel_diff_trt_1 / dNP;
		}
		else
		{
			sel_diff_trt_0 = 0;
			sel_diff_trt_1 = 0;
		}

		cPrevEval[0] = fabs(PrevEval[0] - EigenValue[0]);
		cPrevEval[1] = fabs(PrevEval[1] - EigenValue[1]);
		cPrevAngle = fabs(PrevAngle - LeadAngle);
		cPrevSigma = fabs(PrevSigma - Sigma);
		cPrevEpsilon = fabs(PrevEpsilon - Epsilon);
		cPrevG00 = fabs(PrevG00 - genotypic_variance[0]);
		cPrevG11 = fabs(PrevG11 - genotypic_variance[1]);
		cPrevG01 = fabs(PrevG01 - genotypic_covariance);
		cPrevRg = fabs(PrevRg - genotypic_correlation);

		// Ensure that we are calculating the smallest 
		// possible change in eigenvector angle.
		if (cPrevAngle > 90)
			cPrevAngle = 180 - cPrevAngle;

		PrevEval[0] = EigenValue[0];
		PrevEval[1] = EigenValue[1];
		PrevAngle = LeadAngle;
		PrevSigma = Sigma;
		PrevEpsilon = Epsilon;
		PrevG00 = genotypic_variance[0];
		PrevG11 = genotypic_variance[1];
		PrevG01 = genotypic_covariance;
		PrevRg = genotypic_correlation;

	}

	void output_population_variables(int generation)
	{
		if (generation == 0) // Output the header in generation zero
		{
			std::cout << "\nGen\tzbar0\tzbar1\tP00\tP11\tP12\tr(P)\tgbar0\tgbar1\tG00\tG11\tG01\tr(G)";
		}

		std::cout << "\n" << generation;
		std::cout << std::setprecision(3) << std::fixed;
		std::cout << "\t" << phenotypic_mean[0];
		std::cout << "\t" << phenotypic_mean[1];
		std::cout << "\t" << phenotypic_variance[0];
		std::cout << "\t" << phenotypic_variance[1];
		std::cout << "\t" << phenotypic_covariance;
		std::cout << "\t" << phenotypic_correlation;
		std::cout << "\t" << genotypic_mean[0];
		std::cout << "\t" << genotypic_mean[1];
		std::cout << "\t" << genotypic_variance[0];
		std::cout << "\t" << genotypic_variance[1];
		std::cout << "\t" << genotypic_covariance;
		std::cout << "\t" << genotypic_correlation;
	}

	int getNumberOfGenerations()
	{
		return NumberOfGenerations;
	}

	int getNumberInterveningGens()
	{
		return NumberOfInterveningGens;
	}

	void save_population_variables(int generation)
	{
		std::ofstream outfile;
		if (generation == 0) // Output the header in generation zero
		{
			outfile.open("output.csv");
			outfile << "Gen,N,zbar0,zbar1,P00,P11,P12,r(P),gbar0,gbar1,G00,G11,G01,r(G)";
			outfile << ",Lambda1,Lambda2,EvecX,EvecY,Angle,Size,Eccen,Strt0,Strt1";
			outfile << ",~Lmbd1,~Lmbd2,~Ang,~Size,~Ecc,~G00,~G11,~G01,~r(g)";
			outfile << ",ASR,Im,If,MdifM,MdifF";
			outfile.close();
		}

		outfile.open("output.csv", std::fstream::app);
		outfile << "\n" << generation;
		outfile << "," << PopulationSize;
		outfile << "," << phenotypic_mean[0];
		outfile << "," << phenotypic_mean[1];
		outfile << "," << phenotypic_variance[0];
		outfile << "," << phenotypic_variance[1];
		outfile << "," << phenotypic_covariance;
		outfile << "," << phenotypic_correlation;
		outfile << "," << genotypic_mean[0];
		outfile << "," << genotypic_mean[1];
		outfile << "," << genotypic_variance[0];
		outfile << "," << genotypic_variance[1];
		outfile << "," << genotypic_covariance;
		outfile << "," << genotypic_correlation;
		outfile << "," << EigenValue[0];
		outfile << "," << EigenValue[1];
		outfile << "," << EigenVector1[0];
		outfile << "," << EigenVector1[1];
		outfile << "," << LeadAngle;
		outfile << "," << Sigma;
		outfile << "," << Epsilon;
		outfile << "," << sel_diff_trt_0;
		outfile << "," << sel_diff_trt_1;
		outfile << "," << cPrevEval[0];
		outfile << "," << cPrevEval[1];
		outfile << "," << cPrevAngle;
		outfile << "," << cPrevSigma;
		outfile << "," << cPrevEpsilon;
		outfile << "," << cPrevG00;
		outfile << "," << cPrevG11;
		outfile << "," << cPrevG01;
		outfile << "," << cPrevRg;
		outfile << "," << SexRatio;
		outfile << "," << Im;
		outfile << "," << If;
		outfile << "," << MdiffMales;
		outfile << "," << MdiffFemales;
		outfile.close();
	}

	void gaussian_mating()
	{
		// This function implements gaussian mate choice in a
		// polygynous mating system. It should be used instead
		// of other mating functions (like polygynous_mating).
		// In this function, each female mates at most once, and
		// each male can mate an unlimited number of times, so
		// the mating system is polygynous.

		// Females choose males based on their trait values.
		// Mating preferences are Gaussian in the sense that each female
		// has an ideal preferred male phenotype and her preferences
		// fall off as males depart from her preferred phenotype.
		// The drop in mating probability as the male departs from
		// the preferred phenotype is modeled as a Gaussian-shaped
		// function.

		// Absolute or relative preferences? This function uses relative
		// preferences, and they are relative to the phenotypic mean of 
		// the males. Consequently, if a female's preference trait has
		// a value of 0.72, then that particular female's ideal male
		// has a trait value 0.72 units greater than the male mean trait
		// value. If she encountered such a male, she would mate with
		// probability 1. Her probability of mating would drop off for
		// males with trait values larger or smaller than the preferred
		// value.

		int i, m, n;
		int iPC;
		bool mate_found;
		int mateID, counter, rnum;
		double dRnum;
		double mate_prob;
		double dNmales;
		double mean_male_trait0;
		double ideal;
		int iNmales;
		int iMaleIndex;
		int *iMaleList = new int[PopulationSize];

		// Check to make sure at least one male is present in the population.
		// Also calculate the mean of the male ornament trait (trait 0). 

		bool males_present = false;
		iNmales = 0;
		mean_male_trait0 = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (!adult[i].Female)
			{
				males_present = true;
				iMaleList[iNmales] = i;
				iNmales++;
				mean_male_trait0 = mean_male_trait0 + adult[i].Phenotype[0];
			}
		}
		dNmales = iNmales;

		if (dNmales > 0)
		{
			mean_male_trait0 = mean_male_trait0 / dNmales;
		}
		else
		{
			mean_male_trait0 = 0;
			PopulationExtinct = true;
		}

		iPC = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female && males_present)
			{
				// Find a mate for this female
				// She only gets MaxMatingEncounters tries to find
				// someone. 

				// She only encounters males.

				// If she doesn't find a mate in the allotted number
				// of tries, then she produces no progeny.

				// If she doesn't find a mate in the allotted number
				// of tries, then she produces no progeny.

				// For each female, we first have to determine her
				// ideal mate phenotype. Trait 1 is the female preference
				// trait and it describes the deviation from the male
				// mean of her preferred mate at Trait 0, the ornament. 

				ideal = adult[i].Phenotype[1] + mean_male_trait0;

				mate_found = false;
				counter = 0;
				while (counter < MaxMatingEncounters && !mate_found)
				{
					rnum = randnum(iNmales);
					iMaleIndex = iMaleList[rnum];

					if (!adult[rnum].Female)
					{
						// Calculate the focal female's (adult[i]) probability 
						// of mating with this random male (adult[iMaleIndex]). The
						// probability is calculated using the Gaussian 
						// probability density function without the normalization
						// term.
						if (GaussianPreferenceVariance > 0)
						{
							mate_prob = exp(-1 * (adult[iMaleIndex].Phenotype[0] - ideal)*
								(adult[iMaleIndex].Phenotype[0] - ideal)
								/ (2 * GaussianPreferenceVariance));
						}
						else
						{
							mate_prob = 1;
						}

						// Generate a random number (roll the dice!)
						dRnum = genrand();

						if (dRnum < mate_prob)
						{
							mateID = iMaleIndex;
							mate_found = true;
						}
					}
					counter++;
				} // end of while

				// If a mate is found, produce progeny
				if (mate_found)
				{
					for (m = 0; m < Fecundity; m++)
					{
						if (iPC >= NprogMax)
							iPC = NprogMax - 1;

						for (n = 0; n < NumberChromosomes; n++)
						{
							ProduceRecombinedChromosome(progeny[iPC].MaternalChromosome[n], adult[i], n, ExpRecombPerChromosome);
							ProduceRecombinedChromosome(progeny[iPC].PaternalChromosome[n], adult[mateID], n, ExpRecombPerChromosome);
						}

						progeny[iPC].set_sex();
						progeny[iPC].calculate_genotypic_values(NumberChromosomes, NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P, epistasis_allowed, ep_par_collection);
						progeny[iPC].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
						progeny[iPC].Alive = true;
						iPC++;

					} // end of m loop
					adult[i].MatingSuccess++;
					adult[mateID].MatingSuccess++;
				} // end of if (mate_found)

			} // end of if (adult[i].Female && males_present)
		} // end of i loop
		Nprogeny = iPC;
		delete[] iMaleList;
	}

	bool is_extinct()
	{
		return PopulationExtinct;
	}

	void sex_limited_selection(bool init_gens)
	{
		int i;
		double survival_prob, dRnum1, optimum[2];
		double local_sel_str[2];

		if (init_gens)
		{
			local_sel_str[0] = InitialSelectionStrength[0];
			local_sel_str[1] = InitialSelectionStrength[1];
			optimum[0] = InitialOptimum[0];
			optimum[1] = InitialOptimum[1];
		}
		else
		{
			local_sel_str[0] = SelectionStrength[0];
			local_sel_str[1] = SelectionStrength[1];
			optimum[0] = Optimum[0];
			optimum[1] = Optimum[1];
		}

		for (i = 0; i < Nprogeny; i++)
		{
			// Calculate survival probabilities for males and females separately
			if (!progeny[i].Female)
			{
				if (local_sel_str[0] > 0)
				{
					survival_prob = exp(-1.0 *
						(progeny[i].Phenotype[0] - optimum[0])*
						(progeny[i].Phenotype[0] - optimum[0])
						/ (2 * local_sel_str[0]));
				}
				else
				{
					survival_prob = 1; // SelectionStrength = 0 means no selection.
				}
			}
			else
			{
				if (local_sel_str[1] > 0)
				{
					survival_prob = exp(-1.0 *
						(progeny[i].Phenotype[1] - optimum[1])*
						(progeny[i].Phenotype[1] - optimum[1])
						/ (2 * local_sel_str[1]));
				}
				else
				{
					survival_prob = 1; // SelectionStrength = 0 means no selection.
				}
			}

			// Generate a random number to see if the individual survives
			dRnum1 = genrand();
			if (dRnum1 < survival_prob)
				progeny[i].Alive = true;
			else
				progeny[i].Alive = false;
		} // end of i loop
	}

	int getNumberOfInitialGenerations()
	{
		return NumberOfInitialGenerations;
	}

	void calculate_values_adults()
	{
		int i;
		double MeanMSmales, MeanMSfemales, MeanTrait0males, MeanTrait1females;
		double VarMSmales, VarMSfemales, Nma, Nfe;

		Nma = 0;
		Nfe = 0;
		MeanMSmales = 0;
		MeanMSfemales = 0;
		MeanTrait0males = 0;
		MeanTrait1females = 0;

		// Calculate mean mating success and mean trait values
		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female)
			{
				Nfe++;
				MeanMSfemales = MeanMSfemales + adult[i].MatingSuccess;
				MeanTrait1females = MeanTrait1females + adult[i].Phenotype[1];
			}
			else
			{
				Nma++;
				MeanMSmales = MeanMSmales + adult[i].MatingSuccess;
				MeanTrait0males = MeanTrait0males + adult[i].Phenotype[0];
			}
		} // end of i

		if (Nma > 0)
		{
			MeanMSmales = MeanMSmales / Nma;
			MeanTrait0males = MeanTrait0males / Nma;
		}
		else
		{
			MeanMSmales = 0;
			MeanTrait0males = 0;
		}

		if (Nfe > 0)
		{
			MeanMSfemales = MeanMSfemales / Nfe;
			MeanTrait1females = MeanTrait1females / Nfe;
		}
		else {
			MeanMSfemales = 0;
			MeanTrait1females = 0;
		}

		// Calculate the Sex Ratio.
		if (Nma + Nfe > 0)
			SexRatio = Nma / (Nma + Nfe);
		else
			SexRatio = 0;

		// Calculate opportunities for sexual selection
		// and mating differentials

		VarMSmales = 0;
		VarMSfemales = 0;
		MdiffMales = 0;
		MdiffFemales = 0;

		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female)
			{
				VarMSfemales = VarMSfemales + (adult[i].MatingSuccess - MeanMSfemales)
					*(adult[i].MatingSuccess - MeanMSfemales);
				if (MeanMSfemales > 0)
					MdiffFemales = MdiffFemales + (adult[i].Phenotype[1] - MeanTrait1females)
					*(adult[i].MatingSuccess / MeanMSfemales - 1);
			}
			else
			{
				VarMSmales = VarMSmales + (adult[i].MatingSuccess - MeanMSmales)
					*(adult[i].MatingSuccess - MeanMSmales);
				if (MeanMSmales > 0)
					MdiffMales = MdiffMales + (adult[i].Phenotype[0] - MeanTrait0males)
					*(adult[i].MatingSuccess / MeanMSmales - 1);
			}
		} // end of i

		if (Nfe > 0)
		{
			VarMSfemales = VarMSfemales / Nfe;
			MdiffFemales = MdiffFemales / Nfe;
		}
		else
		{
			VarMSfemales = 0;
			MdiffFemales = 0;
		}

		if (Nma > 0)
		{
			VarMSmales = VarMSmales / Nma;
			MdiffMales = MdiffMales / Nma;
		}
		else
		{
			VarMSmales = 0;
			MdiffMales = 0;
		}

		if (MeanMSmales > 0)
			Im = VarMSmales / (MeanMSmales*MeanMSmales);
		else
			Im = 0;

		if (MeanMSfemales > 0)
			If = VarMSfemales / (MeanMSfemales*MeanMSfemales);
		else
			If = 0;


	}

	void selection(bool initial_generations)
	{
		if (initial_generations)
		{
			if (InitialSelectionSexLimited)
				sex_limited_selection(true);
			else
				natural_selection(true);
		}
		else
		{
			if (ExperimentalSelectionSexLimited)
				sex_limited_selection(false);
			else
				natural_selection(false);
		}
	}

	void monogamy()
	{
		// This function implements strict monogamy.
		// Under this mating system, each female mates once,
		// and each male also mates once. Sex is assigned
		// on the fly, and the sex ratio is equalized.
		// If there's an uneven number of individuals,
		// one individual does not mate.

		int i, m, n;
		int iPC;
		int maleID, femaleID, rnum;
		int number_females;

		bool *already_chosen = new bool[PopulationSize];
		int *male_list = new int[PopulationSize];
		int *female_list = new int[PopulationSize];

		for (i = 0; i < PopulationSize; i++)
		{
			adult[i].Female = true;
			already_chosen[i] = false;
		}

		// We need to make half the population male
		// If the population contains an odd number of
		// individuals, then we need one more male than
		// female, so below when all the females mate
		// they do not run out of males.

		int number_males = PopulationSize / 2;
		if (PopulationSize % 2 != 0)
			number_males++;

		// Now make a list of unique males

		for (i = 0; i < number_males; i++)
		{
			rnum = randnum(PopulationSize);
			while (already_chosen[rnum])
				rnum = randnum(PopulationSize);
			adult[rnum].Female = false;
			already_chosen[rnum] = true;
			male_list[i] = rnum;
		}

		// Make a list of females

		number_females = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female)
			{
				female_list[number_females] = i;
				number_females++;
			}
		}

		iPC = 0;
		for (i = 0; i < number_females; i++)
		{

			maleID = male_list[i];
			femaleID = female_list[i];

			// Produce progeny
			// If a mate is found, produce progeny
		
			for (m = 0; m < Fecundity; m++)
			{
				if (iPC >= NprogMax)
					iPC = NprogMax - 1;

				for (n = 0; n < NumberChromosomes; n++)
				{
					ProduceRecombinedChromosome(progeny[iPC].MaternalChromosome[n], adult[femaleID], n, ExpRecombPerChromosome);
					ProduceRecombinedChromosome(progeny[iPC].PaternalChromosome[n], adult[maleID], n, ExpRecombPerChromosome);
				}

				progeny[iPC].set_sex();
				progeny[iPC].calculate_genotypic_values(NumberChromosomes, NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P, epistasis_allowed, ep_par_collection);
				progeny[iPC].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
				iPC++;
			} // end of m loop
			adult[femaleID].MatingSuccess++;
			adult[maleID].MatingSuccess++;
		} // end of i loop
		Nprogeny = iPC;

		delete[] already_chosen;
		delete[] male_list;
		delete[] female_list;
	}

	void store_variables_in_memory()
	{
		sim_data.vN.push_back(PopulationSize);
		sim_data.vZbar0.push_back(phenotypic_mean[0]);
		sim_data.vZbar1.push_back(phenotypic_mean[1]);
		sim_data.vP00.push_back(phenotypic_variance[0]);
		sim_data.vP11.push_back(phenotypic_variance[1]);
		sim_data.vP01.push_back(phenotypic_covariance);
		sim_data.vRp.push_back(phenotypic_correlation);
		sim_data.vGbar0.push_back(genotypic_mean[0]);
		sim_data.vGbar1.push_back(genotypic_mean[1]);
		sim_data.vG00.push_back(genotypic_variance[0]);
		sim_data.vG11.push_back(genotypic_variance[1]);
		sim_data.vG01.push_back(genotypic_covariance);
		sim_data.vRg.push_back(genotypic_correlation);
		sim_data.vLambda1.push_back(EigenValue[0]);
		sim_data.vLambda2.push_back(EigenValue[1]);
		sim_data.vEvecX.push_back(EigenVector1[0]);
		sim_data.vEvecY.push_back(EigenVector1[1]);
		sim_data.vAngle.push_back(LeadAngle);
		sim_data.vSize.push_back(Sigma);
		sim_data.vEccen.push_back(Epsilon);
		sim_data.vStrt0.push_back(sel_diff_trt_0);
		sim_data.vStrt1.push_back(sel_diff_trt_1);
		sim_data.vcLmbd1.push_back(cPrevEval[0]);
		sim_data.vcLmbd2.push_back(cPrevEval[1]);
		sim_data.vcAng.push_back(cPrevAngle);
		sim_data.vcSize.push_back(cPrevSigma);
		sim_data.vcEcc.push_back(cPrevEpsilon);
		sim_data.vcG00.push_back(cPrevG00);
		sim_data.vcG11.push_back(cPrevG11);
		sim_data.vcG01.push_back(cPrevG01);
		sim_data.vcRg.push_back(cPrevRg);
		sim_data.vASR.push_back(SexRatio);
		sim_data.vIm.push_back(Im);
		sim_data.vIf.push_back(If);
		sim_data.vMdifM.push_back(MdiffMales);
		sim_data.vMdifF.push_back(MdiffFemales);
		sim_data.vV0.push_back(est_va00);
		sim_data.vV1.push_back(est_va11);
		sim_data.vV01.push_back(est_va01);
	}

	void save_stored_variables(std::string filename)
	{
		size_t i;

		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);

		outfile << "Gen,N,zbar0,zbar1,P00,P11,P12,r(P),gbar0,gbar1,G00,G11,G01,r(G)";
		outfile << ",Lambda1,Lambda2,EvecX,EvecY,Angle,Size,Eccen,Strt0,Strt1";
		outfile << ",~Lmbd1,~Lmbd2,~Ang,~Size,~Ecc,~G00,~G11,~G01,~r(g)";
		outfile << ",ASR,Im,If,MdifM,MdifF,estG00,estG11,estG01";

		for (i = 0; i < sim_data.vN.size(); i++)
		{
			outfile << "\n" << i;
			outfile << "," << sim_data.vN[i];
			outfile << "," << sim_data.vZbar0[i];
			outfile << "," << sim_data.vZbar1[i];
			outfile << "," << sim_data.vP00[i];
			outfile << "," << sim_data.vP11[i];
			outfile << "," << sim_data.vP01[i];
			outfile << "," << sim_data.vRp[i];
			outfile << "," << sim_data.vGbar0[i];
			outfile << "," << sim_data.vGbar1[i];
			outfile << "," << sim_data.vG00[i];
			outfile << "," << sim_data.vG11[i];
			outfile << "," << sim_data.vG01[i];
			outfile << "," << sim_data.vRg[i];
			outfile << "," << sim_data.vLambda1[i];
			outfile << "," << sim_data.vLambda2[i];
			outfile << "," << sim_data.vEvecX[i];
			outfile << "," << sim_data.vEvecY[i];
			outfile << "," << sim_data.vAngle[i];
			outfile << "," << sim_data.vSize[i];
			outfile << "," << sim_data.vEccen[i];
			outfile << "," << sim_data.vStrt0[i];
			outfile << "," << sim_data.vStrt1[i];
			outfile << "," << sim_data.vcLmbd1[i];
			outfile << "," << sim_data.vcLmbd2[i];
			outfile << "," << sim_data.vcAng[i];
			outfile << "," << sim_data.vcSize[i];
			outfile << "," << sim_data.vcEcc[i];
			outfile << "," << sim_data.vcG00[i];
			outfile << "," << sim_data.vcG11[i];
			outfile << "," << sim_data.vcG01[i];
			outfile << "," << sim_data.vcRg[i];
			outfile << "," << sim_data.vASR[i];
			outfile << "," << sim_data.vIm[i];
			outfile << "," << sim_data.vIf[i];
			outfile << "," << sim_data.vMdifM[i];
			outfile << "," << sim_data.vMdifF[i];
			outfile << "," << sim_data.vV0[i];
			outfile << "," << sim_data.vV1[i];
			outfile << "," << sim_data.vV01[i];
		}
		outfile.close();
	}

	void calc_run_means()
	{
		sim_data.calculate_means();
	}

	void append_means_to_output_file(std::string filename)
	{
		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn, std::fstream::app);
		outfile << "\n" << "Mean:";
		outfile << "," << sim_data.meanN;
		outfile << "," << sim_data.meanZbar0;
		outfile << "," << sim_data.meanZbar1;
		outfile << "," << sim_data.meanP00;
		outfile << "," << sim_data.meanP11;
		outfile << "," << sim_data.meanP01;
		outfile << "," << sim_data.meanRp;
		outfile << "," << sim_data.meanGbar0;
		outfile << "," << sim_data.meanGbar1;
		outfile << "," << sim_data.meanG00;
		outfile << "," << sim_data.meanG11;
		outfile << "," << sim_data.meanG01;
		outfile << "," << sim_data.meanRg;
		outfile << "," << sim_data.meanLambda1;
		outfile << "," << sim_data.meanLambda2;
		outfile << "," << sim_data.meanEvecX;
		outfile << "," << sim_data.meanEvecY;
		outfile << "," << sim_data.meanAngle;
		outfile << "," << sim_data.meanSize;
		outfile << "," << sim_data.meanEccen;
		outfile << "," << sim_data.meanStrt0;
		outfile << "," << sim_data.meanStrt1;
		outfile << "," << sim_data.meancLmbd1;
		outfile << "," << sim_data.meancLmbd2;
		outfile << "," << sim_data.meancAng;
		outfile << "," << sim_data.meancSize;
		outfile << "," << sim_data.meancEcc;
		outfile << "," << sim_data.meancG00;
		outfile << "," << sim_data.meancG11;
		outfile << "," << sim_data.meancG01;
		outfile << "," << sim_data.meancRg;
		outfile << "," << sim_data.meanASR;
		outfile << "," << sim_data.meanIm;
		outfile << "," << sim_data.meanIf;
		outfile << "," << sim_data.meanMdifM;
		outfile << "," << sim_data.meanMdifF;
		outfile << "," << sim_data.meanV0;
		outfile << "," << sim_data.meanV1;
		outfile << "," << sim_data.meanV01;
		outfile.close();
	}

	mean_recorder report_means()
	{
		m_rec.mean_list.push_back(sim_data.meanN);
		m_rec.mean_list.push_back(sim_data.meanZbar0);
		m_rec.mean_list.push_back(sim_data.meanZbar1);
		m_rec.mean_list.push_back(sim_data.meanP00);
		m_rec.mean_list.push_back(sim_data.meanP11);
		m_rec.mean_list.push_back(sim_data.meanP01);
		m_rec.mean_list.push_back(sim_data.meanRp);
		m_rec.mean_list.push_back(sim_data.meanGbar0);
		m_rec.mean_list.push_back(sim_data.meanGbar1);
		m_rec.mean_list.push_back(sim_data.meanG00);
		m_rec.mean_list.push_back(sim_data.meanG11);
		m_rec.mean_list.push_back(sim_data.meanG01);
		m_rec.mean_list.push_back(sim_data.meanRg);
		m_rec.mean_list.push_back(sim_data.meanLambda1);
		m_rec.mean_list.push_back(sim_data.meanLambda2);
		m_rec.mean_list.push_back(sim_data.meanEvecX);
		m_rec.mean_list.push_back(sim_data.meanEvecY);
		m_rec.mean_list.push_back(sim_data.meanAngle);
		m_rec.mean_list.push_back(sim_data.meanSize);
		m_rec.mean_list.push_back(sim_data.meanEccen);
		m_rec.mean_list.push_back(sim_data.meanStrt0);
		m_rec.mean_list.push_back(sim_data.meanStrt1);
		m_rec.mean_list.push_back(sim_data.meancLmbd1);
		m_rec.mean_list.push_back(sim_data.meancLmbd2);
		m_rec.mean_list.push_back(sim_data.meancAng);
		m_rec.mean_list.push_back(sim_data.meancSize);
		m_rec.mean_list.push_back(sim_data.meancEcc);
		m_rec.mean_list.push_back(sim_data.meancG00);
		m_rec.mean_list.push_back(sim_data.meancG11);
		m_rec.mean_list.push_back(sim_data.meancG01);
		m_rec.mean_list.push_back(sim_data.meancRg);
		m_rec.mean_list.push_back(sim_data.meanASR);
		m_rec.mean_list.push_back(sim_data.meanIm);
		m_rec.mean_list.push_back(sim_data.meanIf);
		m_rec.mean_list.push_back(sim_data.meanMdifM);
		m_rec.mean_list.push_back(sim_data.meanMdifF);
		m_rec.mean_list.push_back(sim_data.meanV0);
		m_rec.mean_list.push_back(sim_data.meanV1);
		m_rec.mean_list.push_back(sim_data.meanV01);
		return m_rec;
	}

	void save_parameter_values(std::string filename, double migr_rate)
	{
		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);
		outfile << "Parameter_Values:\n";

		// Demographic Parameters
		outfile << "Demographic_Parameters:\n";
		outfile << "No_Generations:," << NumberOfGenerations << "\n";
		outfile << "Initial_Pop_Size:," << PopulationSize << "\n";
		outfile << "Carrying_Capacity:," << CarryingCapacity << "\n";
		outfile << "Female_Fecundity:," << Fecundity << "\n";
		outfile << "Migration_Rate:," << migr_rate << "\n";
		outfile << "Adult_Sample_Size:," << SampleSizeAdults << "\n";

		// Mating Parameters
		outfile << "Mating_Parameters:\n";
		outfile << "Max_Mating_Enc.:," << MaxMatingEncounters << "\n";
		outfile << "Gaussian_Pref_Var:," << GaussianPreferenceVariance << "\n";

		// Quantitative Genetic Parameters
		outfile << "Quantitative_Genetic_Parameters:\n";
		outfile << "Number_Chromosomes:," << NumberChromosomes << "\n";
		outfile << "QTL_per_Chr_Trt0:," << NumberQTLsPerChromosome << "\n";
		outfile << "QTL_per_Chr_Trt1:," << NumberQTLsPerChromosome_1 << "\n";
		outfile << "QTL_per_Chr_Pleio:," << NumberQTLsPerChromosome_P << "\n";
		outfile << "Env_Variance_Trt0:," << EnvironmentalVariance[0] << "\n";
		outfile << "Env_Variance_Trt1:," << EnvironmentalVariance[1] << "\n";
		outfile << "Exp_Recomb_Per_Chr:," << ExpRecombPerChromosome << "\n";

		// Marker Loci
		outfile << "Marker_Loci_Parameters:\n";
		outfile << "Markers_per_Chr:," << NumberMarkerLociPerChromosome << "\n";
		outfile << "Marker_Mut_Rate:," << MutationRatePerMarker << "\n";
		outfile << "Max_Alleles_Per_Mark:," << MaxNumAllelesPerMarker << "\n";
		outfile << "Fst_Weighting_Variance:," << fst_weighting_variance << "\n";

		// Mutational Parameters
		outfile << "Mutational_Parameters:\n";
		outfile << "Mut_Var_Trait0:," << MutationalVariance[0] << "\n";
		outfile << "Mut_Var_Trait1:," << MutationalVariance[1] << "\n";
		outfile << "Mut_Correlation:," << MutationalCorrelation << "\n";
		outfile << "Mutation_Rate:," << MutationRatePerLocus << "\n";

		// Selection Parameters
		outfile << "Selection_Parameters:\n";
		outfile << "Omega_Trait0:," << SelectionStrength[0] << "\n";
		outfile << "Omega_Trait1:," << SelectionStrength[1] << "\n";
		outfile << "Selection_Corr:," << SelectionalCorrelation << "\n";
		outfile << "Optimum_Trait0:," << Optimum[0] << "\n";
		outfile << "Optimum_Trait1:," << Optimum[1] << "\n";
		if (ExperimentalSelectionSexLimited)
			outfile << "Sex_Limited_Sel:,true\n";
		else
			outfile << "Sex_Limited_Sel:,false\n";

		// Initial Generations
		outfile << "Initial_Generations_Parameters:\n";
		outfile << "No_Initial_Gens:," << NumberOfInitialGenerations << "\n";
		outfile << "No_Intervening_Gens:," << NumberOfInterveningGens << "\n";
		outfile << "Init_Omega_Trait0:," << InitialSelectionStrength[0] << "\n";
		outfile << "Init_Omega_Trait1:," << InitialSelectionStrength[1] << "\n";
		outfile << "Init_Sel_Corr:," << InitialSelectionalCorrelation << "\n";
		outfile << "Init_Opt_Trt0:," << InitialOptimum[0] << "\n";
		outfile << "Init_Opt_Trt1:," << InitialOptimum[1] << "\n";
		if (InitialSelectionSexLimited)
			outfile << "Init_Sex_Lim_Sel:,true\n";
		else
			outfile << "Init_Sex_Lim_Sel:,false\n";

		outfile << "Epistatic_Parameters:\n";
		outfile << "Epistasis_Allowed:,";
		if (epistasis_allowed)
			outfile << "true\n";
		else
			outfile << "false\n";
		outfile << "Nonpleiotropic_Epistatic_Variances:\n";
		outfile << "Ep_0by0_on_trt_0," << ep_var_list.ep_var_0by0_on_0 << "\n";
		outfile << "Ep_1by1_on_trt_1," << ep_var_list.ep_var_1by1_on_1 << "\n";
		outfile << "Ep_0by1_on_trt_0," << ep_var_list.ep_var_0by1_on_0 << "\n";
		outfile << "Ep_0by1_on_trt_1," << ep_var_list.ep_var_0by1_on_1 << "\n";
		outfile << "Ep_1by1_on_trt_0," << ep_var_list.ep_var_1by1_on_0 << "\n";
		outfile << "Ep_0by0_on_trt_1," << ep_var_list.ep_var_0by0_on_1 << "\n";

		outfile << "Pleiotropic_Epistatic_Variances:\n";
		outfile << "Pleio_0by0_on_trt_0," << ep_var_list.pleio_ep_var_0by0_on_0 << "\n";
		outfile << "Pleio_1by1_on_trt_1," << ep_var_list.pleio_ep_var_1by1_on_1 << "\n";
		outfile << "Pleio_0by1_on_trt_0," << ep_var_list.pleio_ep_var_0by1_on_0 << "\n";
		outfile << "Pleio_0by1_on_trt_1," << ep_var_list.pleio_ep_var_0by1_on_1 << "\n";
		outfile << "Pleio_1by1_on_trt_0," << ep_var_list.pleio_ep_var_1by1_on_0 << "\n";
		outfile << "Pleio_0by0_on_trt_1," << ep_var_list.pleio_ep_var_0by0_on_1 << "\n";

		outfile.close();

	}

	void update_parameter_values(parameter_value_set &parm_set, epistatic_parameter_collection e_p_c)
	{
		// Set Parameter Values

		// Initial Generations Parameters
		NumberOfInitialGenerations = parm_set.p_init_gens;
		NumberOfInterveningGens = parm_set.p_intervening_gens;
		InitialSelectionStrength[0] = parm_set.p_init_w00;
		InitialSelectionStrength[1] = parm_set.p_init_w11;
		InitialSelectionalCorrelation = parm_set.p_init_sel_corr;
		InitialOptimum[0] = parm_set.p_init_opt0;
		InitialOptimum[1] = parm_set.p_init_opt1;
		InitialSelectionSexLimited = parm_set.p_init_sex_lim;
		StartingQTLAllelicEffectStdDev = parm_set.p_starting_QTL_stdev;

		// Demographic Parameters
		NumberOfGenerations = parm_set.p_gens;
		CarryingCapacity = parm_set.p_carry_cap;
		PopulationSize = CarryingCapacity;
		Fecundity = parm_set.p_fecund;
		PopulationExtinct = false;

		// Mating Parameters
		MaxMatingEncounters = parm_set.p_mate_enc;
		GaussianPreferenceVariance = parm_set.p_pref_var;

		// Quantitative Genetic Parameters
		NumberChromosomes = parm_set.p_number_chromosomes;
		NumberQTLsPerChromosome = parm_set.p_loci_trt0;
		NumberQTLsPerChromosome_1 = parm_set.p_loci_trt1;
		NumberQTLsPerChromosome_P = parm_set.p_loci_pleio;
		EnvironmentalVariance[0] = parm_set.p_env_var0;
		EnvironmentalVariance[1] = parm_set.p_env_var1;
		EnvironmentalStDev[0] = sqrt(EnvironmentalVariance[0]);
		EnvironmentalStDev[1] = sqrt(EnvironmentalVariance[1]);
		ExpRecombPerChromosome = parm_set.p_exp_recomb_rate;

		// Mutational Parameters
		MutationalVariance[0] = parm_set.p_mut_var0;
		MutationalVariance[1] = parm_set.p_mut_var1;
		MutationalCorrelation = parm_set.p_mut_corr;
		MutationRatePerLocus = parm_set.p_mut_rate;

		// Selection Parameters
		SelectionStrength[0] = parm_set.p_exp_w00;
		SelectionStrength[1] = parm_set.p_exp_w11;
		SelectionalCorrelation = parm_set.p_exp_sel_corr;
		Optimum[0] = parm_set.p_exp_opt0;
		Optimum[1] = parm_set.p_exp_opt1;
		ExperimentalSelectionSexLimited = parm_set.p_exp_sex_lim;

		// Marker Loci
		NumberMarkerLociPerChromosome = parm_set.p_markers_per_chromosome;
		MutationRatePerMarker = parm_set.p_marker_mut_rate;
		MaxNumAllelesPerMarker = parm_set.p_max_alleles_per_marker;

		SampleSizeAdults = parm_set.p_sample_size;
		all_calc_interval = parm_set.p_all_freq_calc_interval;
		fst_weighting_variance = parm_set.p_fst_weighting_variance;

		// Epistatic parameters
		epistasis_allowed = parm_set.p_permit_epistasis;
		ep_par_collection = e_p_c;
		ep_var_list = parm_set.p_epistatic_vars;

		// Set some housekeeping variables
		double td;
		td = NumberQTLsPerChromosome;
		ExpectedQTLMutationsPerChromosome = MutationRatePerLocus * td;
		td = NumberQTLsPerChromosome_1;
		ExpectedQTLMutationsPerChromosome_1 = MutationRatePerLocus * td;
		td = NumberQTLsPerChromosome_P;
		ExpectedQTLMutationsPerChromosome_P = MutationRatePerLocus * td;
		sqrtEQTLMPC = sqrt(ExpectedQTLMutationsPerChromosome);
		sqrtEQTLMPC_1 = sqrt(ExpectedQTLMutationsPerChromosome_1);
		sqrtEQTLMPC_P = sqrt(ExpectedQTLMutationsPerChromosome_P);

		td = NumberMarkerLociPerChromosome;
		ExpectedMarkerMutationsPerChromosome = td * MutationRatePerMarker;
		sqrtEMMPC = sqrt(ExpectedMarkerMutationsPerChromosome);

	}

	int getNumberofProgenyAlive()
	{
		int prog_alive = 0;
		for (int i = 0; i < Nprogeny; i++)
		{
			if (progeny[i].Alive)
			{
				prog_alive++;
			}
		}
		return prog_alive;
	}

	std::vector<int> generate_migrant_list(int N_migrants)
	{
		std::vector<int> id_list;
		int i;
		double migrant_list_unfilled, progeny_left, keep_prob;
		int number_migrants_chosen;
		double drnd1;

		progeny_left = 0;
		for (i = 0; i < Nprogeny; i++)
			if (progeny[i].Alive)
				progeny_left++;

		migrant_list_unfilled = N_migrants;
		number_migrants_chosen = 0;
		for (i = 0; i < Nprogeny; i++)
		{
			if (progeny[i].Alive)
			{
				keep_prob = migrant_list_unfilled / progeny_left;
				drnd1 = genrand();
				if (drnd1 < keep_prob)
				{
					id_list.push_back(i);
					migrant_list_unfilled = migrant_list_unfilled - 1;
					number_migrants_chosen++;
				}
				progeny_left = progeny_left - 1;
			} // end of if (progeny[i].Alive)
		} // end of i
		return id_list;
	}

	void send_progeny_migrant(int id_number, individual &recipient)
	{
		deep_copy_individuals(recipient, progeny[id_number]);
	}

	void replace_progeny_with_migrant(int id_number, individual &migrant)
	{
		deep_copy_individuals(progeny[id_number], migrant);
	}

	void save_all_genotypes(std::string filename)
	{
		int i, j, k;

		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);

		outfile << "Individual,";
		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberMarkerLociPerChromosome; j++)
				outfile << "C" << i << "M" << j << ",";
		}

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberQTLsPerChromosome; j++)
				outfile << "C" << i << "Q0" << QTLlocus[i].QTLlocation[j] << ",";
			for (j = 0; j < NumberQTLsPerChromosome_1; j++)
				outfile << "C" << i << "Q1" << QTLlocus[i].QTLlocation_1[j] << ",";
			for (j = 0; j < NumberQTLsPerChromosome_P; j++)
				outfile << "C" << i << "QP" << QTLlocus[i].QTLlocation_P[j] << ",";
		}
		outfile << "Sex\n";

		for (k = 0; k < PopulationSize; k++)
		{
			outfile << "A" << k << ",";
			for (i = 0; i < NumberChromosomes; i++)
			{
				for (j = 0; j < NumberMarkerLociPerChromosome; j++)
				{
					outfile << adult[k].MaternalChromosome[i].MarkerLoci[j] << "x";
					outfile << adult[k].PaternalChromosome[i].MarkerLoci[j] << ",";
				}
			}
			for (i = 0; i < NumberChromosomes; i++)
			{
				for (j = 0; j < NumberQTLsPerChromosome; j++)
				{
					outfile << adult[k].MaternalChromosome[i].QTLeffect[j] << "x";
					outfile << adult[k].PaternalChromosome[i].QTLeffect[j] << ",";
				}
				for (j = 0; j < NumberQTLsPerChromosome_1; j++)
				{
					outfile << adult[k].MaternalChromosome[i].QTLeffect_1[j] << "x";
					outfile << adult[k].PaternalChromosome[i].QTLeffect_1[j] << ",";
				}
				for (j = 0; j < NumberQTLsPerChromosome_P; j++)
				{
					outfile << adult[k].MaternalChromosome[i].QTLeffect_P0[j] << "|";
					outfile << adult[k].MaternalChromosome[i].QTLeffect_P1[j] << "x";
					outfile << adult[k].PaternalChromosome[i].QTLeffect_P0[j] << "|";
					outfile << adult[k].PaternalChromosome[i].QTLeffect_P1[j] << ",";
				}
			}
			if (adult[k].Female)
				outfile << "Female";
			else
				outfile << "Male";
			outfile << "\n";
		}

		for (k = 0; k < Nprogeny; k++)
		{
			outfile << "P" << k << ",";
			for (i = 0; i < NumberChromosomes; i++)
			{
				for (j = 0; j < NumberMarkerLociPerChromosome; j++)
				{
					outfile << progeny[k].MaternalChromosome[i].MarkerLoci[j] << "x";
					outfile << progeny[k].PaternalChromosome[i].MarkerLoci[j] << ",";
				}
			}
			for (i = 0; i < NumberChromosomes; i++)
			{
				for (j = 0; j < NumberQTLsPerChromosome; j++)
				{
					outfile << progeny[k].MaternalChromosome[i].QTLeffect[j] << "x";
					outfile << progeny[k].PaternalChromosome[i].QTLeffect[j] << ",";
				}
				for (j = 0; j < NumberQTLsPerChromosome_1; j++)
				{
					outfile << progeny[k].MaternalChromosome[i].QTLeffect_1[j] << "x";
					outfile << progeny[k].PaternalChromosome[i].QTLeffect_1[j] << ",";
				}
				for (j = 0; j < NumberQTLsPerChromosome_P; j++)
				{
					outfile << progeny[k].MaternalChromosome[i].QTLeffect_P0[j] << "|";
					outfile << progeny[k].MaternalChromosome[i].QTLeffect_P1[j] << "x";
					outfile << progeny[k].PaternalChromosome[i].QTLeffect_P0[j] << "|";
					outfile << progeny[k].PaternalChromosome[i].QTLeffect_P1[j] << ",";
				}
			}
			if (progeny[k].Female)
				outfile << "Female";
			else
				outfile << "Male";
			outfile << "\n";
		}
	}

	std::vector<int> report_qtl_locations()
	{
		std::vector<int> locations;
		int i, j;
		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberQTLsPerChromosome; j++)
			{
				locations.push_back(QTLlocus[i].QTLlocation[j]);
			}
			for (j = 0; j < NumberQTLsPerChromosome_1; j++)
			{
				locations.push_back(QTLlocus[i].QTLlocation_1[j]);
			}
			for (j = 0; j < NumberQTLsPerChromosome_P; j++)
			{
				locations.push_back(QTLlocus[i].QTLlocation_P[j]);
			}
		}
		return locations;
	}

	void set_qtl_locations(std::vector<int> locations)
	{
		int counter = 0;
		int i, j;
		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberQTLsPerChromosome; j++)
			{
				QTLlocus[i].QTLlocation[j] = locations[counter];
				counter++;
			}
			for (j = 0; j < NumberQTLsPerChromosome_1; j++)
			{
				QTLlocus[i].QTLlocation_1[j] = locations[counter];
				counter++;
			}
			for (j = 0; j < NumberQTLsPerChromosome_P; j++)
			{
				QTLlocus[i].QTLlocation_P[j] = locations[counter];
				counter++;
			}
		}
	}

	void sample_adults()
	{
		int iii;
		double dNumberAdultsLeft, dSampleUnfilled;
		double dSampleProb;
		int iNumSampledSoFar;

		dNumberAdultsLeft = PopulationSize;
		dSampleUnfilled = SampleSizeAdults;
		iNumSampledSoFar = 0;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			dSampleProb = dSampleUnfilled / dNumberAdultsLeft;
			if (genrand() < dSampleProb)
			{
				adult_sample[iNumSampledSoFar] = iii;
				iNumSampledSoFar++;
				dSampleUnfilled = dSampleUnfilled - 1;
			}
			dNumberAdultsLeft = dNumberAdultsLeft - 1;
		}
		ActualAdultSampleSize = iNumSampledSoFar;
	}

	marker_allele_frequencies calculate_marker_freqs_adult_sample(int whichchromosome, int whichmarker)
	{
		marker_allele_frequencies temp_maf;
		int iii, tempindex;
		double *tempallelefreq = new double[MaxNumAllelesPerMarker];

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = 0;

		double allelecounter = 0;
		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempindex = adult_sample[iii];
			tempallelefreq[adult[tempindex].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
			tempallelefreq[adult[tempindex].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
			allelecounter = allelecounter + 2;
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			tempallelefreq[iii] = tempallelefreq[iii] / allelecounter;
			temp_maf.freq.push_back(tempallelefreq[iii]);
		}

		delete[] tempallelefreq;
		return temp_maf;
	}

	qtl_allele_frequencies calculate_qtl_allele_frequencies_trt0_adult_sample(int whichchromosome, int whichqtl)
	{
		qtl_allele_frequencies temp_qtl_af;
		std::vector<double> AlleleList;
		std::vector<double> AlleleFreqAdults;
		int iii, jjj;
		int NumAllelesInList = 0;
		bool Onlist;
		int tempid;

		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempid = adult_sample[iii];
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].MaternalChromosome[whichchromosome].QTLeffect[whichqtl] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(adult[tempid].MaternalChromosome[whichchromosome].QTLeffect[whichqtl]);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].PaternalChromosome[whichchromosome].QTLeffect[whichqtl] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(adult[tempid].PaternalChromosome[whichchromosome].QTLeffect[whichqtl]);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
		}

		// Now we have a list of all of the alleles for this qtl
		// Figure out allele frequencies

		double allelecounter = 0;
		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempid = adult_sample[iii];
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].MaternalChromosome[whichchromosome].QTLeffect[whichqtl] == AlleleList[jjj])
					AlleleFreqAdults[jjj]++;
				if (adult[tempid].PaternalChromosome[whichchromosome].QTLeffect[whichqtl] == AlleleList[jjj])
					AlleleFreqAdults[jjj]++;
			}
			allelecounter = allelecounter + 2;
		}

		for (jjj = 0; jjj < NumAllelesInList; jjj++)
		{
			AlleleFreqAdults[jjj] = AlleleFreqAdults[jjj] / allelecounter;
		}

		// put the data in the allele freq class
		for (jjj = 0; jjj < NumAllelesInList; jjj++)
		{
			temp_qtl_af.trt0_effect.push_back(AlleleList[jjj]);
			temp_qtl_af.trt1_effect.push_back(0);
			temp_qtl_af.freq.push_back(AlleleFreqAdults[jjj]);
		}

		return temp_qtl_af;
	}

	qtl_allele_frequencies calculate_qtl_allele_frequencies_pleio_adult_sample(int whichchromosome, int whichqtl)
	{
		qtl_allele_frequencies temp_qtl_af;
		std::vector<double> AlleleList0;
		std::vector<double> AlleleList1;
		std::vector<double> AlleleFreqAdults;
		int iii, jjj;
		int NumAllelesInList = 0;
		bool Onlist;
		int tempid;

		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempid = adult_sample[iii];
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_P0[whichqtl] == AlleleList0[jjj]
					&& adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_P1[whichqtl] == AlleleList1[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList0.push_back(adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_P0[whichqtl]);
				AlleleList1.push_back(adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_P1[whichqtl]);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_P0[whichqtl] == AlleleList0[jjj]
					&& adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_P1[whichqtl] == AlleleList1[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList0.push_back(adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_P0[whichqtl]);
				AlleleList1.push_back(adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_P1[whichqtl]);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
		}

		// Now we have a list of all of the alleles for this qtl
		// Figure out allele frequencies

		double allelecounter = 0;
		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempid = adult_sample[iii];
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_P0[whichqtl] == AlleleList0[jjj]
					&& adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_P1[whichqtl] == AlleleList1[jjj])
					AlleleFreqAdults[jjj]++;
				if (adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_P0[whichqtl] == AlleleList0[jjj]
					&& adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_P1[whichqtl] == AlleleList1[jjj])
					AlleleFreqAdults[jjj]++;
			}
			allelecounter = allelecounter + 2;
		}

		for (jjj = 0; jjj < NumAllelesInList; jjj++)
		{
			AlleleFreqAdults[jjj] = AlleleFreqAdults[jjj] / allelecounter;
		}

		// put the data in the allele freq class
		for (jjj = 0; jjj < NumAllelesInList; jjj++)
		{
			temp_qtl_af.trt0_effect.push_back(AlleleList0[jjj]);
			temp_qtl_af.trt1_effect.push_back(AlleleList1[jjj]);
			temp_qtl_af.freq.push_back(AlleleFreqAdults[jjj]);
		}

		return temp_qtl_af;
	}

	qtl_allele_frequencies calculate_qtl_allele_frequencies_trt1_adult_sample(int whichchromosome, int whichqtl)
	{
		qtl_allele_frequencies temp_qtl_af;
		std::vector<double> AlleleList;
		std::vector<double> AlleleFreqAdults;
		int iii, jjj;
		int NumAllelesInList = 0;
		bool Onlist;
		int tempid;

		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempid = adult_sample[iii];
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_1[whichqtl] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_1[whichqtl]);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_1[whichqtl] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_1[whichqtl]);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
		}

		// Now we have a list of all of the alleles for this qtl
		// Figure out allele frequencies

		double allelecounter = 0;
		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempid = adult_sample[iii];
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (adult[tempid].MaternalChromosome[whichchromosome].QTLeffect_1[whichqtl] == AlleleList[jjj])
					AlleleFreqAdults[jjj]++;
				if (adult[tempid].PaternalChromosome[whichchromosome].QTLeffect_1[whichqtl] == AlleleList[jjj])
					AlleleFreqAdults[jjj]++;
			}
			allelecounter = allelecounter + 2;
		}

		for (jjj = 0; jjj < NumAllelesInList; jjj++)
		{
			AlleleFreqAdults[jjj] = AlleleFreqAdults[jjj] / allelecounter;
		}

		// put the data in the allele freq class
		for (jjj = 0; jjj < NumAllelesInList; jjj++)
		{
			temp_qtl_af.trt0_effect.push_back(0);
			temp_qtl_af.trt1_effect.push_back(AlleleList[jjj]);
			temp_qtl_af.freq.push_back(AlleleFreqAdults[jjj]);
		}

		return temp_qtl_af;
	}

	void calculate_allele_frequencies_all_chromosomes()
	{
		int i, j;

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberMarkerLociPerChromosome; j++)
			{
				chr_data[i].marker_afs[j] = calculate_marker_freqs_adult_sample(i, j);
			}

			for (j = 0; j < NumberQTLsPerChromosome; j++)
			{
				chr_data[i].qtl_afs_trt0[j] = calculate_qtl_allele_frequencies_trt0_adult_sample(i, j);
				chr_data[i].qtl_loc_trt0[j] = QTLlocus[i].QTLlocation[j];
			}

			for (j = 0; j < NumberQTLsPerChromosome_1; j++)
			{
				chr_data[i].qtl_afs_trt1[j] = calculate_qtl_allele_frequencies_trt1_adult_sample(i, j);
				chr_data[i].qtl_loc_trt1[j] = QTLlocus[i].QTLlocation_1[j];
			}

			for (j = 0; j < NumberQTLsPerChromosome_P; j++)
			{
				chr_data[i].qtl_afs_pleio[j] = calculate_qtl_allele_frequencies_pleio_adult_sample(i, j);
				chr_data[i].qtl_loc_pleio[j] = QTLlocus[i].QTLlocation_P[j];
			}
		}
	}

	void save_marker_freqs_to_file(std::string filename)
	{
		int i, j, k;

		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);

		outfile << "Locus";
		for (j = 0; j < MaxNumAllelesPerMarker; j++)
		{
			outfile << ",A" << j;
		}
		outfile << "\n";

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberMarkerLociPerChromosome; j++)
			{
				outfile << "C" << i << "M" << j;
				for (k = 0; k < MaxNumAllelesPerMarker; k++)
				{
					outfile << "," << chr_data[i].marker_afs[j].freq[k];
				}
				outfile << "\n";
			}
		}

		outfile.close();

	}

	void save_qtl_freqs_to_file(std::string filename)
	{
		int i, j, k;

		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberQTLsPerChromosome; j++)
			{
				outfile << "C" << i << "Q" << j;
				outfile << "_Location_" << QTLlocus[i].QTLlocation[j];
				outfile << "\ntrt0effect,trt1effect,freq,";
				outfile << "trt0mean,trt0var,trt1mean,trt1var";
				for (k = 0; k < chr_data[i].qtl_afs_trt0[j].freq.size(); k++)
				{
					outfile << "\n" << chr_data[i].qtl_afs_trt0[j].trt0_effect[k];
					outfile << "," << chr_data[i].qtl_afs_trt0[j].trt1_effect[k];
					outfile << "," << chr_data[i].qtl_afs_trt0[j].freq[k];
					if (k == 0)
					{
						outfile << "," << chr_data[i].qtl_mean_trt0[j];
						outfile << "," << chr_data[i].qtl_var_trt0[j];
						outfile << ",0,0";
					}
				}
				outfile << "\n";
			}
		}

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberQTLsPerChromosome_1; j++)
			{
				outfile << "C" << i << "Q" << j;
				outfile << "_Location_" << QTLlocus[i].QTLlocation_1[j];
				outfile << "\ntrt0effect,trt1effect,freq,";
				outfile << "trt0mean,trt0var,trt1mean,trt1var";
				for (k = 0; k < chr_data[i].qtl_afs_trt1[j].freq.size(); k++)
				{
					outfile << "\n" << chr_data[i].qtl_afs_trt1[j].trt0_effect[k];
					outfile << "," << chr_data[i].qtl_afs_trt1[j].trt1_effect[k];
					outfile << "," << chr_data[i].qtl_afs_trt1[j].freq[k];
					if (k == 0)
					{
						outfile << ",0,0";
						outfile << "," << chr_data[i].qtl_mean_trt1[j];
						outfile << "," << chr_data[i].qtl_var_trt1[j];
					}
				}
				outfile << "\n";
			}
		}

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberQTLsPerChromosome_P; j++)
			{
				outfile << "C" << i << "Q" << j;
				outfile << "_Location_" << QTLlocus[i].QTLlocation_P[j];
				outfile << "\ntrt0effect,trt1effect,freq,";
				outfile << "trt0mean,trt0var,trt1mean,trt1var";
				for (k = 0; k < chr_data[i].qtl_afs_pleio[j].freq.size(); k++)
				{
					outfile << "\n" << chr_data[i].qtl_afs_pleio[j].trt0_effect[k];
					outfile << "," << chr_data[i].qtl_afs_pleio[j].trt1_effect[k];
					outfile << "," << chr_data[i].qtl_afs_pleio[j].freq[k];
					if (k == 0)
					{
						outfile << "," << chr_data[i].qtl_mean_pleio_0[j];
						outfile << "," << chr_data[i].qtl_var_pleio_0[j];
						outfile << "," << chr_data[i].qtl_mean_pleio_1[j];
						outfile << "," << chr_data[i].qtl_var_pleio_1[j];
					}
				}
				outfile << "\n";
			}
		}

		outfile.close();

	}

	int report_number_chromosomes()
	{
		return NumberChromosomes;
	}

	int report_markers_per_chromosome()
	{
		return NumberMarkerLociPerChromosome;
	}

	int report_qtl0_per_chromosome()
	{
		return NumberQTLsPerChromosome;
	}

	int report_qtl1_per_chromosome()
	{
		return NumberQTLsPerChromosome_1;
	}

	int report_qtlP_per_chromosome()
	{
		return NumberQTLsPerChromosome_P;
	}

	chromosomedata report_chromosome_data(int whichchromosome)
	{
		return chr_data[whichchromosome];
	}

	int report_all_calc_int()
	{
		return all_calc_interval;
	}

	double report_fst_weight_var()
	{
		return fst_weighting_variance;
	}

	int report_sample_size()
	{
		return ActualAdultSampleSize;
	}

	void calc_per_qtl_means_and_vars_progeny(int whichchromosome)
	{
		int i, j;
		double temp;
		double tmean;
		double nprog;
		double tvar;

		nprog = static_cast<double>(Nprogeny);

		// QTLs for Trait 0

		chr_data[whichchromosome].qtl_mean_trt0.clear();
		chr_data[whichchromosome].qtl_var_trt0.clear();

		for (i = 0; i < NumberQTLsPerChromosome; i++)
		{
			tmean = 0;
			tvar = 0;
			for (j = 0; j < Nprogeny; j++)
			{

				temp = progeny[j].MaternalChromosome[whichchromosome].QTLeffect[i] + progeny[j].PaternalChromosome[whichchromosome].QTLeffect[i];
				tmean = tmean + temp;
				tvar = tvar + temp * temp;
			}
			if (nprog > 0)
			{
				tmean = tmean / nprog;
				tvar = tvar / nprog - tmean * tmean;
			}
			chr_data[whichchromosome].qtl_mean_trt0.push_back(tmean);
			chr_data[whichchromosome].qtl_var_trt0.push_back(tvar);
		}


		// QTLs for Trait 1

		chr_data[whichchromosome].qtl_mean_trt1.clear();
		chr_data[whichchromosome].qtl_var_trt1.clear();

		for (i = 0; i < NumberQTLsPerChromosome_1; i++)
		{
			tmean = 0;
			tvar = 0;
			for (j = 0; j < Nprogeny; j++)
			{

				temp = progeny[j].MaternalChromosome[whichchromosome].QTLeffect_1[i] + progeny[j].PaternalChromosome[whichchromosome].QTLeffect_1[i];
				tmean = tmean + temp;
				tvar = tvar + temp * temp;
			}
			if (nprog > 0)
			{
				tmean = tmean / nprog;
				tvar = tvar / nprog - tmean * tmean;
			}
			chr_data[whichchromosome].qtl_mean_trt1.push_back(tmean);
			chr_data[whichchromosome].qtl_var_trt1.push_back(tvar);
		}

		// QTLs for Pleiotropic Loci

		chr_data[whichchromosome].qtl_mean_pleio_0.clear();
		chr_data[whichchromosome].qtl_var_pleio_0.clear();
		chr_data[whichchromosome].qtl_mean_pleio_1.clear();
		chr_data[whichchromosome].qtl_var_pleio_1.clear();

		for (i = 0; i < NumberQTLsPerChromosome_P; i++)
		{
			tmean = 0;
			tvar = 0;
			for (j = 0; j < Nprogeny; j++)
			{

				temp = progeny[j].MaternalChromosome[whichchromosome].QTLeffect_P0[i] + progeny[j].PaternalChromosome[whichchromosome].QTLeffect_P0[i];
				tmean = tmean + temp;
				tvar = tvar + temp * temp;
			}
			if (nprog > 0)
			{
				tmean = tmean / nprog;
				tvar = tvar / nprog - tmean * tmean;
			}
			chr_data[whichchromosome].qtl_mean_pleio_0.push_back(tmean);
			chr_data[whichchromosome].qtl_var_pleio_0.push_back(tvar);

			tmean = 0;
			tvar = 0;
			for (j = 0; j < Nprogeny; j++)
			{

				temp = progeny[j].MaternalChromosome[whichchromosome].QTLeffect_P1[i] + progeny[j].PaternalChromosome[whichchromosome].QTLeffect_P1[i];
				tmean = tmean + temp;
				tvar = tvar + temp * temp;
			}
			if (nprog > 0)
			{
				tmean = tmean / nprog;
				tvar = tvar / nprog - tmean * tmean;
			}
			chr_data[whichchromosome].qtl_mean_pleio_1.push_back(tmean);
			chr_data[whichchromosome].qtl_var_pleio_1.push_back(tvar);
		}

	}

	void calc_qtl_means_vars_all_chromosomes()
	{
		int i;

		for (i = 0; i < NumberChromosomes; i++)
			calc_per_qtl_means_and_vars_progeny(i);
	}

	void estimate_additive_variances_progeny()
	{
		// Basically just choose two offspring at random
		// and have them produce hypothetical progeny.
		// Use parent-offspring resemblance to estimate the
		// additive genetic variance.
		int parent_1, parent_2;
		int i, m, n;
		int n_families;
		int temp_fecund;
		double d_tf, d_nf;

		n_families = 500; // Create 500 families
		temp_fecund = 4; // must be 10 or less (memory is allocated for only 10 temp progeny)
		d_tf = temp_fecund;
		d_nf = n_families;

		double *midparent00 = new double[n_families];
		double *midparent11 = new double[n_families];
		double *midoffspring00 = new double[n_families];
		double *midoffspring11 = new double[n_families];

		for (i = 0; i < n_families; i++) 
		{
			parent_1 = randnum(Nprogeny);
			parent_2 = randnum(Nprogeny);
			
			// Produce hypothetical progeny
			
			for (m = 0; m < temp_fecund; m++)
			{
				for (n = 0; n < NumberChromosomes; n++)
				{
					ProduceRecombinedChromosome(temp_progeny[m].MaternalChromosome[n], progeny[parent_1], n, ExpRecombPerChromosome);
					ProduceRecombinedChromosome(temp_progeny[m].PaternalChromosome[n], progeny[parent_2], n, ExpRecombPerChromosome);
				}

				temp_progeny[m].calculate_genotypic_values(NumberChromosomes, NumberQTLsPerChromosome, NumberQTLsPerChromosome_1, NumberQTLsPerChromosome_P, epistasis_allowed, ep_par_collection);
				temp_progeny[m].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
			} // end of m loop
			
			midparent00[i] = (progeny[parent_1].Phenotype[0] + progeny[parent_2].Phenotype[0])/2;
			midparent11[i] = (progeny[parent_1].Phenotype[1] + progeny[parent_2].Phenotype[1])/2;
			midoffspring00[i] = 0;
			midoffspring11[i] = 0;
			for (m = 0; m < temp_fecund; m++)
			{
				midoffspring00[i] = midoffspring00[i] + temp_progeny[m].Phenotype[0];
				midoffspring11[i] = midoffspring11[i] + temp_progeny[m].Phenotype[1];
			}
			midoffspring00[i] = midoffspring00[i] / d_tf;
			midoffspring11[i] = midoffspring11[i] / d_tf;

		} // end of i loop

		// Calculate the covariances as estimates of additive genetic variances and covariances
		double cov00, cov11, cov01, cov10;
		double mean_parent0, mean_parent1, mean_off0, mean_off1;

		mean_parent0 = 0;
		mean_parent1 = 0;
		mean_off0 = 0;
		mean_off1 = 0;
		for (i = 0; i < n_families; i++)
		{
			mean_parent0 = mean_parent0 + midparent00[i];
			mean_parent1 = mean_parent1 + midparent11[i];
			mean_off0 = mean_off0 + midoffspring00[i];
			mean_off1 = mean_off1 + midoffspring11[i];
		}
		mean_parent0 = mean_parent0 / d_nf;
		mean_parent1 = mean_parent1 / d_nf;
		mean_off0 = mean_off0 / d_nf;
		mean_off1 = mean_off1 / d_nf;

		cov00 = 0;
		cov11 = 0;
		cov01 = 0;
		cov10 = 0;
		for (i = 0; i < n_families; i++)
		{
			cov00 = cov00 + (midparent00[i] - mean_parent0)*(midoffspring00[i] - mean_off0);
			cov11 = cov11 + (midparent11[i] - mean_parent1)*(midoffspring11[i] - mean_off1);
			cov01 = cov01 + (midparent00[i] - mean_parent0)*(midoffspring11[i] - mean_off1);
			cov10 = cov10 + (midparent11[i] - mean_parent1)*(midoffspring00[i] - mean_off0);
		}
		cov00 = cov00 / (d_nf - 1);
		cov11 = cov11 / (d_nf - 1);
		cov01 = cov01 / (d_nf - 1);
		cov10 = cov10 / (d_nf - 1);

		est_va00 = cov00 * 2;
		est_va11 = cov11 * 2;
		est_va01 = (cov01 + cov10);

		delete[] midparent00;
		delete[] midparent11;
		delete[] midoffspring00;
		delete[] midoffspring11;

	}

	void report_means(std::vector<double> &variables)
	{
		variables.push_back(sim_data.meanZbar0);
		variables.push_back(sim_data.meanZbar1);
		variables.push_back(sim_data.meanG00);
		variables.push_back(sim_data.meanG11);
		variables.push_back(sim_data.meanG01);
		variables.push_back(sim_data.meanV0);
		variables.push_back(sim_data.meanV1);
		variables.push_back(sim_data.meanV01);
	}

	void report_parameter_vals(std::vector<double> &parameters)
	{
		parameters.push_back(NumberOfGenerations);
		parameters.push_back(PopulationSize);
		parameters.push_back(SampleSizeAdults);
		parameters.push_back(NumberChromosomes);
		parameters.push_back(NumberMarkerLociPerChromosome);
		parameters.push_back(MutationRatePerMarker);
		parameters.push_back(NumberQTLsPerChromosome);
		parameters.push_back(NumberQTLsPerChromosome_1);
		parameters.push_back(NumberQTLsPerChromosome_P);
		parameters.push_back(MutationRatePerLocus);
		parameters.push_back(EnvironmentalVariance[0]);
		parameters.push_back(EnvironmentalVariance[1]);
		parameters.push_back(ExpRecombPerChromosome);
		parameters.push_back(SelectionStrength[0]);
		parameters.push_back(SelectionStrength[1]);
		parameters.push_back(SelectionalCorrelation);
		parameters.push_back(Optimum[0]);
		parameters.push_back(Optimum[1]);
	}

	void report_optima(std::vector<double> &parameters)
	{
		parameters.push_back(Optimum[0]);
		parameters.push_back(Optimum[1]);
	}

	double CalculateLDprogeny(int whichchromosome, int markerA, int markerB)
	{
		int iii, jjj;
		double *allelefreqA = new double[MaxNumAllelesPerMarker];
		double *allelefreqB = new double[MaxNumAllelesPerMarker];
		double **jointAB;
		jointAB = new double *[MaxNumAllelesPerMarker];
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			jointAB[iii] = new double[MaxNumAllelesPerMarker];
		double count;
		double **Dij;
		Dij = new double *[MaxNumAllelesPerMarker];
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			Dij[iii] = new double[MaxNumAllelesPerMarker];
		double **Dmax;
		Dmax = new double *[MaxNumAllelesPerMarker];
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			Dmax[iii] = new double[MaxNumAllelesPerMarker];

		double Dprime;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			allelefreqA[iii] = 0;
			allelefreqB[iii] = 0;
			for (jjj = 0; jjj < MaxNumAllelesPerMarker; jjj++)
			{
				jointAB[iii][jjj] = 0;
			}
		}

		count = 0;
		int maternalalleleA, maternalalleleB, paternalalleleA, paternalalleleB;
		for (iii = 0; iii < Nprogeny; iii++)
		{
			maternalalleleA = progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[markerA];
			maternalalleleB = progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[markerB];
			paternalalleleA = progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[markerA];
			paternalalleleB = progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[markerB];

			allelefreqA[maternalalleleA]++;
			allelefreqA[paternalalleleA]++;
			allelefreqB[maternalalleleB]++;
			allelefreqB[paternalalleleB]++;
			jointAB[maternalalleleA][maternalalleleB]++;
			jointAB[paternalalleleA][paternalalleleB]++;

			count = count + 2;
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			allelefreqA[iii] = allelefreqA[iii] / count;
			allelefreqB[iii] = allelefreqB[iii] / count;
			for (jjj = 0; jjj < MaxNumAllelesPerMarker; jjj++)
				jointAB[iii][jjj] = jointAB[iii][jjj] / count;
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			for (jjj = 0; jjj < MaxNumAllelesPerMarker; jjj++)
			{
				if (allelefreqA[iii] > 0 && allelefreqB[jjj] > 0)
				{
					Dij[iii][jjj] = jointAB[iii][jjj] - allelefreqA[iii] * allelefreqB[jjj];
					if (Dij[iii][jjj] < 0)
						Dmax[iii][jjj] = std::min(allelefreqA[iii] * allelefreqB[jjj], (1 - allelefreqA[iii])*(1 - allelefreqB[jjj]));
					else
						Dmax[iii][jjj] = std::min((1 - allelefreqA[iii])*allelefreqB[jjj], allelefreqA[iii] * (1 - allelefreqB[jjj]));

				}
			}
		}

		Dprime = 0;
		bool decentDmax = true;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			for (jjj = 0; jjj < MaxNumAllelesPerMarker; jjj++)
			{
				if (allelefreqA[iii] > 0 && allelefreqB[jjj] > 0)
				{
					if (Dmax[iii][jjj] > 0)
						Dprime = Dprime + allelefreqA[iii] * allelefreqB[jjj] * fabs(Dij[iii][jjj]) / Dmax[iii][jjj];
					else
						decentDmax = false;

				}
			}
		}
		if (!decentDmax)
			Dprime = -5;

		// output for checking

		//cout << "\n";
		//for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		//{
		//	cout << "Allele " << iii << ":\t" << allelefreqA[iii] << "\t" << allelefreqB[iii] << "\n";
		//}
		//
		//for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		//{
		//	for (jjj = 0; jjj < MaxNumAllelesPerMarker; jjj++)
		//	{
		//		cout << iii << ", " << jjj << ":\t" << jointAB[iii][jjj] << "\n";
		//	}
		//}


		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			delete[] jointAB[iii];
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			delete[] Dij[iii];
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			delete[] Dmax[iii];
		delete[] jointAB;
		delete[] Dij;
		delete[] Dmax;
		delete[] allelefreqA;
		delete[] allelefreqB;

		return Dprime;
	}

	void save_pairwise_LD(std::string filename)
	{
		int i, j;

		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);

		outfile << "LocusA,LocusB,Genome_Loc,LD\n";

		for (i = 0; i < NumberChromosomes; i++)
		{
			for (j = 0; j < NumberMarkerLociPerChromosome - 1; j++)
			{
				outfile << j << "," << j + 1;
				outfile << "," <<  i * NumberMarkerLociPerChromosome + j;
				outfile << "," << CalculateLDprogeny(i, j, j + 1);
				outfile << "\n";
			}
		}

		outfile.close();
	}

};

class meta_population
{
	// The meta-population is a collection of populations.
	// Each population has its own simulation engine.
	// Populations are completely independent, except for migration.
	// And their chromosomes have to be the same -- same numbers
	// of markers, same numbers of QTLs, same QTL locations.

private:
	int Number_of_populations;
	double migration_rate;
	simulation_engine *population;
	int r_num_seed;
	std::vector<chromosomedata> chrom_data_pop0;
	std::vector<chromosomedata> chrom_data_pop1;
	std::vector<chromosomedata> chrom_data_combined;

public:
	meta_population()
	{

	}

	~meta_population()
	{

	}

	void initialize_metapop(double migr_rate)
	{
		// Seed the random number generator
		r_num_seed = static_cast<int>(time(NULL));
		sgenrand(r_num_seed);

		Number_of_populations = 2;
		migration_rate = migr_rate;
		population = new simulation_engine[Number_of_populations];

	}

	void deinitialize()
	{
		delete[] population;
	}

	void set_parameter_values(parameter_value_set parms)
	{
		// Establish meta-population-wide epistatic parameter
		// values.

		epistatic_parameter_collection metapop_epistatic_parms;
		metapop_epistatic_parms.set_epistatic_effects(parms.p_number_chromosomes, parms.p_loci_trt0, parms.p_loci_trt1, parms.p_loci_pleio, parms.p_epistatic_vars);

		int m;
		for (m = 0; m < Number_of_populations; m++)
		{
			population[m].update_parameter_values(parms, metapop_epistatic_parms);
		}

		// make any population-specific parameter value changes here

		// For now just change population two to be different
		parameter_value_set parms_2 = parms;
		parms_2.p_exp_opt0 = parms.pop2_opt0;
		parms_2.p_exp_opt1 = parms.pop2_opt1;
		parms_2.p_exp_w00 = parms.pop2_w00;
		parms_2.p_exp_w11 = parms.pop2_w11;
		parms_2.p_exp_sel_corr = parms.pop2_sel_corr;
		population[1].update_parameter_values(parms_2, metapop_epistatic_parms);

	}

	void save_parameter_values(std::string outfile_name)
	{
		std::stringstream sss;
		std::string big_fn;
		int m;
		for (m = 0; m < Number_of_populations; m++)
		{
			sss.str("");
			sss << outfile_name << "_parameters_pop_" << m << ".csv";
			big_fn = sss.str();
			population[m].save_parameter_values(big_fn, migration_rate);
		}
	}

	void progeny_stepping_stone_migration(simulation_engine &pop_a, simulation_engine &pop_b)
	{
		// Figure out the number of migrants
		int i, iNmigrants;
		double pop_a_N, pop_b_N, expected_N_migrants;

		pop_a_N = static_cast<double>(pop_a.getNumberofProgenyAlive());
		pop_b_N = static_cast<double>(pop_b.getNumberofProgenyAlive());
		expected_N_migrants = migration_rate * (pop_a_N + pop_b_N) / 2;
		if (expected_N_migrants > pop_a_N)
			expected_N_migrants = pop_a_N;
		if (expected_N_migrants > pop_b_N)
			expected_N_migrants = pop_b_N;
		
		if (genrand() < expected_N_migrants - floor(expected_N_migrants))
		{
			iNmigrants = static_cast<int>(floor(expected_N_migrants)) + 1;
		}
		else
		{
			iNmigrants = static_cast<int>(floor(expected_N_migrants));
		}

		std::vector<int> migrant_id_list_a, migrant_id_list_b;
		migrant_id_list_a = pop_a.generate_migrant_list(iNmigrants);
		migrant_id_list_b = pop_b.generate_migrant_list(iNmigrants);

		if (static_cast<int>(migrant_id_list_a.size()) < iNmigrants)
			iNmigrants = static_cast<int>(migrant_id_list_a.size());
		if (static_cast<int>(migrant_id_list_b.size()) < iNmigrants)
			iNmigrants = static_cast<int>(migrant_id_list_b.size());

		// Now swap the migrants
		for (i = 0; i < iNmigrants; i++)
		{
			pop_a.send_progeny_migrant(migrant_id_list_a[i], pop_a.migrant);
			pop_b.send_progeny_migrant(migrant_id_list_b[i], pop_b.migrant);
			pop_a.replace_progeny_with_migrant(migrant_id_list_a[i], pop_b.migrant);
			pop_b.replace_progeny_with_migrant(migrant_id_list_b[i], pop_a.migrant);
		}

	}

	qtl_allele_frequencies combine_two_qtl_allele_frequency_lists(qtl_allele_frequencies &qtl0, qtl_allele_frequencies &qtl1)
	{
		qtl_allele_frequencies qtl_both;
		int i, j;
		// First make a list of all allele names and add it to qtl_both
		
		bool on_list;
		
		for (i = 0; i < qtl0.freq.size(); i++)
		{
			on_list = false;
			for (j = 0; j < qtl_both.freq.size(); j++)
			{
				if (qtl0.trt0_effect[i] == qtl_both.trt0_effect[j] &&
					qtl0.trt1_effect[i] == qtl_both.trt1_effect[j])
					on_list = true;
			}
			
			if (!on_list)
			{
				qtl_both.trt0_effect.push_back(qtl0.trt0_effect[i]);
				qtl_both.trt1_effect.push_back(qtl0.trt1_effect[i]);
				qtl_both.freq.push_back(0);
			}
		}

		for (i = 0; i < qtl1.freq.size(); i++)
		{
			on_list = false;
			for (j = 0; j < qtl_both.freq.size(); j++)
			{
				if (qtl1.trt0_effect[i] == qtl_both.trt0_effect[j] &&
					qtl1.trt1_effect[i] == qtl_both.trt1_effect[j])
					on_list = true;
			}

			if (!on_list)
			{
				qtl_both.trt0_effect.push_back(qtl1.trt0_effect[i]);
				qtl_both.trt1_effect.push_back(qtl1.trt1_effect[i]);
				qtl_both.freq.push_back(0);
			}
		}

		// Now, given the list, we need to calculate the average allele frequency across pops

		double freq0, freq1;
		for (i = 0; i < qtl_both.freq.size(); i++)
		{
			freq0 = 0;
			for (j = 0; j < qtl0.freq.size(); j++)
			{
				if (qtl0.trt0_effect[j] == qtl_both.trt0_effect[i] &&
					qtl0.trt1_effect[j] == qtl_both.trt1_effect[i])
					freq0 = qtl0.freq[j];
			}
			freq1 = 0;
			for (j = 0; j < qtl1.freq.size(); j++)
			{
				if (qtl1.trt0_effect[j] == qtl_both.trt0_effect[i] &&
					qtl1.trt1_effect[j] == qtl_both.trt1_effect[i])
					freq1 = qtl1.freq[j];
			}

			qtl_both.freq[i] = (freq0 + freq1) / 2.0;
		}

		return qtl_both;

	}

	chromosomedata combine_chromosome_data(chromosomedata &p0, chromosomedata &p1)
	{
		chromosomedata combo;
		marker_allele_frequencies temp_maf;

		// Calculate average allele frequencies across all marker loci

		int i, j;
		double freq0, freq1, freqcomb;

		for (i = 0; i < p0.marker_afs.size(); i++)
		{
			temp_maf.freq.clear();
			for (j = 0; j < p0.marker_afs[i].freq.size(); j++)
			{
				freq0 = p0.marker_afs[i].freq[j];
				freq1 = p1.marker_afs[i].freq[j];
				freqcomb = (freq0 + freq1) / 2.0;
				temp_maf.freq.push_back(freqcomb);
			}
			combo.marker_afs.push_back(temp_maf);
		}

		// Calculate average allele frequencies across all qtls
		qtl_allele_frequencies temp_qtl_af;
		for (i = 0; i < p0.qtl_afs_trt0.size(); i++)
		{
			temp_qtl_af = combine_two_qtl_allele_frequency_lists(p0.qtl_afs_trt0[i], p1.qtl_afs_trt0[i]);
			combo.qtl_afs_trt0.push_back(temp_qtl_af);
		}
		for (i = 0; i < p0.qtl_afs_trt1.size(); i++)
		{
			temp_qtl_af = combine_two_qtl_allele_frequency_lists(p0.qtl_afs_trt1[i], p1.qtl_afs_trt1[i]);
			combo.qtl_afs_trt1.push_back(temp_qtl_af);
		}
		for (i = 0; i < p0.qtl_afs_pleio.size(); i++)
		{
			temp_qtl_af = combine_two_qtl_allele_frequency_lists(p0.qtl_afs_pleio[i], p1.qtl_afs_pleio[i]);
			combo.qtl_afs_pleio.push_back(temp_qtl_af);
		}

		combo.qtl_loc_trt0 = p0.qtl_loc_trt0;
		combo.qtl_loc_trt1 = p0.qtl_loc_trt1;
		combo.qtl_loc_pleio = p0.qtl_loc_pleio;

		return combo;
	}

	void save_combined_marker_data(chromosomedata &p0, chromosomedata &p1, chromosomedata &combo, std::string filename)
	{
		int i, j;
		
		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);

		outfile << "Locus,Allele,Pop0,Pop1,Both";
		for (i = 0; i < combo.marker_afs.size(); i++)
		{
			for (j = 0; j < combo.marker_afs[i].freq.size(); j++)
			{
				outfile << "\nLoc_" << i << ",";
				outfile << "All_" << j << ",";
				outfile << p0.marker_afs[i].freq[j] << ",";
				outfile << p1.marker_afs[i].freq[j] << ",";
				outfile << combo.marker_afs[i].freq[j];
			}
		}

		outfile.close();
		
	}

	void save_combined_qtl_data(chromosomedata &p0, chromosomedata &p1, chromosomedata &combo, std::string filename)
	{
		int i, j, k;
		bool onlist;

		std::ofstream outfile;
		char temp_fn[256];
		convert_string_to_char(filename, temp_fn);
		outfile.open(temp_fn);
		
		outfile << "Locus,Location,Allele,Trt0_eff,Trt1_eff,Pop0,Pop1,Both";
		for (i = 0; i < combo.qtl_afs_trt0.size(); i++)
		{
			for (j = 0; j < combo.qtl_afs_trt0[i].freq.size(); j++)
			{
				outfile << "\nQ_" << i << "," << combo.qtl_loc_trt0[i] << ",";
				outfile << "All_" << j << ",";
				outfile << combo.qtl_afs_trt0[i].trt0_effect[j] << ",";
				outfile << combo.qtl_afs_trt0[i].trt1_effect[j] << ",";

				onlist = false;
				for (k = 0; k < p0.qtl_afs_trt0[i].freq.size(); k++)
				{
					if (p0.qtl_afs_trt0[i].trt0_effect[k] == combo.qtl_afs_trt0[i].trt0_effect[j] &&
						p0.qtl_afs_trt0[i].trt1_effect[k] == combo.qtl_afs_trt0[i].trt1_effect[j])
					{
						outfile << p0.qtl_afs_trt0[i].freq[k] << ",";
						onlist = true;
					}
				}
				if (!onlist)
					outfile << "0" << ",";

				onlist = false;
				for (k = 0; k < p1.qtl_afs_trt0[i].freq.size(); k++)
				{
					if (p1.qtl_afs_trt0[i].trt0_effect[k] == combo.qtl_afs_trt0[i].trt0_effect[j] &&
						p1.qtl_afs_trt0[i].trt1_effect[k] == combo.qtl_afs_trt0[i].trt1_effect[j])
					{
						outfile << p1.qtl_afs_trt0[i].freq[k] << ",";
						onlist = true;
					}
				}
				if (!onlist)
					outfile << "0" << ",";

				outfile << combo.qtl_afs_trt0[i].freq[j];
			}
		}

		for (i = 0; i < combo.qtl_afs_trt1.size(); i++)
		{
			for (j = 0; j < combo.qtl_afs_trt1[i].freq.size(); j++)
			{
				outfile << "\nQ_" << i << "," << combo.qtl_loc_trt1[i] << ",";
				outfile << "All_" << j << ",";
				outfile << combo.qtl_afs_trt1[i].trt0_effect[j] << ",";
				outfile << combo.qtl_afs_trt1[i].trt1_effect[j] << ",";

				onlist = false;
				for (k = 0; k < p0.qtl_afs_trt1[i].freq.size(); k++)
				{
					if (p0.qtl_afs_trt1[i].trt0_effect[k] == combo.qtl_afs_trt1[i].trt0_effect[j] &&
						p0.qtl_afs_trt1[i].trt1_effect[k] == combo.qtl_afs_trt1[i].trt1_effect[j])
					{
						outfile << p0.qtl_afs_trt1[i].freq[k] << ",";
						onlist = true;
					}
				}
				if (!onlist)
					outfile << "0" << ",";

				onlist = false;
				for (k = 0; k < p1.qtl_afs_trt1[i].freq.size(); k++)
				{
					if (p1.qtl_afs_trt1[i].trt0_effect[k] == combo.qtl_afs_trt1[i].trt0_effect[j] &&
						p1.qtl_afs_trt1[i].trt1_effect[k] == combo.qtl_afs_trt1[i].trt1_effect[j])
					{
						outfile << p1.qtl_afs_trt1[i].freq[k] << ",";
						onlist = true;
					}
				}
				if (!onlist)
					outfile << "0" << ",";

				outfile << combo.qtl_afs_trt1[i].freq[j];
			}
		}

		for (i = 0; i < combo.qtl_afs_pleio.size(); i++)
		{
			for (j = 0; j < combo.qtl_afs_pleio[i].freq.size(); j++)
			{
				outfile << "\nQ_" << i << "," << combo.qtl_loc_pleio[i] << ",";
				outfile << "All_" << j << ",";
				outfile << combo.qtl_afs_pleio[i].trt0_effect[j] << ",";
				outfile << combo.qtl_afs_pleio[i].trt1_effect[j] << ",";

				onlist = false;
				for (k = 0; k < p0.qtl_afs_pleio[i].freq.size(); k++)
				{
					if (p0.qtl_afs_pleio[i].trt0_effect[k] == combo.qtl_afs_pleio[i].trt0_effect[j] &&
						p0.qtl_afs_pleio[i].trt1_effect[k] == combo.qtl_afs_pleio[i].trt1_effect[j])
					{
						outfile << p0.qtl_afs_pleio[i].freq[k] << ",";
						onlist = true;
					}
				}
				if (!onlist)
					outfile << "0" << ",";

				onlist = false;
				for (k = 0; k < p1.qtl_afs_pleio[i].freq.size(); k++)
				{
					if (p1.qtl_afs_pleio[i].trt0_effect[k] == combo.qtl_afs_pleio[i].trt0_effect[j] &&
						p1.qtl_afs_pleio[i].trt1_effect[k] == combo.qtl_afs_pleio[i].trt1_effect[j])
					{
						outfile << p1.qtl_afs_pleio[i].freq[k] << ",";
						onlist = true;
					}
				}
				if (!onlist)
					outfile << "0" << ",";

				outfile << combo.qtl_afs_pleio[i].freq[j];
			}
		}

		outfile.close();
	}

	void calc_and_save_fst_values_markers(chromosomedata &p0, chromosomedata &p1, chromosomedata &combo, std::string filename, std::vector<double> &fst_list, bool save_to_file)
	{
		int i, j;
		double temp_fst, h_t, h_s0, h_s1;
		std::ofstream outfile;

		if (save_to_file)
		{
			char temp_fn[256];
			convert_string_to_char(filename, temp_fn);
			outfile.open(temp_fn);
		}
		// First Calculate Fst for Marker Loci

		combo.marker_fst_values.clear();
		combo.marker_heterozygosity.clear();

		if (save_to_file)
			outfile << "Marker,Fst";

		for (i = 0; i < combo.marker_afs.size(); i++)
		{
			h_t = 1;
			for (j = 0; j < combo.marker_afs[i].freq.size(); j++)
			{
				h_t = h_t - (combo.marker_afs[i].freq[j] * combo.marker_afs[i].freq[j]);
			}

			h_s0 = 1;
			for (j = 0; j < p0.marker_afs[i].freq.size(); j++)
			{
				h_s0 = h_s0 - (p0.marker_afs[i].freq[j] * p0.marker_afs[i].freq[j]);
			}

			h_s1 = 1;
			for (j = 0; j < p1.marker_afs[i].freq.size(); j++)
			{
				h_s1 = h_s1 - (p1.marker_afs[i].freq[j] * p1.marker_afs[i].freq[j]);
			}

			if (h_t > 0)
				temp_fst = (h_t - (h_s0 + h_s1) / 2) / h_t;
			else
				temp_fst = 0;

			if (save_to_file)
				outfile << "\n" << i << "," << temp_fst;
			
			fst_list.push_back(temp_fst);
			combo.marker_fst_values.push_back(temp_fst);
			combo.marker_heterozygosity.push_back(h_t);
		}

		if (save_to_file)
			outfile.close();

	}

	void calc_and_save_fst_values_qtl(chromosomedata &p0, chromosomedata &p1, chromosomedata &combo, std::string filename, std::vector<double> &fst_list, bool save_to_file)
	{
		int i, j;
		double temp_fst, h_t, h_s0, h_s1;

		std::ofstream outfile;

		if (save_to_file)
		{
			char temp_fn[256];
			convert_string_to_char(filename, temp_fn);
			outfile.open(temp_fn);
			outfile << "Type,Locus,Fst";
		}

		combo.qtl_fst_values_trt0.clear();
		combo.qtl_fst_values_trt1.clear();
		combo.qtl_fst_values_pleio.clear();

		// Trait 0 loci
		for (i = 0; i < combo.qtl_afs_trt0.size(); i++)
		{
			h_t = 1;
			for (j = 0; j < combo.qtl_afs_trt0[i].freq.size(); j++)
			{
				h_t = h_t - (combo.qtl_afs_trt0[i].freq[j]*combo.qtl_afs_trt0[i].freq[j]);
			}

			h_s0 = 1;
			for (j = 0; j < p0.qtl_afs_trt0[i].freq.size(); j++)
			{
				h_s0 = h_s0 - (p0.qtl_afs_trt0[i].freq[j] * p0.qtl_afs_trt0[i].freq[j]);
			}

			h_s1 = 1;
			for (j = 0; j < p1.qtl_afs_trt0[i].freq.size(); j++)
			{
				h_s1 = h_s1 - (p1.qtl_afs_trt0[i].freq[j] * p1.qtl_afs_trt0[i].freq[j]);
			}

			if (h_t > 0)
				temp_fst = (h_t - (h_s0 + h_s1) / 2) / h_t;
			else
				temp_fst = 0;

			if (save_to_file)
				outfile << "\nTrt0," << combo.qtl_loc_trt0[i] << "," << temp_fst;
			
			fst_list.push_back(temp_fst);
			combo.qtl_fst_values_trt0.push_back(temp_fst);
		}

		// Trait 1 loci
		for (i = 0; i < combo.qtl_afs_trt1.size(); i++)
		{
			h_t = 1;
			for (j = 0; j < combo.qtl_afs_trt1[i].freq.size(); j++)
			{
				h_t = h_t - (combo.qtl_afs_trt1[i].freq[j] * combo.qtl_afs_trt1[i].freq[j]);
			}

			h_s0 = 1;
			for (j = 0; j < p0.qtl_afs_trt1[i].freq.size(); j++)
			{
				h_s0 = h_s0 - (p0.qtl_afs_trt1[i].freq[j] * p0.qtl_afs_trt1[i].freq[j]);
			}

			h_s1 = 1;
			for (j = 0; j < p1.qtl_afs_trt1[i].freq.size(); j++)
			{
				h_s1 = h_s1 - (p1.qtl_afs_trt1[i].freq[j] * p1.qtl_afs_trt1[i].freq[j]);
			}

			if (h_t > 0)
				temp_fst = (h_t - (h_s0 + h_s1) / 2) / h_t;
			else
				temp_fst = 0;

			if (save_to_file)
				outfile << "\nTrt1," << combo.qtl_loc_trt1[i] << "," << temp_fst;
			fst_list.push_back(temp_fst);
			combo.qtl_fst_values_trt1.push_back(temp_fst);
		}

		// Pleitropic loci
		for (i = 0; i < combo.qtl_afs_pleio.size(); i++)
		{
			h_t = 1;
			for (j = 0; j < combo.qtl_afs_pleio[i].freq.size(); j++)
			{
				h_t = h_t - (combo.qtl_afs_pleio[i].freq[j] * combo.qtl_afs_pleio[i].freq[j]);
			}

			h_s0 = 1;
			for (j = 0; j < p0.qtl_afs_pleio[i].freq.size(); j++)
			{
				h_s0 = h_s0 - (p0.qtl_afs_pleio[i].freq[j] * p0.qtl_afs_pleio[i].freq[j]);
			}

			h_s1 = 1;
			for (j = 0; j < p1.qtl_afs_pleio[i].freq.size(); j++)
			{
				h_s1 = h_s1 - (p1.qtl_afs_pleio[i].freq[j] * p1.qtl_afs_pleio[i].freq[j]);
			}

			if (h_t > 0)
				temp_fst = (h_t - (h_s0 + h_s1) / 2) / h_t;
			else
				temp_fst = 0;

			if (save_to_file)
				outfile << "\nPleio," << combo.qtl_loc_pleio[i] << "," << temp_fst;

			fst_list.push_back(temp_fst);
			combo.qtl_fst_values_pleio.push_back(temp_fst);
		}

		if (save_to_file)
			outfile.close();
	}

	void smooth_fst_list(int n_chrom, int n_markers, double weightingvariance, std::vector<double> &fst_list, std::vector<double> &smoothed_fst_list)
	{
		int iii, jjj, kkk;
		double weightingcoefficient;
		int weightingdistance;
		double weightedFst;
		double numbersummed;

		weightingdistance = static_cast<int>(floor(sqrt(weightingvariance)) * 3);

		int chrom_beg;
		int chrom_end;
		for (iii = 0; iii < n_chrom; iii++)
		{
			chrom_beg = iii * n_markers;
			chrom_end = iii * n_markers + n_markers;
			for (jjj = chrom_beg; jjj < chrom_end; jjj++)
			{
				weightedFst = 0;
				numbersummed = 0;
				for (kkk = (jjj - weightingdistance); kkk < jjj + weightingdistance + 1; kkk++)
				{
					if (kkk >= chrom_beg && kkk < chrom_end)
					{
						if (fst_list[iii] >= 0)
						{
							weightingcoefficient = exp(-0.5*(kkk - jjj)*(kkk - jjj) / weightingvariance);
							weightedFst = weightedFst + weightingcoefficient * fst_list[kkk];
							numbersummed = numbersummed + weightingcoefficient;
						}
					}
				}
				smoothed_fst_list.push_back(weightedFst / numbersummed);
			}
		}
	}

	void smooth_fst_values_per_chromosome(int n_markers, double weightingvariance, chromosomedata &combo)
	{
		int jjj, kkk;
		double weightingcoefficient;
		int weightingdistance;
		double weightedFst;
		double numbersummed;

		weightingdistance = static_cast<int>(floor(sqrt(weightingvariance)) * 3);

		combo.smoothed_marker_fsts.clear();
		for (jjj = 0; jjj < n_markers; jjj++)
		{
			weightedFst = 0;
			numbersummed = 0;
			for (kkk = (jjj - weightingdistance); kkk < jjj + weightingdistance + 1; kkk++)
				{
					if (kkk >= 0 && kkk < n_markers)
					{
						if (combo.marker_fst_values[kkk] >= 0)
						{
							weightingcoefficient = exp(-0.5*(kkk - jjj)*(kkk - jjj) / weightingvariance);
							weightedFst = weightedFst + weightingcoefficient * combo.marker_fst_values[kkk];
							numbersummed = numbersummed + weightingcoefficient;
						}
					}
				}
				combo.smoothed_marker_fsts.push_back(weightedFst / numbersummed);
		}
	}

	void smoothed_fst_critical_values_by_chromosome(int n_markers, chromosomedata &chrom_dat)
	{
		int jjj;
		double meanFst, varFst, stdevFst;
		double count;

		meanFst = 0;
		count = 0;
		
		for (jjj = 0; jjj < n_markers; jjj++)
		{
			if (chrom_dat.smoothed_marker_fsts[jjj] >= 0)
			{
				meanFst = meanFst + chrom_dat.smoothed_marker_fsts[jjj];
				count++;
			}
		}
		meanFst = meanFst / count;

		varFst = 0;
		for (jjj = 0; jjj < n_markers; jjj++)
		{
			if (chrom_dat.smoothed_marker_fsts[jjj] >= 0)
			{
				varFst = varFst + (meanFst - chrom_dat.smoothed_marker_fsts[jjj])*(meanFst - chrom_dat.smoothed_marker_fsts[jjj]);
			}
		}
		
		varFst = varFst / count;
		stdevFst = sqrt(varFst);

		chrom_dat.SmoothedCritValue99 = meanFst + 2.57583 * stdevFst;
		chrom_dat.SmoothedCritValue98 = meanFst + 2.32635 * stdevFst;
		chrom_dat.SmoothedCritValue95 = meanFst + 1.95996 * stdevFst;
		chrom_dat.SmoothedCritValue90 = meanFst + 1.64485 * stdevFst;
		chrom_dat.SmoothedCritValue80 = meanFst + 1.28155 * stdevFst;
	}

	void fst_chi_sqr_per_chromosome(int n_markers, int total_sample_size, chromosomedata &chrom_dat)
	{
		int i, j;
		double actualnumberofalleles;
		double degreesoffreedom;
		double pvalue;
		double chisqrteststat;
		double Ntotalsamplesize;
		Ntotalsamplesize = total_sample_size;

		// Clear the chi square vector
		chrom_dat.fst_chi_sqr_p.clear();

		for (i = 0; i < n_markers; i++) // loop through all loci
		{
			// Need to know the actual number of alleles per locus
			actualnumberofalleles = 0;
			for (j = 0; j < chrom_dat.marker_afs[i].freq.size(); j++)
			{
				if (chrom_dat.marker_afs[i].freq[j] > 0)
					actualnumberofalleles++;
			}

			if (chrom_dat.marker_fst_values[i] > 0)
			{
				chisqrteststat = 2 * Ntotalsamplesize*chrom_dat.marker_fst_values[i] * (actualnumberofalleles - 1);
				degreesoffreedom = (actualnumberofalleles - 1);
				pvalue = 1 - chisqr(chisqrteststat, degreesoffreedom);
			}
			else
				pvalue = 1.01;

			chrom_dat.fst_chi_sqr_p.push_back(pvalue);
		}

	}

	void FDRchisqrPval(double fdr_alpha, int n_markers, chromosomedata &chrom_dat)
	{
		int iii, jjj;
		double totalvariablemarkers;
		std::vector<double> pvals;
		pvals.reserve(n_markers);

		for (jjj = 0; jjj < n_markers; jjj++)
		{
			pvals.push_back(chrom_dat.fst_chi_sqr_p[jjj]);
		}

		sort(pvals.begin(), pvals.end());

		totalvariablemarkers = 0;
		for (iii = 0; iii < n_markers; iii++)
		{
			if (pvals[iii] != 1.01)
				totalvariablemarkers++;
		}

		double pvalindex, tempFDRp;
		for (iii = static_cast<int>(totalvariablemarkers) - 1; iii >= 0; iii = iii - 1)
		{
			pvalindex = iii + 1;
			tempFDRp = (pvalindex / totalvariablemarkers)*fdr_alpha;
			if (pvals[iii] <= tempFDRp)
				break;
		}

		chrom_dat.FDRp = tempFDRp;
	}

	void determine_significant_loci_chisqr(chromosomedata &chrom_dat)
	{
		int i;
		chrom_dat.is_significant.clear();
		for (i = 0; i < chrom_dat.fst_chi_sqr_p.size(); i++)
		{
			if (chrom_dat.fst_chi_sqr_p[i] > chrom_dat.FDRp)
				chrom_dat.is_significant.push_back(false);
			else
				chrom_dat.is_significant.push_back(true);
		}
	}

	void calc_FSTprime_per_chrom(chromosomedata &chrom_dat)
	{
		// Use only loci with a heterozygosity > 0.05 for this analysis.
		// We're using an H(T) of 0.05 or greater.
		// The regular FST has to be calculated already, because that
		// is also when the HT is calculated.

		int i;
		double mean_fst;
		double n_markers;

		mean_fst = 0;
		n_markers = 0;
		for (i = 0; i < chrom_dat.marker_fst_values.size(); i++)
		{
			if (chrom_dat.marker_heterozygosity[i] >= 0.05)
			{
				mean_fst = mean_fst + chrom_dat.marker_fst_values[i];
				n_markers++;
			}
		}
		mean_fst = mean_fst / n_markers;

		chrom_dat.marker_fst_prime_value.clear();
		for (i = 0; i < chrom_dat.marker_fst_values.size(); i++)
		{
			chrom_dat.marker_fst_prime_value.push_back(chrom_dat.marker_fst_values[i]/mean_fst);
		}

		// Calc p-values
		chrom_dat.fst_prime_p_value.clear();
		for (i = 0; i < chrom_dat.marker_fst_prime_value.size(); i++)
		{
			chrom_dat.fst_prime_p_value.push_back(1-chisqr(chrom_dat.marker_fst_prime_value[i], 1));
		}

		// Calc FDR crit val
		int iii, jjj;
		double totalvariablemarkers;
		std::vector<double> pvals;

		for (jjj = 0; jjj < chrom_dat.fst_prime_p_value.size(); jjj++)
			pvals.push_back(chrom_dat.fst_prime_p_value[jjj]);

		sort(pvals.begin(), pvals.end());

		totalvariablemarkers = 0;
		for (iii = 0; iii < chrom_dat.marker_heterozygosity.size(); iii++)
		{
			if (chrom_dat.marker_heterozygosity[iii] >= 0.05)
				totalvariablemarkers++;
		}

		double pvalindex, tempFDRp;
		for (iii = static_cast<int>(totalvariablemarkers) - 1; iii >= 0; iii = iii - 1)
		{
			pvalindex = iii + 1;
			tempFDRp = (pvalindex / totalvariablemarkers)*0.05;
			if (pvals[iii] <= tempFDRp)
				break;
		}

		// Determine which loci are significant at FDR 0.05
		chrom_dat.fst_prime_is_significant_point05.clear();
		chrom_dat.fst_prime_is_significant_point01.clear();
		chrom_dat.fst_prime_is_significant_fdrpoint05.clear();
		for (i = 0; i < chrom_dat.fst_prime_p_value.size(); i++)
		{
			if (chrom_dat.fst_prime_p_value[i] > 0.05)
				chrom_dat.fst_prime_is_significant_point05.push_back(false);
			else
				chrom_dat.fst_prime_is_significant_point05.push_back(true);

			if (chrom_dat.fst_prime_p_value[i] > 0.01)
				chrom_dat.fst_prime_is_significant_point01.push_back(false);
			else
				chrom_dat.fst_prime_is_significant_point01.push_back(true);

			if (chrom_dat.fst_prime_p_value[i] > tempFDRp)
				chrom_dat.fst_prime_is_significant_fdrpoint05.push_back(false);
			else
				chrom_dat.fst_prime_is_significant_fdrpoint05.push_back(true);
		}

	}

	int run_simulation(std::string outfile_name)
	{
		std::stringstream sss;
		int ig, i, generations;
		int num_chro, num_markers_per, num_qtl0_per;
		int num_qtl1_per, num_qtlp_per;
		int j, k;
		double fst_wv;
		
		// Initialize all of the populations
		for (i = 0; i < Number_of_populations; i++)
		{
			population[i].initialize_population();
		}

		// Set the positions of the QTLs to be the same in the two pops
		std::vector<int> qtl_locations;
		qtl_locations = population[0].report_qtl_locations();
		population[1].set_qtl_locations(qtl_locations);

		// Initialize Chromosome Data Vectors
		num_chro = population[0].report_number_chromosomes();
		num_markers_per = population[0].report_markers_per_chromosome();
		num_qtl0_per = population[0].report_qtl0_per_chromosome();
		num_qtl1_per = population[0].report_qtl1_per_chromosome();
		num_qtlp_per = population[0].report_qtlP_per_chromosome();
		fst_wv = population[0].report_fst_weight_var();

		// Prepare vectors for fst data
		std::vector<marker_fst_data> markerfsts;
		std::vector<qtl_fst_data> qtlfsts;
		std::vector<double> marker_fst_list;
		std::vector<double> qtl_fst_list;
		std::vector<double> smoothed;

		k = 0;
		markerfsts.resize(num_chro*num_markers_per);
		for (i = 0; i < num_chro; i++)
		{
			for (j = 0; j < num_markers_per; j++)
			{
				markerfsts[k].chromosome = i;
				markerfsts[k].location = j;
				k++;
			}
		}

		int total_qtl_per = num_qtl0_per + num_qtl1_per + num_qtlp_per;
		qtlfsts.resize(num_chro*total_qtl_per);
		k = 0;
		for (i = 0; i < num_chro; i++)
		{
			for (j = 0; j < num_qtl0_per; j++)
			{
				qtlfsts[k].chromosome = i;
				qtlfsts[k].type = 0;
				qtlfsts[k].location = qtl_locations[k];
				k++;
			}
			for (j = 0; j < num_qtl1_per; j++)
			{
				qtlfsts[k].chromosome = i;
				qtlfsts[k].type = 1;
				qtlfsts[k].location = qtl_locations[k];
				k++;
			}
			for (j = 0; j < num_qtlp_per; j++)
			{
				qtlfsts[k].chromosome = i;
				qtlfsts[k].type = 2;
				qtlfsts[k].location = qtl_locations[k];
				k++;
			}
		}


		// Run initial generations
		for (ig = 0; ig < population[0].getNumberOfInitialGenerations(); ig++)
		{
			for (i = 0; i < Number_of_populations; i++)
				population[i].polygynous_mating();
			for (i = 0; i < Number_of_populations; i++)
				population[i].mutation();
			progeny_stepping_stone_migration(population[0], population[1]);
			for (i = 0; i < Number_of_populations; i++)
				population[i].selection(true);
			for (i = 0; i < Number_of_populations; i++)
				population[i].population_regulation();
			/*for (i = 0; i < Number_of_populations; i++)
			{
				std::string temp_s;
				sss.str("");
				sss << outfile_name << "_genotypes_pop_" << i + 1 << "_init_gen_" << ig << ".csv";
				temp_s = sss.str();
				population[i].save_all_genotypes(temp_s);
			}*/
		}

		// Run the intermediate generations
		for (ig = 0; ig < population[0].getNumberInterveningGens(); ig++)
		{
			for (i = 0; i < Number_of_populations; i++)
				population[i].polygynous_mating();
			for (i = 0; i < Number_of_populations; i++)
				population[i].mutation();
			progeny_stepping_stone_migration(population[0], population[1]);
			for (i = 0; i < Number_of_populations; i++)
				population[i].selection(false);
			for (i = 0; i < Number_of_populations; i++)
				population[i].population_regulation();
		}

		// Turn the lifecycle halfway
		for (i = 0; i < Number_of_populations; i++)
			population[i].calculate_values_progeny();
		for (i = 0; i < Number_of_populations; i++)
			population[i].polygynous_mating();
		for (i = 0; i < Number_of_populations; i++)
			population[i].mutation();
		progeny_stepping_stone_migration(population[0], population[1]);

		// Start the experimental generations
		for (generations = 0; generations < population[0].getNumberOfGenerations(); generations++)
		{
			for (i = 0; i < Number_of_populations; i++)
				population[i].selection(false);
			for (i = 0; i < Number_of_populations; i++)
				population[i].calculate_values_progeny();
			for (i = 0; i < Number_of_populations; i++)
				population[i].estimate_additive_variances_progeny();
			for (i = 0; i < Number_of_populations; i++)
				population[i].population_regulation();

			/*for (i = 0; i < Number_of_populations; i++)
			{
				std::string temp_s;
				sss.str("");
				sss << outfile_name << "_genotypes_pop_" << i + 1 << "_exp_gen_" << generations << ".csv";
				temp_s = sss.str();
				population[i].save_all_genotypes(temp_s);
			}*/

			for (i = 0; i < Number_of_populations; i++)
				population[i].gaussian_mating();
			for (i = 0; i < Number_of_populations; i++)
				population[i].mutation();
			for (i = 0; i < Number_of_populations; i++)
				population[i].calculate_values_adults();
			for (i = 0; i < Number_of_populations; i++)
				population[i].store_variables_in_memory();
			progeny_stepping_stone_migration(population[0], population[1]);

			// Every "all_calc_int" experimental generations,
			// take a sample of adults and calculate allele frequencies
			// and Fst values

			if ((generations + 1) % population[0].report_all_calc_int() == 0 || generations == 0 || generations + 1 == population[0].getNumberOfGenerations())
			{
				chrom_data_pop0.clear();
				chrom_data_pop1.clear();
				for (i = 0; i < Number_of_populations; i++)
				{
					population[i].sample_adults();
					population[i].calculate_allele_frequencies_all_chromosomes();
					population[i].calc_qtl_means_vars_all_chromosomes();

					for (j = 0; j < num_chro; j++)
					{
						if (i == 0)
							chrom_data_pop0.push_back(population[0].report_chromosome_data(j));
						if (i == 1)
							chrom_data_pop1.push_back(population[1].report_chromosome_data(j));
					}

					std::string temp_fn;
					sss.str("");
					sss << outfile_name << "_marker_af_pop_" << i << "_exp_gen_" << generations + 1 << ".csv";
					temp_fn = sss.str();

					population[i].save_marker_freqs_to_file(temp_fn);

					sss.str("");
					sss << outfile_name << "_qtl_af_pop_" << i << "_exp_gen_" << generations + 1 << ".csv";
					temp_fn = sss.str();

					population[i].save_qtl_freqs_to_file(temp_fn);

					sss.str("");
					sss << outfile_name << "_LD_pop_" << i << "_exp_gen_" << generations + 1 << ".csv";
					temp_fn = sss.str();

					population[i].save_pairwise_LD(temp_fn);
				}

				// Combine chromosome data and output it
				
				chrom_data_combined.clear();
				for (j = 0; j < num_chro; j++)
				{
					chrom_data_combined.push_back(combine_chromosome_data(chrom_data_pop0[j], chrom_data_pop1[j]));
				}

				// Save to file
				marker_fst_list.clear();
				qtl_fst_list.clear();
				smoothed.clear();
				for (j = 0; j < num_chro; j++)
				{
					std::string temp_fn;
					sss.str("");
					sss << outfile_name << "_all_marker_allelefreqs_chrom_" << j << "_exp_gen_" << generations + 1 << ".csv";
					temp_fn = sss.str();

					save_combined_marker_data(chrom_data_pop0[j], chrom_data_pop1[j], chrom_data_combined[j], temp_fn);
				}

				for (j = 0; j < num_chro; j++)
				{
					std::string temp_fn;
					sss.str("");
					sss << outfile_name << "_all_qtl_allelefreqs_chrom_" << j << "_exp_gen_" << generations + 1 << ".csv";
					temp_fn = sss.str();

					save_combined_qtl_data(chrom_data_pop0[j], chrom_data_pop1[j], chrom_data_combined[j], temp_fn);
				}

				for (j = 0; j < num_chro; j++)
				{
					std::string temp_fn;
					sss.str("");
					sss << outfile_name << "_marker_fst_chrom_" << j << "_exp_gen_" << generations + 1 << ".csv";
					temp_fn = sss.str();

					calc_and_save_fst_values_markers(chrom_data_pop0[j], chrom_data_pop1[j], chrom_data_combined[j], temp_fn, marker_fst_list, false);
				}

				for (j = 0; j < num_chro; j++)
				{
					std::string temp_fn;
					sss.str("");
					sss << outfile_name << "_qtl_fst_chrom_" << j << "_exp_gen_" << generations + 1 << ".csv";
					temp_fn = sss.str();

					calc_and_save_fst_values_qtl(chrom_data_pop0[j], chrom_data_pop1[j], chrom_data_combined[j], temp_fn, qtl_fst_list, false);
				}

				// Generate smoothed Fst values

				smooth_fst_list(num_chro, num_markers_per, fst_wv, marker_fst_list, smoothed);
				for (j = 0; j < num_chro; j++)
					smooth_fst_values_per_chromosome(num_markers_per, fst_wv, chrom_data_combined[j]);

				// Move Fst Lists into Fst storage classes
				
				for (j = 0; j < num_chro*num_markers_per; j++)
				{
					markerfsts[j].fst = marker_fst_list[j];
					markerfsts[j].smoothed = smoothed[j];
				}
				for (j = 0; j < num_chro*total_qtl_per; j++)
					qtlfsts[j].fst = qtl_fst_list[j];

				// Determine FST smoothed critical values for each chromosome
				for (i = 0; i < num_chro; i++)
				{
					smoothed_fst_critical_values_by_chromosome(num_markers_per, chrom_data_combined[i]);
				}

				// Calculate chi-squared p-value for each marker locus
				int tot_sampl_size = population[0].report_sample_size() + population[1].report_sample_size();
				for (i = 0; i < num_chro; i++)
				{
					fst_chi_sqr_per_chromosome(num_markers_per, tot_sampl_size, chrom_data_combined[i]);
				}

				// Determine FDR chi-sqr p-value cutoff
				for (i = 0; i < num_chro; i++)
					FDRchisqrPval(0.05, num_markers_per, chrom_data_combined[i]);

				// Determine Chi-Square significant loci
				for (i = 0; i < num_chro; i++)
					determine_significant_loci_chisqr(chrom_data_combined[i]);

				// Calc FSTprime
				for (i = 0; i < num_chro; i++)
				{
					calc_FSTprime_per_chrom(chrom_data_combined[i]);
				}

				// Find smoothed Fst peaks
				for (i = 0; i < num_chro; i++)
					chrom_data_combined[i].find_all_local_smoothed_fst_maxima(30);

				// Compile Peak Data
				for (i = 0; i < num_chro; i++)
					chrom_data_combined[i].compile_peak_data(10);

				// Save these combined fst lists to files
				std::string temp_fn;
				char temp_fn_ch[256];
				sss.str("");
				sss << outfile_name << "_all_marker_fsts_exp_gen_" << generations + 1 << ".csv";
				temp_fn = sss.str();

				std::ofstream ootfile;
				convert_string_to_char(temp_fn, temp_fn_ch);
				ootfile.open(temp_fn_ch);
				ootfile << "Chrom,Locus,Genome_Location,Fst,Smoothed,Pvalue,Significant,H_T,Genome_Location,FSTprime,FSTprimeP,FSTprimeSig05,FSTprimeSig01,FDR05sigFSTprime";
				int counter = 0;
				for (j = 0; j < num_chro; j++)
				{
					for (k = 0; k < num_markers_per; k++)
					{
						ootfile << "\n";
						ootfile << j << ",";
						ootfile << k << ",";
						ootfile << k + j * num_markers_per << ",";
						ootfile << chrom_data_combined[j].marker_fst_values[k] << ",";
						ootfile << chrom_data_combined[j].smoothed_marker_fsts[k] << ",";
						ootfile << chrom_data_combined[j].fst_chi_sqr_p[k] << ",";
						ootfile << chrom_data_combined[j].is_significant[k] << ",";
						ootfile << chrom_data_combined[j].marker_heterozygosity[k] << ",";
						ootfile << k + j * num_markers_per << ",";
						ootfile << chrom_data_combined[j].marker_fst_prime_value[k] << ",";
						ootfile << chrom_data_combined[j].fst_prime_p_value[k] << ",";
						ootfile << chrom_data_combined[j].fst_prime_is_significant_point05[k] << ",";
						ootfile << chrom_data_combined[j].fst_prime_is_significant_point01[k] << ",";
						ootfile << chrom_data_combined[j].fst_prime_is_significant_fdrpoint05[k];
						counter++;
					}
				}
				ootfile.close();

				sss.str("");
				sss << outfile_name << "_all_qtl_fsts_exp_gen_" << generations + 1 << ".csv";
				temp_fn = sss.str();
				convert_string_to_char(temp_fn, temp_fn_ch);
				ootfile.open(temp_fn_ch);
				ootfile << "Trait,Chrom,Locus,Genome_Location,Fst,Pop0trt0mean,Pop1trt0mean,Pop0trt1mean,Pop1trt1mean,Pop0trt0var,Pop1trt0var,Pop0trt1var,Pop1trt1var";
				
				for (j = 0; j < num_chro; j++)
				{
					for (k = 0; k < num_qtl0_per; k++)
					{
						ootfile << "\n";
						ootfile << "Trait_0,";
						ootfile << j << ",";
						ootfile << chrom_data_combined[j].qtl_loc_trt0[k] << ",";
						ootfile << j * num_markers_per + chrom_data_combined[j].qtl_loc_trt0[k] << ",";
						ootfile << chrom_data_combined[j].qtl_fst_values_trt0[k] << ",";
						ootfile << chrom_data_pop0[j].qtl_mean_trt0[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_mean_trt0[k] << ",";
						ootfile << "0,0,";
						ootfile << chrom_data_pop0[j].qtl_var_trt0[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_var_trt0[k] << ",";
						ootfile << "0,0";
					}
					for (k = 0; k < num_qtl1_per; k++)
					{
						ootfile << "\n";
						ootfile << "Trait_1,";
						ootfile << j << ",";
						ootfile << chrom_data_combined[j].qtl_loc_trt1[k] << ",";
						ootfile << j * num_markers_per + chrom_data_combined[j].qtl_loc_trt1[k] << ",";
						ootfile << chrom_data_combined[j].qtl_fst_values_trt1[k] << ",";
						ootfile << "0,0,";
						ootfile << chrom_data_pop0[j].qtl_mean_trt1[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_mean_trt1[k] << ",";
						ootfile << "0,0,";
						ootfile << chrom_data_pop0[j].qtl_var_trt1[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_var_trt1[k];
					}
					for (k = 0; k < num_qtlp_per; k++)
					{
						ootfile << "\n";
						ootfile << "Trait_Pleio,";
						ootfile << j << ",";
						ootfile << chrom_data_combined[j].qtl_loc_pleio[k] << ",";
						ootfile << j * num_markers_per + chrom_data_combined[j].qtl_loc_pleio[k] << ",";
						ootfile << chrom_data_combined[j].qtl_fst_values_pleio[k] << ",";
						ootfile << chrom_data_pop0[j].qtl_mean_pleio_0[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_mean_pleio_0[k] << ",";
						ootfile << chrom_data_pop0[j].qtl_mean_pleio_1[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_mean_pleio_1[k] << ",";
						ootfile << chrom_data_pop0[j].qtl_var_pleio_0[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_var_pleio_0[k] << ",";
						ootfile << chrom_data_pop0[j].qtl_var_pleio_1[k] << ",";
						ootfile << chrom_data_pop1[j].qtl_var_pleio_1[k];
					}
				}
				ootfile.close();

				// Save Peak data
				sss.str("");
				sss << outfile_name << "_all_peak_data_exp_gen_" << generations + 1 << ".csv";
				temp_fn = sss.str();
				convert_string_to_char(temp_fn, temp_fn_ch);
				ootfile.open(temp_fn_ch);
				ootfile << "Chrom,Locus,Genome_Location,smoothedFst,highestFst,highestFstLoc,sigFDR,sig80CI,sig90CI,";
				ootfile << "sig95CI,sig98CI,sig99CI,sigFSTprime05,sigFSTprime01,sigFSTfdr05,nearestQTL0,distQTL0,nearestQTL1,distQTL1,nearestQTLP,distQTLP";

				for (i = 0; i < num_chro; i++)
				{
					for (j = 0; j < chrom_data_combined[i].peakdata.numberofpeaks; j++)
					{
						ootfile << "\n" << i << ",";
						ootfile << chrom_data_combined[i].peakdata.peaklocation[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.peaklocation[j] + i * num_markers_per << ",";
						ootfile << chrom_data_combined[i].peakdata.peaksmoothedFst[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.highestFstonplateau[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.highestFstIndex[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sigFDR[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sig80CI[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sig90CI[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sig95CI[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sig98CI[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sig99CI[j] << ",";

						ootfile << chrom_data_combined[i].peakdata.sig05FSTprime[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sig01FSTprime[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.sigFDRFSTprime[j] << ",";

						ootfile << chrom_data_combined[i].peakdata.nearestQTL_0[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.distancetoQTL_0[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.nearestQTL_1[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.distancetoQTL_1[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.nearestQTL_P[j] << ",";
						ootfile << chrom_data_combined[i].peakdata.distancetoQTL_P[j] << ",";
					}
				}

				ootfile.close();
	
			}

		} // generations loop

		for (i = 0; i < Number_of_populations; i++)
		{
			std::string temp_fn;
			sss.str("");
			sss << outfile_name << "_output_pop_" << i << ".csv";
			temp_fn = sss.str();
			population[i].save_stored_variables(temp_fn);
			population[i].calc_run_means();
			population[i].append_means_to_output_file(temp_fn);
		}

		// Compile variables of interest


		std::string out_header;
		out_header = "Migr_Rate,No_gens,N,Sample_Size,N_Chrom,N_markers_per_chrom,Mut_Rate_per_Marker,n_qtl_trt0,n_qtl_trt1,n_qtl_p,mut_rate_qtl,env_var0,env_var1,";
		out_header = out_header + "recomb_per_chrom,omega0,omega1,r(omega),opt0_pop0,opt1_pop0,opt0_pop1,opt1_pop1,";
		out_header = out_header + "zbar0_pop0,zbar1_pop0,G00_pop0,G11_pop0,G01_pop0,V00_pop0,V11_pop0,V01_pop0,";
		out_header = out_header + "zbar0_pop1,zbar1_pop1,G00_pop1,G11_pop1,G01_pop1,V00_pop1,V11_pop1,V01_pop1,";
		out_header = out_header + "mean_marker_Fst,n_signif_peaks,n_within_25_loci_qtl_trt0,n_within_25_loci_qtl_trt1,n_within_25_loci_qtl_p,mean_qtl0_fst,mean_qtl1_fst,mean_qtlp_fst,";
		out_header = out_header + "n_sign_pks_fstprime01,n_near_qtl0_fstprime,n_near_qtl1_fst_prime,n_near_qtlp_fstprime";
		std::vector<double> output_line;
		output_line.push_back(migration_rate);
		population[0].report_parameter_vals(output_line);
		population[1].report_optima(output_line);
		population[0].report_means(output_line);
		population[1].report_means(output_line);

		// Calc_mean_marker_fst
		double mean_marker_fst = 0;
		double d_count = 0;
		for (i = 0; i < num_chro; i++)
		{ 
			for (j = 0; j < num_markers_per; j++)
			{
				mean_marker_fst = mean_marker_fst + chrom_data_combined[i].marker_fst_values[j];
				d_count++;
			}
		}
		mean_marker_fst = mean_marker_fst / d_count;
		output_line.push_back(mean_marker_fst);

		// Calc peak info
		double n_sig_peaks = 0;
		double n_near_qtl0 = 0;
		double n_near_qtl1 = 0;
		double n_near_qtlp = 0;
		for (i = 0; i < num_chro; i++)
		{
			for (j = 0; j < chrom_data_combined[i].peakdata.numberofpeaks; j++)
			{
				if (chrom_data_combined[i].peakdata.sig99CI[j] == true)
				{
					n_sig_peaks++;
					if (chrom_data_combined[i].peakdata.distancetoQTL_0[j] >= 0 && chrom_data_combined[i].peakdata.distancetoQTL_0[j] <= 25)
					{
						n_near_qtl0++;
					}
					if (chrom_data_combined[i].peakdata.distancetoQTL_1[j] >= 0 && chrom_data_combined[i].peakdata.distancetoQTL_1[j] <= 25)
					{
						n_near_qtl1++;
					}
					if (chrom_data_combined[i].peakdata.distancetoQTL_P[j] >= 0 && chrom_data_combined[i].peakdata.distancetoQTL_P[j] <= 25)
					{
						n_near_qtlp++;
					}
				}
			}
		}
		output_line.push_back(n_sig_peaks);
		output_line.push_back(n_near_qtl0);
		output_line.push_back(n_near_qtl1);
		output_line.push_back(n_near_qtlp);

		// Calc QTL fsts
		double mean_qtl0_fst = 0;
		double mean_qtl1_fst = 0;
		double mean_qtlp_fst = 0;
		double qtl1_count = 0;
		double qtl0_count = 0;
		double qtlp_count = 0;

		for (j = 0; j < num_chro; j++)
		{
			for (k = 0; k < num_qtl0_per; k++)
			{
				qtl0_count++;
				mean_qtl0_fst = mean_qtl0_fst + chrom_data_combined[j].qtl_fst_values_trt0[k];
			}
			for (k = 0; k < num_qtl1_per; k++)
			{
				qtl1_count++;
				mean_qtl1_fst = mean_qtl1_fst + chrom_data_combined[j].qtl_fst_values_trt1[k];
			}
			for (k = 0; k < num_qtlp_per; k++)
			{
				qtlp_count++;
				mean_qtlp_fst = mean_qtlp_fst + chrom_data_combined[j].qtl_fst_values_pleio[k];
			}
		}

		if (qtl0_count > 0)
			mean_qtl0_fst = mean_qtl0_fst / qtl0_count;
		else
			mean_qtl0_fst = -1;

		if (qtl1_count > 0)
			mean_qtl1_fst = mean_qtl1_fst / qtl1_count;
		else
			mean_qtl1_fst = -1;

		if (qtlp_count > 0)
			mean_qtlp_fst = mean_qtlp_fst / qtlp_count;
		else
			mean_qtlp_fst = -1;

		output_line.push_back(mean_qtl0_fst);
		output_line.push_back(mean_qtl1_fst);
		output_line.push_back(mean_qtlp_fst);

		// Calc peak info for FSTprime
		n_sig_peaks = 0;
		n_near_qtl0 = 0;
		n_near_qtl1 = 0;
		n_near_qtlp = 0;
		for (i = 0; i < num_chro; i++)
		{
			for (j = 0; j < chrom_data_combined[i].peakdata.numberofpeaks; j++)
			{
				if (chrom_data_combined[i].peakdata.sig01FSTprime[j] == true)
				{
					n_sig_peaks++;
					if (chrom_data_combined[i].peakdata.distancetoQTL_0[j] >= 0 && chrom_data_combined[i].peakdata.distancetoQTL_0[j] <= 25)
					{
						n_near_qtl0++;
					}
					if (chrom_data_combined[i].peakdata.distancetoQTL_1[j] >= 0 && chrom_data_combined[i].peakdata.distancetoQTL_1[j] <= 25)
					{
						n_near_qtl1++;
					}
					if (chrom_data_combined[i].peakdata.distancetoQTL_P[j] >= 0 && chrom_data_combined[i].peakdata.distancetoQTL_P[j] <= 25)
					{
						n_near_qtlp++;
					}
				}
			}
		}
		output_line.push_back(n_sig_peaks);
		output_line.push_back(n_near_qtl0);
		output_line.push_back(n_near_qtl1);
		output_line.push_back(n_near_qtlp);


		std::ofstream ootfile;
		std::string temp_fn;
		char temp_fn_ch[256];
		sss.str("");
		sss << outfile_name << "__oneline_summary_gen_" << generations << ".csv";
		temp_fn = sss.str();
		convert_string_to_char(temp_fn, temp_fn_ch);
		ootfile.open(temp_fn_ch);

		for (i = 0; i < output_line.size(); i++)
		{
			if (i != 0)
				ootfile << ",";
			ootfile << output_line[i];
		}

		ootfile.close();
			
		sss.str("");
		sss << outfile_name << "__oneline_summary_header.csv";
		temp_fn = sss.str();
		convert_string_to_char(temp_fn, temp_fn_ch);
		ootfile.open(temp_fn_ch);

		ootfile << out_header;

		ootfile.close();

		for (i = 0; i < Number_of_populations; i++)
			population[i].deinitialize_population();

		return 1;
	}
};
