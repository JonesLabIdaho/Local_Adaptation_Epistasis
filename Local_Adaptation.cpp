#include <iostream>
#include "MTwisterFunctions.h"
#include "simulation_engine.h"
#include <sstream>

int main(int argc, char* argv[])
{
	int flag;
	parameter_value_set parameter_values;
	flag = parse_command_line_arguments(argc, argv, parameter_values);
	if (flag == 1)
	{
		// output help file and exit
		return 0;
	}

	std::stringstream ss;
	int rep, number_of_replicates;
	number_of_replicates = parameter_values.p_reps;

	std::string outfile_name, long_outfile_name;
	outfile_name = parameter_values.file_name;

	meta_population *mypops = new meta_population[1];

	for (rep = 0; rep < number_of_replicates; rep++)
	{
		ss.str("");
		ss << "rep_" << rep + 1 << "_" << outfile_name;
		long_outfile_name = ss.str();

		mypops[0].initialize_metapop(parameter_values.p_migration_rate);
		mypops[0].set_parameter_values(parameter_values);
		mypops[0].save_parameter_values(long_outfile_name);
		mypops[0].run_simulation(long_outfile_name);
		mypops[0].deinitialize();
	}
	
	delete[] mypops;

	return 0;
}


