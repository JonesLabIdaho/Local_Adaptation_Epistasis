
#pragma once
#include <algorithm>
#include <sstream>
#include <ctime>
#include "MTwisterFunctions.h"
#include <vector>
#include "ChiSquare.h"
using namespace std;

class Chromosome
{
public:
	int* MarkerLoci;
	double* QTLeffect;
};

class QTLloc
{
public:
	int* QTLlocation;
};

class Gamete
{
public:
	Chromosome* GameteChromosome;
};

class Individual
{
public:
	bool IsFemale;
	bool IsAlive;
	Chromosome* PaternalChromosome;
	Chromosome* MaternalChromosome;
	double Phenotype;
	double BreedingValue;
	double matingsuccess;
	int NumberOfOffspring;
	vector<int> OffspringIndex;
	int FatherIndex;
	int MotherIndex;
};

class MarkerAlleleFrequencies
{
public:
	double* m_allelefrequency;
};

class QTLAlleleFrequencies
{
public:
	int m_numberalleles;
	vector<double> m_qtlallelefreq;
	vector<double> m_qtlallelename;
};

class ChromosomeData
{
public:
	double* m_Fst;
	double* m_FstSmooth;
	double* m_ChiSqrPval;
	bool* m_significant;
	bool* m_smoothedFstpeak;
	MarkerAlleleFrequencies* m_markerallelefrequenciesadults;
	QTLAlleleFrequencies* m_qtlallelefrequenciesadults;
	MarkerAlleleFrequencies* m_markerallelefrequenciesprogeny;
	QTLAlleleFrequencies* m_qtlallelefrequenciesprogeny;
};

class PeakData
{
public:
	int numberofpeaks;
	vector<int> peaklocation;
	vector<double> peaksmoothedFst;
	vector<double> highestFstonplateau;
	vector<int> highestFstIndex;
	vector<bool> sigFDR;
	vector<bool> sig80CI;
	vector<bool> sig90CI;
	vector<bool> sig95CI;
	vector<bool> sig98CI;
	vector<bool> sig99CI;
	vector<bool> sig95BS;
	vector<bool> sig99BS;
	vector<int> nearestQTL;
	vector<int> distancetoQTL;

};

class Population
{
public:
	int PopulationSize;
	int CarryingCapacity;
	int NumberOfProgeny;
	int MaxNumberOffspring;
	int MaximumFecundity;
	int NumberChromosomes;
	int NumberMarkerLociPerChromosome;
	int NumberQTLsPerChromosome;
	Individual* Adult;
	Individual* Progeny;
	double StartingP;
	double Trait0MutationalVariance, Trait0MutationalStdDev;
	double Trait0optimum, Omega;
	QTLloc* QTLlocus;  // need an array that covers all chromosomes
	int MaxMateEncounters;
	double AttractivenessThresh;
	double MateProbBelowThresh;
	double MateProbAboveThresh;
	double MutationRatePerMarker;
	double MutationRatePerQTL;
	double ExpectedQTLMutationsPerChromosome;
	double ExpectedMarkerMutationsPerChromosome;
	double sqrtEMMPC, sqrtEQTLMPC;
	double EnvironmentalVariance, EnvStdDev;
	int MaxNumAllelesPerMarker;
	double ExpRecombPerChromosome;
	int SampleSizeAdults;
	int SampleSizeProgeny;
	int ActualAdultSampleSize;
	int ActualProgenySampleSize;
	int *AdultSample;
	int *ProgenySample;
	double *tempallelefreq;
	double SmoothedCritValue99, SmoothedCritValue98, SmoothedCritValue95, SmoothedCritValue90, SmoothedCritValue80;
	double bs95cutoff, bs99cutoff;
	double FDRcriticalvalue;
	double FDRalpha;
	int peakplateauwidth;  // How many marker loci do we look in each direction for significance?
	double weightingvariance;
	int bootstrapreps;
	int preliminarygenerations;
	int experimentalgenerations;

	double LinearPreferenceSlope;
	double ProbabilityOfMatingWithZero;
	double GaussianPreferenceMean;
	double GaussianPreferenceVariance;
	int MaximumMatingEncounters;
	double StartingQTLAllelicEffectStdDev;

	enum PreferenceFunction{pf_random,pf_linear,pf_threshold,pf_gaussian};
	PreferenceFunction pf_type;
	PreferenceFunction pf_type_initial;
	PreferenceFunction pf_type_experimental;

	// Declare the variables that need to be calculated here

	double MeanPhenotype, MeanBreedingValue;
	ChromosomeData* chromosomedata;
	PeakData* peakdata;

public: // Population Functions
	void InitializeDefaultParameterValues() // do this before initialize or program will crash
	{
		preliminarygenerations = 0;
		experimentalgenerations = 2;
		CarryingCapacity = 5000;
		MaximumFecundity = 4;
		NumberChromosomes = 4;
		NumberMarkerLociPerChromosome = 1000;
		NumberQTLsPerChromosome = 2;
		Trait0MutationalVariance = 1;
		MutationRatePerMarker = 0.000;
		MutationRatePerQTL = 0.000;
		EnvironmentalVariance = 0.0;
		Trait0optimum = 0;
		Omega = 0;
		MaxNumAllelesPerMarker = 4;
		ExpRecombPerChromosome = 0.5;
		StartingQTLAllelicEffectStdDev = 0.5;
		FDRalpha = 0.05;
		peakplateauwidth = 5;
		pf_type_initial = pf_random;
		pf_type_experimental = pf_gaussian;
		weightingvariance = 225;
		bootstrapreps = 10000;

		//Set default type of preference function

		pf_type = pf_gaussian;
		MaximumMatingEncounters = 20;

		// Mating System Parameters (Threshold)

		AttractivenessThresh = 2;  
		MateProbBelowThresh = 0.20;  // Probability of mating with a male below threshold
		MateProbAboveThresh = 0.80;  // Prob. of mating with a male above threshold

		// Mating System Parameters (Linear Preference)

		LinearPreferenceSlope = 0.5;
		ProbabilityOfMatingWithZero = 0.0;

		// Mating System Parameters (Gaussian Preference)

		GaussianPreferenceMean = 2.0;
		GaussianPreferenceVariance = 0.25;

		// Sampling Parameters

		SampleSizeAdults = 2500;
		SampleSizeProgeny = SampleSizeAdults;  // For now, we sample one progeny per female.

	}

	void Initialize() 
	{
		int iii, jjj, kkk;
		double dubrnd;

		// Initialize random number generator

		sgenrand(static_cast<unsigned long>(time(NULL)));

		// Calculate some variables based on the parameters
		
		MaxNumberOffspring = MaximumFecundity*CarryingCapacity;
		double tempdub;
		Trait0MutationalStdDev = sqrt(Trait0MutationalVariance);
		PopulationSize = CarryingCapacity;
		tempdub = NumberMarkerLociPerChromosome;
		ExpectedMarkerMutationsPerChromosome = MutationRatePerMarker*tempdub;
		tempdub = NumberQTLsPerChromosome;
		ExpectedQTLMutationsPerChromosome = MutationRatePerQTL*tempdub;
		sqrtEMMPC = sqrt(ExpectedMarkerMutationsPerChromosome);
		sqrtEQTLMPC = sqrt(ExpectedQTLMutationsPerChromosome);
		EnvStdDev = sqrt(EnvironmentalVariance);

		AdultSample = new int[SampleSizeAdults];
		ProgenySample = new int[SampleSizeProgeny];

		// class level variables for calculating variables of interest
		tempallelefreq = new double[MaxNumAllelesPerMarker];
		

		// Allocate all the memory that will be needed for the population
		Adult = new Individual[CarryingCapacity];
		Progeny = new Individual[MaxNumberOffspring];
		QTLlocus = new QTLloc[NumberChromosomes];
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			QTLlocus[iii].QTLlocation = new int[NumberQTLsPerChromosome];
		}

		for (iii = 0; iii < CarryingCapacity; iii++)
		{
			Adult[iii].MaternalChromosome = new Chromosome[NumberChromosomes];
			Adult[iii].PaternalChromosome = new Chromosome[NumberChromosomes];
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				Adult[iii].MaternalChromosome[jjj].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				Adult[iii].MaternalChromosome[jjj].QTLeffect = new double[NumberQTLsPerChromosome];

				Adult[iii].PaternalChromosome[jjj].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				Adult[iii].PaternalChromosome[jjj].QTLeffect = new double[NumberQTLsPerChromosome];
			}
		}

		for (iii = 0; iii < MaxNumberOffspring; iii++)
		{
			Progeny[iii].MaternalChromosome = new Chromosome[NumberChromosomes];
			Progeny[iii].PaternalChromosome = new Chromosome[NumberChromosomes];
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				Progeny[iii].MaternalChromosome[jjj].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				Progeny[iii].MaternalChromosome[jjj].QTLeffect = new double[NumberQTLsPerChromosome];

				Progeny[iii].PaternalChromosome[jjj].MarkerLoci = new int[NumberMarkerLociPerChromosome];
				Progeny[iii].PaternalChromosome[jjj].QTLeffect = new double[NumberQTLsPerChromosome];
			}
		}

		chromosomedata = new ChromosomeData[NumberChromosomes];
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			chromosomedata[iii].m_Fst = new double[NumberMarkerLociPerChromosome];
			chromosomedata[iii].m_FstSmooth = new double[NumberMarkerLociPerChromosome];
			chromosomedata[iii].m_ChiSqrPval = new double[NumberMarkerLociPerChromosome];
			chromosomedata[iii].m_significant = new bool[NumberMarkerLociPerChromosome];
			chromosomedata[iii].m_smoothedFstpeak = new bool[NumberMarkerLociPerChromosome];
			chromosomedata[iii].m_markerallelefrequenciesadults = new MarkerAlleleFrequencies[NumberMarkerLociPerChromosome];
			chromosomedata[iii].m_qtlallelefrequenciesadults = new QTLAlleleFrequencies[NumberQTLsPerChromosome];
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
				chromosomedata[iii].m_markerallelefrequenciesadults[jjj].m_allelefrequency = new double[MaxNumAllelesPerMarker];
			chromosomedata[iii].m_markerallelefrequenciesprogeny = new MarkerAlleleFrequencies[NumberMarkerLociPerChromosome];
			chromosomedata[iii].m_qtlallelefrequenciesprogeny = new QTLAlleleFrequencies[NumberQTLsPerChromosome];
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
				chromosomedata[iii].m_markerallelefrequenciesprogeny[jjj].m_allelefrequency = new double[MaxNumAllelesPerMarker];
		}

		peakdata = new PeakData[NumberChromosomes];

		// Set the initial conditions for the population

		double* tempQTLallele = new double[MaxNumAllelesPerMarker];
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempQTLallele[iii] = randnorm(0,StartingQTLAllelicEffectStdDev);

		for (iii = 0; iii < CarryingCapacity; iii++) 
		{
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberMarkerLociPerChromosome; kkk++)
				{
					Adult[iii].MaternalChromosome[jjj].MarkerLoci[kkk] = iii%MaxNumAllelesPerMarker;
					Adult[iii].PaternalChromosome[jjj].MarkerLoci[kkk] = iii%MaxNumAllelesPerMarker;
				} // kkk loop

				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
				{
					Adult[iii].MaternalChromosome[jjj].QTLeffect[kkk] = tempQTLallele[iii%MaxNumAllelesPerMarker];
					Adult[iii].PaternalChromosome[jjj].QTLeffect[kkk] = tempQTLallele[iii%MaxNumAllelesPerMarker];
				} // kkk loop

			} // jjj loop

			dubrnd = genrand();
			if (dubrnd < 0.5)
				Adult[iii].IsFemale = true;
			else
				Adult[iii].IsFemale = false;

			Adult[iii].IsAlive = true;

			Adult[iii].BreedingValue = 0;
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
				{
					Adult[iii].BreedingValue = Adult[iii].BreedingValue + Adult[iii].MaternalChromosome[jjj].QTLeffect[kkk] + Adult[iii].PaternalChromosome[jjj].QTLeffect[kkk];
				}
			}
			Adult[iii].Phenotype = Adult[iii].BreedingValue + randnorm(0,EnvStdDev);
		} // iii loop

		// Initialize Progeny with zeros for everything

		for (iii = 0; iii < MaxNumberOffspring; iii++)
		{
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberMarkerLociPerChromosome; kkk++)
				{
					Progeny[iii].MaternalChromosome[jjj].MarkerLoci[kkk] = 0;
					Progeny[iii].PaternalChromosome[jjj].MarkerLoci[kkk] = 0;
				}

				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
				{
					Progeny[iii].MaternalChromosome[jjj].QTLeffect[kkk] = 0;
					Progeny[iii].PaternalChromosome[jjj].QTLeffect[kkk] = 0;
				}

			} // jjj
			Progeny[iii].IsAlive = false;
			Progeny[iii].IsFemale = false;
		}

		// Set the location of the QTLs relative to the marker loci
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for (jjj = 0; jjj < NumberQTLsPerChromosome; jjj++)
			{
				QTLlocus[iii].QTLlocation[jjj] = randnum(NumberMarkerLociPerChromosome);
			}
		}

		PopulationSize = CarryingCapacity;

		delete[] tempQTLallele;

	}	// End of Initialize Function

	void Deinitialize()
	{
		int iii, jjj;
		// Free up all of the memory declared during the initialization
		for (iii = 0; iii < CarryingCapacity; iii++)
		{
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				delete[] Adult[iii].MaternalChromosome[jjj].MarkerLoci;
				delete[] Adult[iii].MaternalChromosome[jjj].QTLeffect;

				delete[] Adult[iii].PaternalChromosome[jjj].MarkerLoci;
				delete[] Adult[iii].PaternalChromosome[jjj].QTLeffect;
			}
			delete[] Adult[iii].MaternalChromosome;
			delete[] Adult[iii].PaternalChromosome;
		}

		for (iii = 0; iii < MaxNumberOffspring; iii++)
		{
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				delete[] Progeny[iii].MaternalChromosome[jjj].MarkerLoci;
				delete[] Progeny[iii].MaternalChromosome[jjj].QTLeffect;

				delete[] Progeny[iii].PaternalChromosome[jjj].MarkerLoci;
				delete[] Progeny[iii].PaternalChromosome[jjj].QTLeffect;
			}
			delete[] Progeny[iii].MaternalChromosome;
			delete[] Progeny[iii].PaternalChromosome;
		}

		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			delete[] QTLlocus[iii].QTLlocation;
		}
		delete[] QTLlocus;
		delete[] Adult;
		delete[] Progeny;

		delete[] AdultSample;
		delete[] ProgenySample;

				
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			delete[] chromosomedata[iii].m_Fst;
			delete[] chromosomedata[iii].m_FstSmooth;
			delete[] chromosomedata[iii].m_ChiSqrPval;
			delete[] chromosomedata[iii].m_significant;
			delete[] chromosomedata[iii].m_smoothedFstpeak;
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
				delete[] chromosomedata[iii].m_markerallelefrequenciesadults[jjj].m_allelefrequency;
			delete[] chromosomedata[iii].m_markerallelefrequenciesadults;
			delete[] chromosomedata[iii].m_qtlallelefrequenciesadults;
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
				delete[] chromosomedata[iii].m_markerallelefrequenciesprogeny[jjj].m_allelefrequency;
			delete[] chromosomedata[iii].m_markerallelefrequenciesprogeny;
			delete[] chromosomedata[iii].m_qtlallelefrequenciesprogeny;
		}
		delete[] chromosomedata;

		delete[] peakdata;

		delete[] tempallelefreq;


	} // End of deinitialize function

	void OutputAdultsToFile(int gen)
	{
		int iii, jjj, kkk;
		string a_filename;
		string tempstring;
		stringstream tsstr;
		tsstr << "Adults_Gen_";
		tsstr << gen;
		tsstr << ".xls";
		a_filename = tsstr.str();
		ofstream a_file;
		a_file.open(a_filename);
		
		a_file << "Adult\tSex\tAlive\tPheno\tGeno\t";
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
				a_file << "C" << iii << "M" << jjj << "\t";
			for (jjj = 0; jjj < NumberQTLsPerChromosome; jjj++)
				a_file << "C" << iii << "Q" << jjj << "\t";
		}
		a_file << "\n";

		for (iii = 0; iii < PopulationSize; iii++)
		{
			a_file << iii << '\t';
			if (Adult[iii].IsFemale)
				a_file << "fem\t";
			else
				a_file << "mal\t";
			if (Adult[iii].IsAlive)
				a_file << "live\t";
			else
				a_file << "dead\t";
			a_file << Adult[iii].Phenotype << "\t";
			a_file << Adult[iii].BreedingValue << "\t";

			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberMarkerLociPerChromosome; kkk++)
					a_file << Adult[iii].MaternalChromosome[jjj].MarkerLoci[kkk] << "x" << Adult[iii].PaternalChromosome[jjj].MarkerLoci[kkk] << "\t";
				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
					a_file << Adult[iii].MaternalChromosome[jjj].QTLeffect[kkk] << "x" << Adult[iii].PaternalChromosome[jjj].QTLeffect[kkk] << "\t";
			} // jjj
			a_file << "\n";
		} // iii
	} // Output Adults End

	void OutputProgenyToFile(int gen)
	{
		int iii, jjj, kkk;
		string a_filename;
		string tempstring;
		stringstream tsstr;
		tsstr << "Progeny_Gen_";
		tsstr << gen;
		tsstr << ".xls";
		a_filename = tsstr.str();
		ofstream a_file;
		a_file.open(a_filename);
		
		a_file << "Prog\tSex\tAlive\tPheno\tGeno\t";
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
				a_file << "C" << iii << "M" << jjj << "\t";
			for (jjj = 0; jjj < NumberQTLsPerChromosome; jjj++)
				a_file << "C" << iii << "Q" << jjj << "\t";
		}
		a_file << "\n";

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			a_file << iii << '\t';
			if (Progeny[iii].IsFemale)
				a_file << "fem\t";
			else
				a_file << "mal\t";
			if (Progeny[iii].IsAlive)
				a_file << "live\t";
			else
				a_file << "dead\t";
			a_file << Progeny[iii].Phenotype << "\t";
			a_file << Progeny[iii].BreedingValue << "\t";

			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberMarkerLociPerChromosome; kkk++)
					a_file << Progeny[iii].MaternalChromosome[jjj].MarkerLoci[kkk] << "x" << Progeny[iii].PaternalChromosome[jjj].MarkerLoci[kkk] << "\t";
				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
					a_file << Progeny[iii].MaternalChromosome[jjj].QTLeffect[kkk] << "x" << Progeny[iii].PaternalChromosome[jjj].QTLeffect[kkk] << "\t";
			} // jjj
			a_file << "\n";
		} // iii
	} // end of output progeny

	void ProduceRecombinedChromosome(Chromosome &RecombinedChromosome, Individual &Parent, int WhichChromosome, double ExpectedRecombinationEvents) // Chromosome Result, Adult Parent, Chromosome Number, Max of 8 events
	{
		int PRi, PRj;
		int NumberRecombinationEvents = 0;
		int SegmentStart[22], SegmentEnd[22];
		int BreakPoint[20];

		if (ExpectedRecombinationEvents < 6)
			NumberRecombinationEvents = poissonrand(ExpectedRecombinationEvents);
		if (ExpectedRecombinationEvents >= 6)
			NumberRecombinationEvents = positiveroundnorm(ExpectedRecombinationEvents,sqrt(ExpectedRecombinationEvents));

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
			sort(begin(BreakPoint),end(BreakPoint));

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
				SegmentStart[PRi+1] = BreakPoint[PRi];
				if (SegmentMaternal[PRi])
					SegmentMaternal[PRi+1] = false;
				else
					SegmentMaternal[PRi+1] = true;
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

			// Next do it for the QTLs
			for (PRi = 0; PRi < NumberQTLsPerChromosome; PRi++)
			{
				for (PRj = 0; PRj < NumberSegments; PRj++)
				{
					if (QTLlocus[WhichChromosome].QTLlocation[PRi] >= SegmentStart[PRj] && QTLlocus[WhichChromosome].QTLlocation[PRi] < SegmentEnd[PRj])
					{
						if(SegmentMaternal[PRj])
							RecombinedChromosome.QTLeffect[PRi] = Parent.MaternalChromosome[WhichChromosome].QTLeffect[PRi];
						else
							RecombinedChromosome.QTLeffect[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect[PRi];
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
			}
			else
			{
				for (PRi = 0; PRi < NumberMarkerLociPerChromosome; PRi++)
					RecombinedChromosome.MarkerLoci[PRi] = Parent.PaternalChromosome[WhichChromosome].MarkerLoci[PRi];
				for (PRi = 0; PRi < NumberQTLsPerChromosome; PRi++)
					RecombinedChromosome.QTLeffect[PRi] = Parent.PaternalChromosome[WhichChromosome].QTLeffect[PRi];
			}
		} // else
	}  // End of ProduceRecombinedChromosome

	void PolygynousMating()
	{
		int iii, mmm, nnn;
		int MaleID, FemaleID;
		double drand;

		int ProgenyCounter;
		int NumberMales;
		bool MateFound;
		double dNumberMales, dMaleMeanPhenotype;

		int *LivingMaleIDList = new int[PopulationSize];

		//int NumberFemales;
		//double dNumberFemales;
		//bool SampleCurrentFemale;
		//double dFemaleSampleUnfilled;
		//double dNumberUnsampledFemalesLeft;
		//double dSampleProb;
		//int iCurrentSampleSize;

		// Count up the number of males in the population
		// and make a list of them
		// Also need the mean male phenotype
		NumberMales = 0;
		dMaleMeanPhenotype = 0;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			if (Adult[iii].IsFemale == false)
			{
				LivingMaleIDList[NumberMales] = iii;
				dMaleMeanPhenotype = dMaleMeanPhenotype + Adult[iii].Phenotype;
				Adult[iii].matingsuccess = 0;
				NumberMales++;
			} // end of if
			
			// Regardless of whether the individual is male or female,
			// we have to clear the vector of their offspring

			Adult[iii].OffspringIndex.clear();
			Adult[iii].NumberOfOffspring = 0;

		} // end of iii
		dNumberMales = NumberMales;
		dMaleMeanPhenotype = dMaleMeanPhenotype/dNumberMales;

		//NumberFemales = PopulationSize - NumberMales;
		//dNumberFemales = NumberFemales;
		//dNumberUnsampledFemalesLeft = dNumberFemales;
		//dFemaleSampleUnfilled = SampleSizeAdults;
		//iCurrentSampleSize = 0;

		// Finish this sampling of females and offspring routine
		//
		//if (SampleSizeAdults > NumberFemales)
		//	SampleSizeAdults = NumberFemales;
		//SampleSizeProgeny = SampleSizeAdults;

		double dtempmatingprob;
		int EncounterCounter;
		ProgenyCounter = 0;
		for (iii = 0; iii < PopulationSize; iii++) 
		{
			// Check to see if the individual is female 
			
			if (Adult[iii].IsFemale && NumberMales > 0)
			{
				FemaleID = iii;

				// Threshold mate choice
				if(pf_type == pf_threshold)
				{
					// Find a mate for the female
					MateFound = false;
					EncounterCounter = 0;
					while (!MateFound && EncounterCounter < MaximumMatingEncounters)
					{
						MaleID = LivingMaleIDList[randnum(NumberMales)];
						drand = genrand();
						if (Adult[MaleID].Phenotype >= (AttractivenessThresh) && drand < MateProbAboveThresh)
							MateFound = true;
						if (Adult[MaleID].Phenotype < (AttractivenessThresh) && drand < MateProbBelowThresh)
							MateFound = true;
						EncounterCounter++;
					} // end of while
					
					if (MateFound)
						Adult[MaleID].matingsuccess++;
				}

				// Gaussian Mate Choice

				if (pf_type == pf_gaussian)
				{
					MateFound = false;
					EncounterCounter = 0;
					while (!MateFound && EncounterCounter < MaximumMatingEncounters)
					{
						MaleID = LivingMaleIDList[randnum(NumberMales)];
						dtempmatingprob = exp(-0.5*(Adult[MaleID].Phenotype-GaussianPreferenceMean)*(Adult[MaleID].Phenotype-GaussianPreferenceMean)
												/GaussianPreferenceVariance);
						if (genrand()<dtempmatingprob)
							MateFound = true;
						EncounterCounter++;
					}
					if (MateFound)
						Adult[MaleID].matingsuccess++;
				}

				// Linear Mate Choice

				if (pf_type == pf_linear)
				{
					MateFound = false;
					EncounterCounter = 0;
					while(!MateFound && EncounterCounter < MaximumMatingEncounters)
					{
						MaleID = LivingMaleIDList[randnum(NumberMales)];
						dtempmatingprob = ProbabilityOfMatingWithZero + LinearPreferenceSlope*Adult[MaleID].Phenotype;
						if (genrand()<dtempmatingprob)
							MateFound = true;
						EncounterCounter++;
					}
					if (MateFound)
						Adult[MaleID].matingsuccess++;
				}

				// Random Mating
				if(pf_type == pf_random || !MateFound)
				{
					MateFound = false;
					while (!MateFound)
					{
						MaleID = LivingMaleIDList[randnum(NumberMales)];
						MateFound = true;
					}
					Adult[MaleID].matingsuccess++;
				}

				// Create Progeny with this male and female
				for (mmm = 0; mmm < MaximumFecundity; mmm++)
				{
					// The mother and father each pass a chromosome to the offspring
					for (nnn = 0; nnn < NumberChromosomes; nnn++)
					{
						ProduceRecombinedChromosome(Progeny[ProgenyCounter].MaternalChromosome[nnn],Adult[FemaleID],nnn,ExpRecombPerChromosome); 
						ProduceRecombinedChromosome(Progeny[ProgenyCounter].PaternalChromosome[nnn],Adult[MaleID],nnn,ExpRecombPerChromosome);
					} // end of nnn

					// Keep track of parental identities
					// Also keep track of offspring of parents
					Progeny[ProgenyCounter].FatherIndex = MaleID;
					Progeny[ProgenyCounter].MotherIndex = FemaleID;
					Adult[MaleID].OffspringIndex.push_back(ProgenyCounter);
					Adult[MaleID].NumberOfOffspring++;
					Adult[FemaleID].OffspringIndex.push_back(ProgenyCounter);
					Adult[FemaleID].NumberOfOffspring++;

					ProgenyCounter++;
				} // end of mmm

				// Check to see if this female is sampled

				//dSampleProb = dFemaleSampleUnfilled/dNumberUnsampledFemalesLeft;
				//dsamprand = genrand();
				//if (dsamprand < dSampleProb)
				//{
				//	SampleCurrentFemale = true;
				//	dFemaleSampleUnfilled = dFemaleSampleUnfilled - 1;
				//}
				//dNumberUnsampledFemalesLeft = dNumberUnsampledFemalesLeft - 1;
				//
				//if (SampleCurrentFemale)
				//{
				//	AdultSample[iCurrentSampleSize] = iii;
				//	ProgenySample[iCurrentSampleSize] = ProgenyCounter - 1;
				//	iCurrentSampleSize++;
				//}

			} // end of if (Adult[iii].IsFemale && NumberMales > 0)
		}  // end of iii

		NumberOfProgeny = ProgenyCounter;

		delete[] LivingMaleIDList;
	}  // end of PolygynousMating()

	double CalcMatingDifferential()
	{
		int iii;
		double matingdiff;
		double meanMSmales = 0;
		double meanTraitmales = 0;
		double Nmales = 0;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			if (!Adult[iii].IsFemale)
			{
				meanMSmales = meanMSmales + Adult[iii].matingsuccess;
				meanTraitmales = meanTraitmales + Adult[iii].Phenotype;
				Nmales++;
			}
		}
		meanMSmales = meanMSmales/Nmales;
		meanTraitmales = meanTraitmales/Nmales;

		matingdiff = 0;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			if (!Adult[iii].IsFemale)
			{
				matingdiff = matingdiff + (Adult[iii].Phenotype - meanTraitmales)*(Adult[iii].matingsuccess/meanMSmales-1);
			}
		}
		matingdiff = matingdiff/Nmales;
		return matingdiff;
	}

	void Mutation()
	{
		int iii, jjj, kkk, NumberMutations;
		int whichlocus;

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				// First perform mutations at the marker loci
				// Maternal chromosome
				if (ExpectedMarkerMutationsPerChromosome < 6)
					NumberMutations = poissonrand(ExpectedMarkerMutationsPerChromosome);
				else
					NumberMutations = positiveroundnorm(ExpectedMarkerMutationsPerChromosome,sqrtEMMPC);

				for (kkk = 0; kkk < NumberMutations; kkk++)
				{
					whichlocus = randnum(NumberMarkerLociPerChromosome);
					Progeny[iii].MaternalChromosome[jjj].MarkerLoci[whichlocus] = randnum(MaxNumAllelesPerMarker);
				}

				// Paternal chromosome
				if (ExpectedMarkerMutationsPerChromosome < 6)
					NumberMutations = poissonrand(ExpectedMarkerMutationsPerChromosome);
				else
					NumberMutations = positiveroundnorm(ExpectedMarkerMutationsPerChromosome,sqrtEMMPC);

				for (kkk = 0; kkk < NumberMutations; kkk++)
				{
					whichlocus = randnum(NumberMarkerLociPerChromosome);
					Progeny[iii].PaternalChromosome[jjj].MarkerLoci[whichlocus] = randnum(MaxNumAllelesPerMarker);
				}

				// Next perform mutations at the QTL loci
				// Maternal chromosome
				if (ExpectedQTLMutationsPerChromosome < 6)
					NumberMutations = poissonrand(ExpectedQTLMutationsPerChromosome);
				else
					NumberMutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome,sqrtEQTLMPC);

				for (kkk = 0; kkk < NumberMutations; kkk++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome);
					Progeny[iii].MaternalChromosome[jjj].QTLeffect[whichlocus] = Progeny[iii].MaternalChromosome[jjj].QTLeffect[whichlocus] + randnorm(0,Trait0MutationalStdDev);
				}

				// Paternal chromosome
				if (ExpectedQTLMutationsPerChromosome < 6)
					NumberMutations = poissonrand(ExpectedQTLMutationsPerChromosome);
				else
					NumberMutations = positiveroundnorm(ExpectedQTLMutationsPerChromosome,sqrtEQTLMPC);

				for (kkk = 0; kkk < NumberMutations; kkk++)
				{
					whichlocus = randnum(NumberQTLsPerChromosome);
					Progeny[iii].MaternalChromosome[jjj].QTLeffect[whichlocus] = Progeny[iii].MaternalChromosome[jjj].QTLeffect[whichlocus] + randnorm(0,Trait0MutationalStdDev);
				}

			} // end of jjj


		} // end of iii

		// Now it's time to figure the progeny's phenotype at the QTL
		// This means that the mutation subroutine must always be run, even when the mutation rates are zero

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			Progeny[iii].BreedingValue = 0;
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
				{
					Progeny[iii].BreedingValue = Progeny[iii].BreedingValue + Progeny[iii].MaternalChromosome[jjj].QTLeffect[kkk] + Progeny[iii].PaternalChromosome[jjj].QTLeffect[kkk];
				} // kkk
			} // jjj
			Progeny[iii].Phenotype = Progeny[iii].BreedingValue + randnorm(0,EnvStdDev);
			Progeny[iii].IsAlive = true;
			if (genrand() < 0.5)
				Progeny[iii].IsFemale = true;
			else
				Progeny[iii].IsFemale = false;
		} // end of iii
		// Now the progeny are all fully created

	} // end of Mutation

	void Selection()
	{
		int iii;
		double dSurvProb;

		// Natural selection only affects the males, because this trait is only expressed in males
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (!Progeny[iii].IsFemale && Omega > 0)
			{
				dSurvProb = exp(-1*(Progeny[iii].Phenotype - Trait0optimum)*(Progeny[iii].Phenotype - Trait0optimum)/(2*Omega));
				if (genrand() < dSurvProb)
					Progeny[iii].IsAlive = true;
				else
					Progeny[iii].IsAlive = false;
			}
			else
			{
				Progeny[iii].IsAlive = true;
			}
		}
	} // end of selection

	void PopulationRegulation()
	{
		int iii, jjj, kkk;
		double dCACarCapUnfilled, dCAProgleft, dCAKeepProb;
		int iCANumAdultsChosen;
		double dCArnd1;

		dCAProgleft = 0;
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive)
				dCAProgleft++;
		}

		dCACarCapUnfilled = CarryingCapacity;

		iCANumAdultsChosen = 0;
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive)
			{
				dCAKeepProb = dCACarCapUnfilled/dCAProgleft;
				dCArnd1 = genrand();
				if (dCArnd1 < dCAKeepProb)
				{
					// Copy Progeny Values into Adult Values
					Adult[iCANumAdultsChosen].IsAlive = true;
					Adult[iCANumAdultsChosen].IsFemale = Progeny[iii].IsFemale;
					Adult[iCANumAdultsChosen].BreedingValue = Progeny[iii].BreedingValue;
					Adult[iCANumAdultsChosen].Phenotype = Progeny[iii].Phenotype;
					for (jjj = 0; jjj < NumberChromosomes; jjj++)
					{
						for (kkk = 0; kkk < NumberMarkerLociPerChromosome; kkk++)
						{
							Adult[iCANumAdultsChosen].MaternalChromosome[jjj].MarkerLoci[kkk] = Progeny[iii].MaternalChromosome[jjj].MarkerLoci[kkk];
							Adult[iCANumAdultsChosen].PaternalChromosome[jjj].MarkerLoci[kkk] = Progeny[iii].PaternalChromosome[jjj].MarkerLoci[kkk];
						}
						for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
						{
							Adult[iCANumAdultsChosen].MaternalChromosome[jjj].QTLeffect[kkk] = Progeny[iii].MaternalChromosome[jjj].QTLeffect[kkk];
							Adult[iCANumAdultsChosen].PaternalChromosome[jjj].QTLeffect[kkk] = Progeny[iii].PaternalChromosome[jjj].QTLeffect[kkk];
						}
					}

					dCACarCapUnfilled = dCACarCapUnfilled - 1;
					iCANumAdultsChosen++;

				} // end of if keep prob
				dCAProgleft = dCAProgleft - 1;
			} // end of if IsAlive
		} // end of iii

	} // end of population regulation

	double ParentOffspringFst(int whichchromosome, int whichmarker, bool includedead)
	{
		double *allelefreq = new double[MaxNumAllelesPerMarker];
		double *progenyallelefreq = new double[MaxNumAllelesPerMarker];
		double *adultallelefreq = new double[MaxNumAllelesPerMarker];
		double countalleles;	
		int iii;
		int maternalallele;
		int paternalallele;

		
		// First calculate allele frequencies in the progeny

		countalleles = 0;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			allelefreq[iii] = 0;

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive || includedead)
			{
				maternalallele = Progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				paternalallele = Progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				if (maternalallele >= 0 && maternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[maternalallele]++;
					countalleles++;
				}
				if (paternalallele >= 0 && paternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[paternalallele]++;
					countalleles++;
				}
			}
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			progenyallelefreq[iii] = allelefreq[iii]/countalleles;

		// Second calculate allele frequencies in the adults

		countalleles = 0;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			allelefreq[iii] = 0;

		for (iii = 0; iii < PopulationSize; iii++)
		{
			if (Adult[iii].IsAlive || includedead)
			{
				maternalallele = Adult[iii].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				paternalallele = Adult[iii].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
			
				if (maternalallele >= 0 && maternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[maternalallele]++;
					countalleles++;
				}
				if (paternalallele >= 0 && paternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[paternalallele]++;
					countalleles++;
				}
			}
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			adultallelefreq[iii] = allelefreq[iii]/countalleles;

		// Finally calculate Fst

		double progenyHexp;
		double adultHexp;
		double lumpedHexp;

		progenyHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			progenyHexp = progenyHexp - progenyallelefreq[iii]*progenyallelefreq[iii];
		adultHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			adultHexp = adultHexp - adultallelefreq[iii]*adultallelefreq[iii];

		// Use the allelefreq array to calculate the average allele frequency across adults and progeny.
		// For now, this is an unweighted average.

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			allelefreq[iii] = (progenyallelefreq[iii]+adultallelefreq[iii])/2;
		lumpedHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			lumpedHexp = lumpedHexp - allelefreq[iii]*allelefreq[iii];

		double Fst;

		if (lumpedHexp > 0)
			Fst = 1 - (progenyHexp + adultHexp)/(2*lumpedHexp);
		else
			Fst = -1;

		// Convert to Fst prime just to see what happens
		//if (Fst > 0)
		//	Fst = Fst*(1 + (progenyHexp+adultHexp)/2)/(1-(progenyHexp+adultHexp)/2);

		// Calculate ln(w11/w22)
		/* ln(w11/w22) is no better than Fst
		if (lumpedHexp > 0)
		{
			double pppp;
			double Pprime;
			double w11w22;
			pppp = 0;
			Pprime = 0;

			for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			{
				if (adultallelefreq[iii] > pppp)
					pppp = adultallelefreq[iii];
			}

			for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			{
				if (progenyallelefreq[iii] > Pprime)
					Pprime = progenyallelefreq[iii];
			}

			if (pppp < 1)
			{
				w11w22 = ((1-pppp)/pppp)*((2*Pprime-pppp)/(1+pppp-2*Pprime));
				if (w11w22 > 0)
					Markerlnw11w22 = log(w11w22);
				else
					Markerlnw11w22 = 0;
			}
			else
				Markerlnw11w22 = 0;
		}  // if (Htotal > 0)
		else
			Markerlnw11w22 = 0; */

		delete[] allelefreq;
		delete[] progenyallelefreq;
		delete[] adultallelefreq;

		return Fst;

	}

	double CalcSelectionDifferential()
	{
		int iii;
		double fitnessmean, traitmean, selectiondifferential;
		double dnprog;
		dnprog = NumberOfProgeny;
		fitnessmean = 0;
		traitmean = 0;
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			traitmean = traitmean + Progeny[iii].Phenotype;
			if (Progeny[iii].IsAlive)
				fitnessmean++;
		}
		traitmean = traitmean/dnprog;
		fitnessmean = fitnessmean/dnprog;

		selectiondifferential = 0;
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive)
				selectiondifferential = selectiondifferential + (Progeny[iii].Phenotype - traitmean)*(1.0/fitnessmean - 1.0);
			if (!Progeny[iii].IsAlive)
				selectiondifferential = selectiondifferential + (Progeny[iii].Phenotype - traitmean)*(0.0/fitnessmean - 1.0);
		}
		selectiondifferential=selectiondifferential/dnprog;
		return selectiondifferential;
	}

	double CalcProgenyBreedingValueMean()
	{
		double dnprog = NumberOfProgeny;
		double pbvmean = 0;
		for (int iii = 0; iii < NumberOfProgeny; iii++)
			pbvmean = pbvmean + Progeny[iii].BreedingValue;
		return pbvmean/dnprog;
	}

	double CalcAdultBreedingValueMean()
	{
		double dnadult = PopulationSize;
		double abvmean = 0;
		for (int iii = 0; iii < PopulationSize; iii++)
			abvmean = abvmean + Adult[iii].BreedingValue;
		return abvmean/dnadult;
	}

	double CalcProgenyBVvariance()
	{
		double dnprog = NumberOfProgeny;
		double dprogvar = 0;
		double dprogmean = 0;
		for (int iii = 0; iii < NumberOfProgeny; iii++)
			dprogmean = dprogmean + Progeny[iii].BreedingValue;
		dprogmean = dprogmean /dnprog;
		for (int iii = 0; iii < NumberOfProgeny; iii++)
			dprogvar = dprogvar + (Progeny[iii].BreedingValue-dprogmean)*(Progeny[iii].BreedingValue-dprogmean);
		dprogvar = dprogvar/dnprog;
		return dprogvar;
	}

	double CalcAdultBVvariance()
	{
		double dnadult = PopulationSize;
		double dadultvar = 0;
		double dadultmean = 0;
		for (int iii = 0; iii < PopulationSize; iii++)
			dadultmean = dadultmean + Adult[iii].BreedingValue;
		dadultmean = dadultmean/dnadult;
		for (int iii = 0; iii < PopulationSize; iii++)
			dadultvar = dadultvar + (Adult[iii].BreedingValue - dadultmean)*(Adult[iii].BreedingValue - dadultmean);
		dadultvar = dadultvar/dnadult;
		return dadultvar;
	}

	double PerLocusQTLvarianceOffspring(int whichchromosome, int whichQTL)
	{
		double meanGenotypicEffect = 0;
		double dnprog = NumberOfProgeny;
		for (int iii = 0; iii < NumberOfProgeny; iii++)
			meanGenotypicEffect = meanGenotypicEffect + Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] + Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL]; 
		meanGenotypicEffect = meanGenotypicEffect/dnprog;

		double QTLvariance = 0;
		for (int iii = 0; iii < NumberOfProgeny; iii++)
			QTLvariance = QTLvariance + (Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] + Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] - meanGenotypicEffect)
							* (Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] + Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] - meanGenotypicEffect);
		QTLvariance = QTLvariance/dnprog;
		return QTLvariance;
	}

	double PerLocusQTLvarianceAdults(int whichchromosome, int whichQTL)
	{
		double meanGenotypicEffect = 0;
		double dnadult = PopulationSize;
		for (int iii = 0; iii < PopulationSize; iii++)
			meanGenotypicEffect = meanGenotypicEffect + Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] +
					Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL];
		meanGenotypicEffect = meanGenotypicEffect/dnadult;

		double QTLvariance = 0;
		for (int iii = 0; iii < PopulationSize; iii++)
			QTLvariance = QTLvariance + (Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] + Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] - meanGenotypicEffect)
							* (Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] + Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] - meanGenotypicEffect);
		QTLvariance = QTLvariance/dnadult;
		return QTLvariance;
	}

	double QTLfst(int whichchromosome, int whichQTL)
	{
		// First make a list of alleles at the locus
		vector<double> AlleleList;
		vector<double> AlleleFreqAverage;
		vector<double> AlleleFreqProgeny;
		vector<double> AlleleFreqAdults;
		int iii, jjj;
		int NumAllelesInList = 0;
		bool Onlist;

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
		}

		for (iii = 0; iii < PopulationSize; iii++)
		{
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}

			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
		}

		// Now we should have a list of all the alleles and their number of occurrences

		// Figure out progeny allele frequencies

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					AlleleFreqProgeny[jjj]++;
				if (Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					AlleleFreqProgeny[jjj]++;
			}
		}

		double dNprog = NumberOfProgeny;
		for (jjj = 0; jjj < NumAllelesInList; jjj++)
			AlleleFreqProgeny[jjj] = AlleleFreqProgeny[jjj]/(dNprog*2);

		// Figure out adult allele frequencies

		for (iii = 0; iii < PopulationSize; iii++)
		{
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if (Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					AlleleFreqAdults[jjj]++;
				if (Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					AlleleFreqAdults[jjj]++;
			}
		}

		double dNadults = PopulationSize;
		for (jjj = 0; jjj < NumAllelesInList; jjj++)
			AlleleFreqAdults[jjj] = AlleleFreqAdults[jjj]/(dNadults*2);

		// Figure out average allele frequencies

		for (jjj = 0; jjj < NumAllelesInList; jjj++)
			AlleleFreqAverage[jjj] = (AlleleFreqProgeny[jjj] + AlleleFreqAdults[jjj])/2;

		double dFst, Hprogeny, Hadults, Htotal;

		Hprogeny = 1;
		Hadults = 1;
		Htotal = 1;
		for (jjj = 0; jjj < NumAllelesInList; jjj++)
		{
			Hprogeny = Hprogeny - (AlleleFreqProgeny[jjj]*AlleleFreqProgeny[jjj]);
			Hadults = Hadults - (AlleleFreqAdults[jjj]*AlleleFreqAdults[jjj]);
			Htotal = Htotal - (AlleleFreqAverage[jjj]*AlleleFreqAverage[jjj]);
		}

		if (Htotal > 0)
			dFst = 1 - (Hprogeny+Hadults)/(2*Htotal);
		else
			dFst = -1;
		
		// Calculate ln(w11/w22)
		/*  This turns out not to be especially useful
		if (Htotal > 0)
		{
			double pppp;
			double Pprime;
			double w11w22;
			pppp = 0;
			Pprime = 0;

			for (iii = 0; iii < NumAllelesInList; iii++)
			{
				if (AlleleFreqAdults[iii] > pppp)
					pppp = AlleleFreqAdults[iii];
			}

			for (iii = 0; iii < NumAllelesInList; iii++)
			{
				if (AlleleFreqProgeny[iii] > Pprime)
					Pprime = AlleleFreqProgeny[iii];
			}

			if (pppp < 1 && pppp > 0.05)
			{
				w11w22 = ((1-pppp)/pppp)*((2*Pprime-pppp)/(1+pppp-2*Pprime));
				if (w11w22 > 0)
					QTLlnw11w22 = log(w11w22);
				else
					QTLlnw11w22 = 0;
			}
			else
				QTLlnw11w22 = 0;
		}  // if (Htotal > 0)
		else
			QTLlnw11w22 = 0;*/


		return dFst;
	}

	double MaleOffspringFst(int whichchromosome, int whichmarker, bool includedead)
	{
		double *allelefreq = new double[MaxNumAllelesPerMarker];
		double *progenyallelefreq = new double[MaxNumAllelesPerMarker];
		double *adultallelefreq = new double[MaxNumAllelesPerMarker];
		double countalleles;	
		int iii;
		int maternalallele;
		int paternalallele;

		// First calculate allele frequencies in the progeny

		countalleles = 0;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			allelefreq[iii] = 0;

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive || includedead)
			{
				maternalallele = Progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				paternalallele = Progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				if (maternalallele >= 0 && maternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[maternalallele]++;
					countalleles++;
				}
				if (paternalallele >= 0 && paternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[paternalallele]++;
					countalleles++;
				}
			}
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			progenyallelefreq[iii] = allelefreq[iii]/countalleles;

		// Second calculate allele frequencies in the adult males

		countalleles = 0;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			allelefreq[iii] = 0;

		for (iii = 0; iii < PopulationSize; iii++)
		{
			if ((Adult[iii].IsAlive || includedead) && !Adult[iii].IsFemale)
			{
				maternalallele = Adult[iii].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				paternalallele = Adult[iii].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
			
				if (maternalallele >= 0 && maternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[maternalallele]++;
					countalleles++;
				}
				if (paternalallele >= 0 && paternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[paternalallele]++;
					countalleles++;
				}
			}
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			adultallelefreq[iii] = allelefreq[iii]/countalleles;

		// Finally calculate Fst

		double progenyHexp;
		double adultHexp;
		double lumpedHexp;

		progenyHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			progenyHexp = progenyHexp - progenyallelefreq[iii]*progenyallelefreq[iii];
		adultHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			adultHexp = adultHexp - adultallelefreq[iii]*adultallelefreq[iii];

		// Use the allelefreq array to calculate the average allele frequency across adults and progeny.
		// For now, this is an unweighted average.

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			allelefreq[iii] = (progenyallelefreq[iii]+adultallelefreq[iii])/2;
		lumpedHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			lumpedHexp = lumpedHexp - allelefreq[iii]*allelefreq[iii];

		double Fst;

		if (lumpedHexp > 0)
			Fst = 1 - (progenyHexp + adultHexp)/(2*lumpedHexp);
		else
			Fst = -1;

		delete[] allelefreq;
		delete[] progenyallelefreq;
		delete[] adultallelefreq;

		return Fst;

	}

	double FstProgenyPaternalChromosomesMaternalChromosomes(int whichchromosome, int whichmarker, bool includedead)
	{
		double *allelefreq = new double[MaxNumAllelesPerMarker];
		double *paternalallelefreq = new double[MaxNumAllelesPerMarker];
		double *maternalallelefreq = new double[MaxNumAllelesPerMarker];
		double countalleles;	
		int iii;
		int maternalallele;
		int paternalallele;

		// First calculate allele frequencies in the progeny

		countalleles = 0;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			allelefreq[iii] = 0;
			paternalallelefreq[iii] = 0;
			maternalallelefreq[iii] = 0;
		}

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive || includedead)
			{
				maternalallele = Progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				paternalallele = Progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker];
				if (maternalallele >= 0 && maternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[maternalallele]++;
					maternalallelefreq[maternalallele]++;
					countalleles++;
				}
				if (paternalallele >= 0 && paternalallele < MaxNumAllelesPerMarker)
				{
					allelefreq[paternalallele]++;
					paternalallelefreq[paternalallele]++;
					countalleles++;
				}
			}
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			allelefreq[iii] = allelefreq[iii]/countalleles;
			maternalallelefreq[iii] = 2*maternalallelefreq[iii]/countalleles;
			paternalallelefreq[iii] = 2*paternalallelefreq[iii]/countalleles;
		}

		// Finally calculate Fst

		double paternalHexp;
		double maternalHexp;
		double lumpedHexp;

		paternalHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			paternalHexp = paternalHexp - paternalallelefreq[iii]*paternalallelefreq[iii];
		maternalHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			maternalHexp = maternalHexp - maternalallelefreq[iii]*maternalallelefreq[iii];
		lumpedHexp = 1;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			lumpedHexp = lumpedHexp - allelefreq[iii]*allelefreq[iii];

		double Fst;

		if (lumpedHexp > 0)
			Fst = 1 - (paternalHexp + maternalHexp)/(2*lumpedHexp);
		else
			Fst = -1;

		delete[] allelefreq;
		delete[] paternalallelefreq;
		delete[] maternalallelefreq;

		return Fst;

	}

	void SmoothAllFstValues()
	{
		int iii, jjj, kkk;
		double weightingcoefficient;
		int weightingdistance;
		double weightedFst;
		double numbersummed;

		weightingdistance = static_cast<int>(floor(sqrt(weightingvariance))*3);

		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
			{
				weightedFst = 0;
				numbersummed = 0;
				for (kkk = (jjj-weightingdistance); kkk < jjj+weightingdistance+1; kkk++)
				{
					if (kkk >= 0 && kkk < NumberMarkerLociPerChromosome)
					{
						if (chromosomedata[iii].m_Fst[kkk] >= 0)
						{
							weightingcoefficient = exp(-0.5*(kkk-jjj)*(kkk-jjj)/weightingvariance);
							weightedFst = weightedFst + weightingcoefficient*chromosomedata[iii].m_Fst[kkk];
							numbersummed = numbersummed + weightingcoefficient;
						}
					}
				}
				chromosomedata[iii].m_FstSmooth[jjj] = weightedFst/numbersummed;
			}
		}

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
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			maternalalleleA = Progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[markerA];
			maternalalleleB = Progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[markerB];
			paternalalleleA = Progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[markerA];
			paternalalleleB = Progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[markerB];

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
			allelefreqA[iii] = allelefreqA[iii]/count;
			allelefreqB[iii] = allelefreqB[iii]/count;
			for (jjj = 0; jjj < MaxNumAllelesPerMarker; jjj++)
				jointAB[iii][jjj] = jointAB[iii][jjj]/count;
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			for (jjj = 0; jjj < MaxNumAllelesPerMarker; jjj++)
			{
				if (allelefreqA[iii] > 0 && allelefreqB[jjj] > 0)
				{
					Dij[iii][jjj] = jointAB[iii][jjj] - allelefreqA[iii]*allelefreqB[jjj];
					if (Dij[iii][jjj] < 0)
						Dmax[iii][jjj] = min(allelefreqA[iii]*allelefreqB[jjj],(1-allelefreqA[iii])*(1-allelefreqB[jjj]));
					else
						Dmax[iii][jjj] = min((1-allelefreqA[iii])*allelefreqB[jjj],allelefreqA[iii]*(1-allelefreqB[jjj]));
			
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
						Dprime = Dprime + allelefreqA[iii]*allelefreqB[jjj]*fabs(Dij[iii][jjj])/Dmax[iii][jjj];
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

	void CalculateAdultMarkerAlleleFrequencies(int whichchromosome, int whichmarker, bool includedead)
	{
		int iii;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = 0;

		double allelecounter = 0;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			if (Adult[iii].IsAlive || includedead)
			{
				tempallelefreq[Adult[iii].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
				tempallelefreq[Adult[iii].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
				allelecounter = allelecounter + 2;
			}
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = tempallelefreq[iii]/allelecounter;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii] = tempallelefreq[iii];

		cout << "\n";
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			cout << iii << "\t" << tempallelefreq[iii] << "\n";
	}

	void CalculateProgenyMarkerAlleleFrequencies(int whichchromosome, int whichmarker, bool includedead)
	{
		int iii;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = 0;

		double allelecounter = 0;
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive || includedead)
			{
				tempallelefreq[Progeny[iii].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
				tempallelefreq[Progeny[iii].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
				allelecounter = allelecounter + 2;
			}
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = tempallelefreq[iii]/allelecounter;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii] = tempallelefreq[iii];

		cout << "\n";
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			cout << iii << "\t" << tempallelefreq[iii] << "\n";

	}

	void CalculateAdultProgenyQTLAlleleFrequencies(int whichchromosome, int whichQTL, bool includedead)
	{
		// First make a list of alleles at the locus
		vector<double> AlleleList;
		vector<double> AlleleFreqAverage;
		vector<double> AlleleFreqProgeny;
		vector<double> AlleleFreqAdults;
		int iii, jjj;
		int NumAllelesInList = 0;
		bool Onlist;

		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
		}

		for (iii = 0; iii < PopulationSize; iii++)
		{
			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}

			Onlist = false;
			for (jjj = 0; jjj < NumAllelesInList; jjj++)
			{
				if(Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
					Onlist = true;
			}
			if (!Onlist)
			{
				AlleleList.push_back(Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL]);
				AlleleFreqAverage.push_back(0);
				AlleleFreqProgeny.push_back(0);
				AlleleFreqAdults.push_back(0);
				NumAllelesInList++;
			}
		}

		// Now we should have a list of all the alleles and their number of occurrences

		// Figure out progeny allele frequencies

		double allelecounter;
		allelecounter = 0;
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			if (Progeny[iii].IsAlive || includedead)
			{
				for (jjj = 0; jjj < NumAllelesInList; jjj++)
				{
					if (Progeny[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
						AlleleFreqProgeny[jjj]++;
					if (Progeny[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
						AlleleFreqProgeny[jjj]++;
				}
				allelecounter = allelecounter + 2;
			}
		}

		for (jjj = 0; jjj < NumAllelesInList; jjj++)
			AlleleFreqProgeny[jjj] = AlleleFreqProgeny[jjj]/allelecounter;

		// Figure out adult allele frequencies

		allelecounter = 0;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			if (Adult[iii].IsAlive || includedead)
			{
				for (jjj = 0; jjj < NumAllelesInList; jjj++)
				{
					if (Adult[iii].MaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
						AlleleFreqAdults[jjj]++;
					if (Adult[iii].PaternalChromosome[whichchromosome].QTLeffect[whichQTL] == AlleleList[jjj])
						AlleleFreqAdults[jjj]++;
				}
				allelecounter = allelecounter + 2;
			}
		}

		for (jjj = 0; jjj < NumAllelesInList; jjj++)
			AlleleFreqAdults[jjj] = AlleleFreqAdults[jjj]/allelecounter;

		// Now pass these values into the arrays

		chromosomedata[whichchromosome].m_qtlallelefrequenciesadults[whichQTL].m_numberalleles = NumAllelesInList;
		chromosomedata[whichchromosome].m_qtlallelefrequenciesprogeny[whichQTL].m_numberalleles = NumAllelesInList;
		chromosomedata[whichchromosome].m_qtlallelefrequenciesadults[whichQTL].m_qtlallelefreq.clear();
		chromosomedata[whichchromosome].m_qtlallelefrequenciesprogeny[whichQTL].m_qtlallelefreq.clear();
		chromosomedata[whichchromosome].m_qtlallelefrequenciesadults[whichQTL].m_qtlallelename.clear();
		chromosomedata[whichchromosome].m_qtlallelefrequenciesprogeny[whichQTL].m_qtlallelename.clear();
		
		for (iii = 0; iii < NumAllelesInList; iii++)
		{
			chromosomedata[whichchromosome].m_qtlallelefrequenciesadults[whichQTL].m_qtlallelename.push_back(AlleleList[iii]);
			chromosomedata[whichchromosome].m_qtlallelefrequenciesadults[whichQTL].m_qtlallelefreq.push_back(AlleleFreqAdults[iii]);
			chromosomedata[whichchromosome].m_qtlallelefrequenciesprogeny[whichQTL].m_qtlallelename.push_back(AlleleList[iii]);
			chromosomedata[whichchromosome].m_qtlallelefrequenciesprogeny[whichQTL].m_qtlallelefreq.push_back(AlleleFreqProgeny[iii]);
		}

		cout << "\ncount\tallele\tadult\tprog\n";
		for (iii = 0; iii < NumAllelesInList; iii++)
		{
			cout << iii << "\t" << chromosomedata[whichchromosome].m_qtlallelefrequenciesadults[whichQTL].m_qtlallelename[iii] << "\t" <<
				chromosomedata[whichchromosome].m_qtlallelefrequenciesadults[whichQTL].m_qtlallelefreq[iii] << "\t" <<
				chromosomedata[whichchromosome].m_qtlallelefrequenciesprogeny[whichQTL].m_qtlallelefreq[iii] << "\n";

		}


	} 

	double CalcAdultSexRatio() 
	{
		int iii;
		double numberfemales, numbermales;
		numberfemales = 0;
		numbermales = 0;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			if (Adult[iii].IsFemale)
				numberfemales++;
			else
				numbermales++;
		}
		return numbermales/(numberfemales+numbermales);
	}

	void ChooseProgenySample()
	{
		int iii;
		double dNumberProgenyLeft, dSampleUnfilled;
		double dSampleProb;
		int iNumSampledSoFar;

		dNumberProgenyLeft = NumberOfProgeny;
		dSampleUnfilled = SampleSizeProgeny;
		iNumSampledSoFar = 0;
		for (iii = 0; iii < NumberOfProgeny; iii++)
		{
			dSampleProb = dSampleUnfilled/dNumberProgenyLeft;
			if (genrand() < dSampleProb)
			{
				ProgenySample[iNumSampledSoFar] = iii;
				iNumSampledSoFar++;
				dSampleUnfilled = dSampleUnfilled-1;
			}
			dNumberProgenyLeft = dNumberProgenyLeft-1;
		}
		ActualProgenySampleSize = iNumSampledSoFar;

		//cout << "\nProgeny Sample:";
		//cout << "Number Sampled:\t" << ActualProgenySampleSize << "\n";
		//for (iii = 0; iii < ActualProgenySampleSize; iii++)
		//	cout << ProgenySample[iii] << ", ";
		//cout << "\n";

	}

	void ChooseAdultsSample()
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
			dSampleProb = dSampleUnfilled/dNumberAdultsLeft;
			if (genrand() < dSampleProb)
			{
				AdultSample[iNumSampledSoFar] = iii;
				iNumSampledSoFar++;
				dSampleUnfilled = dSampleUnfilled-1;
			}
			dNumberAdultsLeft = dNumberAdultsLeft-1;
		}
		ActualAdultSampleSize = iNumSampledSoFar;

		//cout << "\nAdult Sample:";
		//cout << "Number Sampled:\t" << ActualAdultSampleSize << "\n";
		//for (iii = 0; iii < ActualAdultSampleSize; iii++)
		//	cout << AdultSample[iii] << ", ";
		//cout << "\n";

	}

	void CalculateAdultMarkerAlleleFrequenciesFromSample(int whichchromosome, int whichmarker)
	{
		int iii, tempindex;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = 0;

		double allelecounter = 0;
		for (iii = 0; iii < ActualAdultSampleSize; iii++)
		{
			tempindex = AdultSample[iii];
			tempallelefreq[Adult[tempindex].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
			tempallelefreq[Adult[tempindex].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
			allelecounter = allelecounter + 2;
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = tempallelefreq[iii]/allelecounter;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii] = tempallelefreq[iii];

		//cout << "\n";
		//for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		//	cout << iii << "\t" << tempallelefreq[iii] << "\n";
	}

	void CalculateProgenyMarkerAlleleFrequenciesFromSample(int whichchromosome, int whichmarker)
	{
		int iii, tempindex;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = 0;

		double allelecounter = 0;
		for (iii = 0; iii < ActualProgenySampleSize; iii++)
		{
			tempindex = ProgenySample[iii];
			tempallelefreq[Progeny[tempindex].MaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
			tempallelefreq[Progeny[tempindex].PaternalChromosome[whichchromosome].MarkerLoci[whichmarker]]++;
			allelecounter = allelecounter + 2;
		}

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			tempallelefreq[iii] = tempallelefreq[iii]/allelecounter;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
			chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii] = tempallelefreq[iii];

		//cout << "\n";
		//for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		//	cout << iii << "\t" << tempallelefreq[iii] << "\n";
	}

	double POfst(int whichchromosome, int whichmarker)
	{
		double* averageallelefreqs;
		averageallelefreqs = new double[MaxNumAllelesPerMarker];
		double HsProgeny, HsAdults, Ht, Fst;
		int iii;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			averageallelefreqs[iii] = (chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii] + 
				chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii])/2;
		}

		HsProgeny = 1;
		HsAdults = 1;
		Ht = 1;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			HsProgeny = HsProgeny - (chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii] *
				chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii]);
			HsAdults = HsAdults - (chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii] *
				chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii]);
			Ht = Ht - (averageallelefreqs[iii]*averageallelefreqs[iii]);
		}

		if (Ht > 0)
			Fst = 1 - (HsProgeny+HsAdults)/(2*Ht);
		else
			Fst = -1.0;

		return Fst;

		delete[] averageallelefreqs;
	}

	double HighestAlleleFreqProgeny(int whichchromosome, int whichmarker)
	{
		double highestfreq = 0;
		int iii;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			if (chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii] > highestfreq)
				highestfreq = chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii];
		}

		return highestfreq;
	}

	double HighestAlleleFreqAdults(int whichchromosome, int whichmarker)
	{
		double highestfreq = 0;
		int iii;

		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			if (chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii] > highestfreq)
				highestfreq = chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii];
		}

		return highestfreq;
	}

	void DetermineSmoothedFstCriticalValue()
	{
		int iii, jjj;
		double meanFst, varFst, stdevFst;
		double count;

		meanFst = 0;
		count = 0;
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
			{
				if (chromosomedata[iii].m_FstSmooth[jjj] >= 0)
				{
					meanFst = meanFst + chromosomedata[iii].m_FstSmooth[jjj];
					count++;
				}
			}
		}
		meanFst = meanFst/count;

		varFst = 0;
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
			{
				if (chromosomedata[iii].m_FstSmooth[jjj] >= 0)
				{
					varFst = varFst + (meanFst - chromosomedata[iii].m_FstSmooth[jjj])*(meanFst - chromosomedata[iii].m_FstSmooth[jjj]);
				}
			}
		}
		varFst = varFst/count;
		stdevFst = sqrt(varFst);

		SmoothedCritValue99 = meanFst + 2.57583 * stdevFst;
		SmoothedCritValue98 = meanFst + 2.32635 * stdevFst;
		SmoothedCritValue95 = meanFst + 1.95996 * stdevFst;
		SmoothedCritValue90 = meanFst + 1.64485 * stdevFst;
		SmoothedCritValue80 = meanFst + 1.28155 * stdevFst;
	}

	double CalculateFstChiSquareValue(int whichchromosome, int whichmarker) 
	{
		int iii;
		// We need the number of alleles per locus

		double actualnumberofalleles;
		double degreesoffreedom;
		double pvalue;

		actualnumberofalleles = 0;
		for (iii = 0; iii < MaxNumAllelesPerMarker; iii++)
		{
			if (chromosomedata[whichchromosome].m_markerallelefrequenciesadults[whichmarker].m_allelefrequency[iii] > 0 ||
				chromosomedata[whichchromosome].m_markerallelefrequenciesprogeny[whichmarker].m_allelefrequency[iii] > 0)
				actualnumberofalleles++;
		}

		double chisqrteststat;
		double Ntotalsamplesize;
		Ntotalsamplesize = SampleSizeAdults + SampleSizeProgeny;

		if (chromosomedata[whichchromosome].m_Fst[whichmarker] >= 0 && chromosomedata[whichchromosome].m_Fst[whichmarker] <= 1)
		{
			chisqrteststat = 2*Ntotalsamplesize*chromosomedata[whichchromosome].m_Fst[whichmarker]*(actualnumberofalleles-1);
			degreesoffreedom = (actualnumberofalleles-1);
			pvalue = 1-chisqr(chisqrteststat, degreesoffreedom);
		}
		else
			pvalue = 1.01;

		return pvalue;
	}

	double FDRchisqrPval(double FDRalpha)
	{
		int iii, jjj;
		int totalmarkers = NumberChromosomes*NumberMarkerLociPerChromosome;
		double totalvariablemarkers;
		vector<double> pvals;
		pvals.reserve(totalmarkers);
		
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for(jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
			{
				pvals.push_back(chromosomedata[iii].m_ChiSqrPval[jjj]);
			}
		}
		
		sort(pvals.begin(),pvals.end());

		totalvariablemarkers = 0;
		for (iii = 0; iii < totalmarkers; iii++)
		{
			if (pvals[iii] != 1.01)
				totalvariablemarkers++;
		}

		double pvalindex, tempFDRp;
		for (iii = static_cast<int>(totalvariablemarkers) - 1; iii >= 0; iii = iii-1)
		{
			pvalindex = iii+1;
			tempFDRp = (pvalindex/totalvariablemarkers)*FDRalpha;
			if (pvals[iii] <= tempFDRp)
				break;
		}

		return tempFDRp;

	}

	void DetermineSignificantLoci()
	{
		int iii, jjj;
		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
			{
				if (chromosomedata[iii].m_ChiSqrPval[jjj] > FDRcriticalvalue)
					chromosomedata[iii].m_significant[jjj] = false;
				else
					chromosomedata[iii].m_significant[jjj] = true;
			}
		}

	}

	void FindAllLocalSmoothedFstMaxima(int halfregressioninterval)
	{
		int iii, jjj, kkk;

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

		for (iii = 0; iii < NumberChromosomes; iii++)
		{
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
					meanY = meanY + chromosomedata[iii].m_FstSmooth[kkk];
					numbersummed++;
				}
				meanX = meanX/numbersummed;
				meanY = meanY/numbersummed;

				covarXY = 0;
				varX = 0;
				for (kkk = jjj-halfregressioninterval; kkk < jjj + halfregressioninterval; kkk++)
				{
					tempdub = kkk;
					covarXY = covarXY + (tempdub - meanX)*(chromosomedata[iii].m_FstSmooth[kkk] - meanY);
					varX = varX + (chromosomedata[iii].m_FstSmooth[kkk] - meanY)*(chromosomedata[iii].m_FstSmooth[kkk] - meanY);
				}
				covarXY = covarXY/(numbersummed-1);
				varX = varX/(numbersummed-1);

				currentregression = covarXY/varX;
				isamax = false;
				if (previousregression > -1000000)
				{
					regressionproduct = currentregression*previousregression;
					if (regressionproduct <= 0 && previousregression > 0)
						isamax = true;
				}

				if (isamax)
				{
					chromosomedata[iii].m_smoothedFstpeak[jjj-1] = true;
					if (chromosomedata[iii].m_FstSmooth[jjj-1] < lowestpeak)
						lowestpeak = chromosomedata[iii].m_FstSmooth[jjj-1];
					if (chromosomedata[iii].m_FstSmooth[jjj-1] > highestpeak)
						highestpeak = chromosomedata[iii].m_FstSmooth[jjj-1];
				}
				else
				{
					chromosomedata[iii].m_smoothedFstpeak[jjj] = false;
				}

				previousregression = currentregression;
			} // jjj

			// Now we have to deal with the ends
			// See if there's a local maximum at the end, and if it's in the 
			// range of the other peaks, count it.  Otherwise don't.

			highestpointinendsegment = 0;
			for (jjj = 0; jjj < halfregressioninterval*2; jjj++)
			{
				chromosomedata[iii].m_smoothedFstpeak[jjj] = false;
				if (chromosomedata[iii].m_FstSmooth[jjj] > highestpointinendsegment)
				{
					highestpointinendsegment = chromosomedata[iii].m_FstSmooth[jjj];
					indexofhighestpoint = jjj;
				}
			}

			if (indexofhighestpoint != halfregressioninterval*2-1)
			{
				if (highestpointinendsegment > lowestpeak)
					chromosomedata[iii].m_smoothedFstpeak[indexofhighestpoint] = true;
			}

			highestpointinendsegment = 0;
			for (jjj = NumberMarkerLociPerChromosome - halfregressioninterval*2; jjj < NumberMarkerLociPerChromosome; jjj++)
			{
				chromosomedata[iii].m_smoothedFstpeak[jjj] = false;
				if (chromosomedata[iii].m_FstSmooth[jjj] > highestpointinendsegment)
				{
					highestpointinendsegment = chromosomedata[iii].m_FstSmooth[jjj];
					indexofhighestpoint = jjj;
				}
			}

			if (indexofhighestpoint != NumberMarkerLociPerChromosome - halfregressioninterval*2)
			{
				if (highestpointinendsegment > lowestpeak)
				{
					chromosomedata[iii].m_smoothedFstpeak[indexofhighestpoint] = true;
				}
			}

			// Remove all peaks that are adjacent to one another (because in this case
			// it's really just one peak).  In each case, keep the one with the higher Fst
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome - 1; jjj++)
			{
				if (chromosomedata[iii].m_smoothedFstpeak[jjj] && chromosomedata[iii].m_smoothedFstpeak[jjj+1])
				{
					if (chromosomedata[iii].m_Fst[jjj] > chromosomedata[iii].m_Fst[jjj+1])
						chromosomedata[iii].m_smoothedFstpeak[jjj+1] = false;
					else
						chromosomedata[iii].m_smoothedFstpeak[jjj] = false;
				}
			}

		} // iii



		


	}

	void CompilePeakData()
	{
		int iii, jjj, kkk;
		double highestFstnearPeak;
		int highestFstnearPeakindex;
		bool nearbysignificance;
		int platstart, platstop;
		int smallestQTLdist;
		int closestQTLloc;

		for (iii = 0; iii < NumberChromosomes; iii++)
		{
			// Clear all previous data
			peakdata[iii].highestFstonplateau.clear();
			peakdata[iii].highestFstIndex.clear();
			peakdata[iii].peaklocation.clear();
			peakdata[iii].peaksmoothedFst.clear();
			peakdata[iii].sig80CI.clear();
			peakdata[iii].sig90CI.clear();
			peakdata[iii].sig95CI.clear();
			peakdata[iii].sig98CI.clear();
			peakdata[iii].sig99CI.clear();
			peakdata[iii].sig95BS.clear();
			peakdata[iii].sig99BS.clear();
			peakdata[iii].sigFDR.clear();
			peakdata[iii].nearestQTL.clear();
			peakdata[iii].distancetoQTL.clear();

			peakdata[iii].numberofpeaks = 0;
			for (jjj = 0; jjj < NumberMarkerLociPerChromosome; jjj++)
			{
				if (chromosomedata[iii].m_smoothedFstpeak[jjj])  // Is this marker a peak?
				{
					peakdata[iii].peaksmoothedFst.push_back(chromosomedata[iii].m_FstSmooth[jjj]);
					peakdata[iii].peaklocation.push_back(jjj);

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
						if (chromosomedata[iii].m_Fst[kkk] > highestFstnearPeak)
						{
							highestFstnearPeak = chromosomedata[iii].m_Fst[kkk];
							highestFstnearPeakindex = kkk;
						}
					}
					peakdata[iii].highestFstonplateau.push_back(highestFstnearPeak);
					peakdata[iii].highestFstIndex.push_back(highestFstnearPeakindex);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_significant[kkk])
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sigFDR.push_back(true);
					else
						peakdata[iii].sigFDR.push_back(false);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_FstSmooth[kkk] >= SmoothedCritValue99)
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sig99CI.push_back(true);
					else
						peakdata[iii].sig99CI.push_back(false);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_FstSmooth[kkk] >= SmoothedCritValue98)
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sig98CI.push_back(true);
					else
						peakdata[iii].sig98CI.push_back(false);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_FstSmooth[kkk] >= SmoothedCritValue95)
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sig95CI.push_back(true);
					else
						peakdata[iii].sig95CI.push_back(false);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_FstSmooth[kkk] >= SmoothedCritValue90)
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sig90CI.push_back(true);
					else
						peakdata[iii].sig90CI.push_back(false);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_FstSmooth[kkk] >= SmoothedCritValue80)
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sig80CI.push_back(true);
					else
						peakdata[iii].sig80CI.push_back(false);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_FstSmooth[kkk] >= bs99cutoff)
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sig99BS.push_back(true);
					else
						peakdata[iii].sig99BS.push_back(false);

					nearbysignificance = false;
					for (kkk = platstart; kkk < platstop; kkk++)
					{
						if (chromosomedata[iii].m_FstSmooth[kkk] >= bs95cutoff)
							nearbysignificance = true;
					}
					if (nearbysignificance)
						peakdata[iii].sig95BS.push_back(true);
					else
						peakdata[iii].sig95BS.push_back(false);



					closestQTLloc = -1;
					smallestQTLdist = NumberMarkerLociPerChromosome;
					for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
					{
						if (abs(jjj - QTLlocus[iii].QTLlocation[kkk]) < smallestQTLdist)
						{
							smallestQTLdist = abs(jjj - QTLlocus[iii].QTLlocation[kkk]);
							closestQTLloc = QTLlocus[iii].QTLlocation[kkk];
						}
					}
					peakdata[iii].distancetoQTL.push_back(smallestQTLdist);
					peakdata[iii].nearestQTL.push_back(closestQTLloc);

					peakdata[iii].numberofpeaks++;

				}

			}
		}

	}

	void BootstrapSmoothedFstValues()
	{
		// This bootstrapping procedure is similar to that of Catchen.
		// However, ours is simplified by the fact that markers are evenly spaced.
		// For each bootstrap rep, pick loci at random for each raw Fst value.
		// Then average those to get the "smoothed" value.

		int iii, reps;
		int chr_rand, mark_rand;
		double weightingcoefficient;
		int weightingdistance;
		double weightedFst;
		double numbersummed;
		vector<double> bsSmoothedFstList;
		int percent95, percent99;

		weightingdistance = static_cast<int>(floor(sqrt(weightingvariance))*3);

		for (reps = 0; reps < bootstrapreps; reps++)
		{
			weightedFst = 0;
			numbersummed = 0;
			for (iii = 0; iii < (2*weightingdistance)+1; iii++)
			{
				chr_rand = randnum(NumberChromosomes);
				mark_rand = randnum(NumberMarkerLociPerChromosome);
				if (chromosomedata[chr_rand].m_Fst[mark_rand] >= 0)
				{
					weightingcoefficient = exp(-0.5*(iii-weightingdistance)*(iii-weightingdistance)/weightingvariance);
					weightedFst = weightedFst + weightingcoefficient*chromosomedata[chr_rand].m_Fst[mark_rand];
					numbersummed = numbersummed + weightingcoefficient;
				}
			}
			bsSmoothedFstList.push_back(weightedFst/numbersummed);
		}

		sort (bsSmoothedFstList.begin(),bsSmoothedFstList.end());
		percent95 = 95*bootstrapreps/100;
		percent99 = 99*bootstrapreps/100;

		bs95cutoff = bsSmoothedFstList[percent95];
		bs99cutoff = bsSmoothedFstList[percent99];

	}

	void StandardizeQTLeffects() 
	{
		int iii, jjj, kkk;
		double BVmean, BVvar, BVstdev, dpopsize, BVmeanAdjustment;
		dpopsize = PopulationSize;

		BVmean = 0;
		BVvar = 0;
		for (iii = 0; iii < PopulationSize; iii++)
			BVmean = BVmean + Adult[iii].BreedingValue;
		BVmean = BVmean/dpopsize;

		for (iii = 0; iii < PopulationSize; iii++)
			BVvar = BVvar + (BVmean - Adult[iii].BreedingValue)*(BVmean - Adult[iii].BreedingValue);
		BVvar = BVvar/dpopsize;

		BVstdev = sqrt(BVvar);

		double dtotQTLs;
		dtotQTLs = NumberChromosomes*NumberQTLsPerChromosome;

		BVmeanAdjustment = BVmean/(dtotQTLs*2);

		for (iii = 0; iii < PopulationSize; iii++)
		{
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
				{
					Adult[iii].MaternalChromosome[jjj].QTLeffect[kkk] = (Adult[iii].MaternalChromosome[jjj].QTLeffect[kkk]-BVmeanAdjustment)/BVstdev;
					Adult[iii].PaternalChromosome[jjj].QTLeffect[kkk] = (Adult[iii].PaternalChromosome[jjj].QTLeffect[kkk]-BVmeanAdjustment)/BVstdev;
				}
			}
		}

		double tempenveffect;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			tempenveffect = Adult[iii].Phenotype - Adult[iii].BreedingValue;
			Adult[iii].BreedingValue = 0;
			for (jjj = 0; jjj < NumberChromosomes; jjj++)
			{
				for (kkk = 0; kkk < NumberQTLsPerChromosome; kkk++)
				{
					Adult[iii].BreedingValue = Adult[iii].BreedingValue + Adult[iii].MaternalChromosome[jjj].QTLeffect[kkk];
					Adult[iii].BreedingValue = Adult[iii].BreedingValue + Adult[iii].PaternalChromosome[jjj].QTLeffect[kkk];
				}

			}
			Adult[iii].Phenotype = Adult[iii].BreedingValue + tempenveffect;
		}

	}



};