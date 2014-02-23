// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:

#include <OpenMS/CONCEPT/ProgressLogger.h>

#include<OpenMS/DATASTRUCTURES/String.h> 

#include<OpenMS/KERNEL/DPeak.h>
#include<OpenMS/KERNEL/MSExperiment.h>

#include<OpenMS/FILTERING/BASELINE/TopHatFilter.h>
#include<OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include<OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include<OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include<OpenMS/FORMAT/FileHandler.h>
#include<OpenMS/FORMAT/MzDataFile.h>
#include<OpenMS/FORMAT/MzXMLFile.h>

#include<OpenMS/SYSTEM/StopWatch.h>

#include<iostream>  
#include<fstream>
#include<algorithm>
#include<cstdlib>
#include<stdio.h>
#include<math.h>
#include<limits.h>

#include<gsl/gsl_statistics.h>

using namespace OpenMS;
using namespace std;

typedef MSExperiment< >::SpectrumType SpectrumType;

// number of quality descriptors
unsigned int NUM_QUAL_DESC = 27; //35;

//#define DB_COMPSTATS

// skip weird scans, but print error message
void writeEmptyScan(unsigned int scan_nr, std::ofstream& out)
{
	cout << "Error in scan " << scan_nr << endl;
	out << scan_nr << " ";
	for (unsigned int i=0; i<NUM_QUAL_DESC;++i)	out << "NA ";

	out << endl; 
}

double calculateXrea(vector<double>& intensities)
{
  vector<double> cumulatedNormalizedIntensities;
  double integral = 0; 
  double tic = 0;
  double cumulatedIntensity = 0;
     
  for (vector<double>::const_iterator citer = intensities.begin(); citer != intensities.end(); ++citer)
  {
    tic += *citer;
  }
        
  // collect cumulated intensities
  for (size_t i = 0; i < intensities.size(); ++i)
  {
    cumulatedIntensity += intensities[i];
    cumulatedNormalizedIntensities.push_back(cumulatedIntensity/tic);
  }
    
    
  unsigned int bins = cumulatedNormalizedIntensities.size();
  for (unsigned int i = 0; i < bins; ++i)
  {
    integral += cumulatedNormalizedIntensities[i] / (bins+1);
  }
    
  double alpha = intensities.at(intensities.size() - 1) / tic;
  double xrea = (0.5 - integral) / ( 0.5 + alpha);
    
  return xrea;
}

double computeDiff2SmoothedScan(SpectrumType& s)
{
	GaussFilter gf;
	Param p;
	p.setValue("gaussian_width",0.2);
	gf.setParameters(p);
	
	SpectrumType new_scan;
	gf.filter(s,new_scan);

	double eucl_diff = 0.0;

	#ifdef DB_COMPSTATS
	String of1 = "smoothed_scan_" + String(s.getRT());
	String of2 = "real_scan_" + String(s.getRT());
	ofstream os(of1.c_str());
	ofstream oo(of2.c_str());
	#endif
	for (unsigned int i=0;i<s.size();++i)
	{
		eucl_diff += pow(s[i].getIntensity() - new_scan[i].getIntensity(),2);
	
		#ifdef DB_COMPSTATS
		os << new_scan[i].getMZ() << " " << new_scan[i].getIntensity() << endl;
		oo << s[i].getMZ() << " " << s[i].getIntensity() << endl;
		#endif	
	}
	#ifdef DB_COMPSTATS
	os.close();
	oo.close();
	#endif

	if (eucl_diff == 0.0) eucl_diff = 0.0001;

	return sqrt(eucl_diff);
}

// compute euclidean distance between original and baseline-removed scan
double computeDiff2Baseline(SpectrumType& s)
{
	// sampling depth
	double spacing = 0.1;
	// size of structuring element 
	double struc_elem_length = 2.5;

	// resampling
	LinearResampler lnr;
	Param p1;
	p1.setValue("spacing",spacing);
	lnr.setParameters(p1);

	SpectrumType res_scan;
	lnr.raster(s,res_scan);

	// tophat
	SpectrumType new_scan;
	TopHatFilter thf;
	Param p2;
	p2.setValue("struc_elem_length",struc_elem_length);
	thf.setParameters(p2);
	thf.filter(res_scan,new_scan);
	
	double eucl_diff = 0.0;
	#ifdef DB_COMPSTATS
	String of1 = "blremoved_scan_" + String(s.getRT());
	String of2 = "resampled_scan_" + String(s.getRT());
	ofstream os(of1.c_str());
	ofstream oo(of2.c_str());
	#endif
	for (unsigned int i=0;i<res_scan.size();++i)
	{
		eucl_diff += pow(res_scan[i].getIntensity() - new_scan[i].getIntensity(),2);
			
		#ifdef DB_COMPSTATS
		os << new_scan[i].getMZ() << " " <<  new_scan[i].getIntensity() << endl;
		oo << res_scan[i].getMZ()  <<  " " << res_scan[i].getIntensity() << endl;
		#endif
	}
	#ifdef DB_COMPSTATS
	os.close();
	oo.close();
	#endif
	if (eucl_diff == 0.0) eucl_diff = 0.0001;
	return sqrt(eucl_diff);
}

int main(int argc, char *argv[])
{        
  if (argc != 2)
  {
    cout << "ComputeScanWiseStats <inputfile>" << endl;
    return (-1);
  }

  String infile = String(argv[1]);
 
  MSExperiment< > exp;
 	
	FileHandler fh;
	FileHandler::Type in_type = fh.getTypeByFileName(infile);

	if (in_type==FileHandler::UNKNOWN)
	{
		in_type = fh.getTypeByContent(infile);
	}
	if (in_type==FileHandler::UNKNOWN)
	{
		cout << "Error: Could not determine input file type!" << endl;
		return (-1);
	}

	if (! fh.loadExperiment(infile,exp,in_type,ProgressLogger::CMD) )
	{
		cout << "Unsupported or corrupt input file. Aborting!" << endl;
		return (-1);			
	}
        
  vector<double> intensities;
  vector<double> mzvals;
	//vector<double> snvals;

  double it_sum = 0; // sum of intensities
	UInt scan_num = 0; // scan number
  
	double min_it    = numeric_limits<double>::max();
	double sn_min_it = 0.0;
	
	double max_it    = 0.0;
	double sn_max_it = 0.0;

  int n = infile.find('.',0);
  String ending = "_scanwise_stats.out";
  String outfile = infile.replace(n,infile.size(),ending);
	
	// write R header
  ofstream out( outfile.c_str() );
  out << "scan rt quantile1_it max_it mean_it median_it min_it quantile3_it variance_it ";
	out << "skewness_it kurtosis_it itsum quantile1_mz max_mz mean_mz median_mz min_mz ";
	out << "quantile3_mz variance_mz skewness_mz kurtosis_mz ";
	out << "sn_min_it sn_max_it num_points xrea bl_diff sm_diff" << endl;
	//out << "quantile1_sn max_sn mean_sn median_sn min_sn quantile3_sn variance_sn skewness_sn ";
	//out << "kurtosis_sn num_points xrea bl_diff sm_diff" << endl;

	ProgressLogger pl;
	pl.setLogType(ProgressLogger::CMD);
	pl.startProgress(0,exp.size() ,"stats");

  for (MSExperiment< >::iterator spec = exp.begin(); spec != exp.end(); ++spec)
  {
		pl.setProgress(scan_num++);
    if (spec->getMSLevel() != 1 || spec->size() < 20) continue;
		
    it_sum = 0;
		intensities.clear();
		mzvals.clear();
		//snvals.clear();

		SignalToNoiseEstimatorMeanIterative< > sn_estim;
		sn_estim.init(spec->begin(),spec->end());

		for (SpectrumType::const_iterator peak = spec->begin(); peak != spec->end(); ++peak)
    {
      intensities.push_back(peak->getIntensity());
      mzvals.push_back( peak->getMZ() );
      it_sum += peak->getIntensity();
			//snvals.push_back(sn_estim.getSignalToNoise(peak));

			if (peak->getIntensity() > max_it)
			{
				max_it    = peak->getIntensity();
				sn_max_it = sn_estim.getSignalToNoise(peak);
			}
			if (peak->getIntensity() < min_it)
			{
				min_it    = peak->getIntensity();
				sn_min_it = sn_estim.getSignalToNoise(peak);
			}
    }
                
    if (it_sum == 0) 
		{
			cout << "Zero intensity sum !! " << endl;
			writeEmptyScan(scan_num,out);
     	continue; 
		} 

		double bl_diff = computeDiff2Baseline(*spec);
		double sm_diff = computeDiff2SmoothedScan(*spec);

    sort(intensities.begin(), intensities.end());
    //sort(snvals.begin(), snvals.end());
    sort(mzvals.begin(), mzvals.end());
    
		double xrea = calculateXrea(intensities);
       
    unsigned int sz_int    = intensities.size();
    unsigned int sz_mzvals = mzvals.size();
    //unsigned int sz_snvals = snvals.size();
        
    double mean_int    = gsl_stats_mean(&intensities[0], 1, sz_int);
    double mean_mzvals = gsl_stats_mean(&mzvals[0], 1, sz_mzvals);
    //double mean_snvals = gsl_stats_mean(&snvals[0], 1, sz_snvals);

		// write summary stats
		// intensity
    out << scan_num << " ";
		out << spec->getRT() << " ";
		out << gsl_stats_quantile_from_sorted_data(&intensities[0], 1, sz_int, 0.25) << " ";
    out << gsl_stats_max(&intensities[0], 1, sz_int) << " ";
    out << gsl_stats_mean(&intensities[0], 1, sz_int) << " ";
    out << gsl_stats_median_from_sorted_data(&intensities[0], 1, sz_int)  << " ";
    out << gsl_stats_min(&intensities[0], 1, sz_int)  << " ";
    out << gsl_stats_quantile_from_sorted_data(&intensities[0], 1, sz_int, 0.75) << " ";
    out << gsl_stats_variance_m(&intensities[0], 1, sz_int, mean_int) << " ";
		out << gsl_stats_skew(&intensities[0], 1, sz_int) << " ";
		out << gsl_stats_kurtosis(&intensities[0], 1, sz_int) << " ";
		out << it_sum << " ";
    
		// m/z
    out << gsl_stats_quantile_from_sorted_data(&mzvals[0], 1, sz_mzvals, 0.25) << " ";
    out << gsl_stats_max(&mzvals[0], 1, sz_mzvals) << " ";
    out << gsl_stats_mean(&mzvals[0], 1, sz_mzvals) << " ";
    out << gsl_stats_median_from_sorted_data(&mzvals[0], 1, sz_mzvals)  << " ";
    out << gsl_stats_min(&mzvals[0], 1, sz_mzvals)  << " ";
    out << gsl_stats_quantile_from_sorted_data(&mzvals[0], 1, sz_mzvals, 0.75) << " ";
    out << gsl_stats_variance_m(&mzvals[0], 1, sz_mzvals, mean_mzvals) << " ";
    out << gsl_stats_skew(&mzvals[0], 1, sz_mzvals) << " ";
    out << gsl_stats_kurtosis(&mzvals[0], 1, sz_mzvals) << " ";

		// s/n
		out << sn_min_it << " ";
		out << sn_max_it << " ";
    
		/*out << gsl_stats_quantile_from_sorted_data(&snvals[0], 1, sz_snvals, 0.25) << " ";
    out << gsl_stats_max(&snvals[0], 1, sz_snvals) << " ";
    out << gsl_stats_mean(&snvals[0], 1, sz_snvals) << " ";
    out << gsl_stats_median_from_sorted_data(&snvals[0], 1, sz_snvals)  << " ";
    out << gsl_stats_min(&snvals[0], 1, sz_snvals)  << " ";
    out << gsl_stats_quantile_from_sorted_data(&snvals[0], 1, sz_snvals, 0.75) << " ";
    out << gsl_stats_variance_m(&snvals[0], 1, sz_snvals, mean_snvals) << " ";
    out << gsl_stats_skew(&snvals[0], 1, sz_snvals) << " ";
    out << gsl_stats_kurtosis(&snvals[0], 1, sz_snvals) << " ";*/

		out << sz_int << " ";
    out << xrea << " ";
		out << bl_diff << " ";
		out << sm_diff << endl;

  } // end of (for all scans)
  
  out.close();  
}
