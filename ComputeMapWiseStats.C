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

//#define DB_COMPSTATS

double calculateXrea(vector<double>& intensities,double rt)
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
	//String fn = "cum_intensities_" + String(rt);
	//ofstream of(fn.c_str());
  for (unsigned int i = 0; i < bins; ++i)
  {
    integral += cumulatedNormalizedIntensities[i] / (bins+1);
  	//of << cumulatedNormalizedIntensities[i] << endl;
	}
	//of.close();
    
  double alpha = intensities.at(intensities.size() - 1) / tic;
	// alpha is a correction term, see Na and Paek (2006)
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
    cout << "ComputeMapWiseStats <inputfile>" << endl;
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

	vector<double> it_means;
	vector<double> it_medians;
	vector<double> it_mins;
	vector<double> it_maxs;
        
	vector<double> mz_means;
	vector<double> mz_medians;
	vector<double> mz_mins;
	vector<double> mz_maxs;

	vector<double> sn_means;
	vector<double> sn_medians;
	vector<double> sn_mins;
	vector<double> sn_maxs;

	vector<double> szs;
	vector<double> xreas;
	vector<double> bldiffs;
	vector<double> smdiffs;

  double it_sum = 0; // sum of intensities
	UInt scan_num = 0; // scan number
	vector<double> tic;

  int n = infile.find('.',0);
  String ending_mw = "_mapwise_stats.out";
  String outfile_mw = infile.replace(n,infile.size(),ending_mw);

	String ending_sw = "_scanwise_stats.out";
  String outfile_sw = infile.replace(n,infile.size(),ending_sw);
	
	// write R header for scanwise stats
  ofstream out( outfile_sw.c_str() );
	out << "it_max it_min it_mean it_median ";
	out << "mz_max mz_min mz_mean mz_median ";
	out << "sn_max sn_min sn_mean sn_median ";
	out << "szs xreas bldiffs smdiffs ";
	out << "tic_min tic_max tic_skew tic_kurtosis" << endl;
	
	ProgressLogger pl;
	pl.setLogType(ProgressLogger::CMD);
	pl.startProgress(0,exp.size() ,"stats");

  for (MSExperiment< >::iterator spec = exp.begin(); spec != exp.end(); ++spec)
  {
		pl.setProgress(scan_num++);
    if (spec->getMSLevel() != 1 || spec->size() < 20) continue;
		
    it_sum = 0;
  	
		vector<double> intensities;
  	vector<double> mzvals;
		vector<double> snvals;

		SignalToNoiseEstimatorMeanIterative< > sn_estim;
		sn_estim.init(spec->begin(),spec->end());

		for (SpectrumType::const_iterator peak = spec->begin(); peak != spec->end(); ++peak)
    {
      intensities.push_back(peak->getIntensity());
      mzvals.push_back( peak->getMZ() );
      it_sum += peak->getIntensity();
			snvals.push_back(sn_estim.getSignalToNoise(peak));
    }
  	tic.push_back(it_sum);
	
    if (it_sum == 0) 
		{
			cout << "Zero intensity sum !! " << endl;
     	continue; 
		} 

		double bl_diff = computeDiff2Baseline(*spec);
		double sm_diff = computeDiff2SmoothedScan(*spec);

    sort(intensities.begin(), intensities.end());
   	sort(snvals.begin(), snvals.end());
    sort(mzvals.begin(), mzvals.end());
    
		double xrea = calculateXrea(intensities,spec->getRT());
       
    unsigned int sz_int    = intensities.size();
    unsigned int sz_mzvals = mzvals.size();
    unsigned int sz_snvals = snvals.size();

		double i = 0.0;    
	
		// write and store scanwise stats
		// intensity
		i = gsl_stats_max(&intensities[0], 1, sz_int);
		out << i << " ";
		it_maxs.push_back(i);

		i = gsl_stats_mean(&intensities[0], 1, sz_int);
		out << i << " ";
		it_means.push_back(i);
		
		i = gsl_stats_median_from_sorted_data(&intensities[0], 1, sz_int);
		out << i << " ";
		it_medians.push_back(i);
		
		i = gsl_stats_min(&intensities[0], 1, sz_int);
		out << i << " ";
		it_mins.push_back(i);

		
		// m/z
		i = gsl_stats_max(&mzvals[0], 1, sz_mzvals);
		out << i << " ";
		mz_maxs.push_back(i);
		
		i = gsl_stats_mean(&mzvals[0], 1, sz_mzvals);
		out << i << " ";
		mz_means.push_back(i);
		
		i = gsl_stats_median_from_sorted_data(&mzvals[0], 1, sz_mzvals);
		out << i << " ";
		mz_medians.push_back(i);
		
		i = gsl_stats_min(&mzvals[0], 1, sz_mzvals);
		out << i << " ";
		mz_mins.push_back(i);

		// s/n
		i = gsl_stats_max(&snvals[0], 1, sz_snvals);
		out << i << " ";
		sn_maxs.push_back(i);

		i = gsl_stats_mean(&snvals[0], 1, sz_snvals);
		out << i << " ";
		sn_means.push_back(i);

		i = gsl_stats_median_from_sorted_data(&snvals[0], 1, sz_snvals);
		out << i << " ";
		sn_medians.push_back(i);
		
		i = gsl_stats_min(&snvals[0], 1, sz_snvals);
		out << i << " ";
		sn_mins.push_back(i);

		// others
		out << sz_int << " ";
		szs.push_back(sz_int);
		out << xrea << " ";
		xreas.push_back(xrea);
		out << bl_diff << " ";
		bldiffs.push_back(bl_diff);
		out << sm_diff << endl;
		smdiffs.push_back(sm_diff);

  } // end of (for all scans)
	
	out.close();	

	sort(it_maxs.begin(),it_maxs.end());
	sort(it_mins.begin(),it_mins.end());
	sort(it_means.begin(),it_means.end());
	sort(it_medians.begin(),it_medians.end());

	sort(mz_maxs.begin(),mz_maxs.end());
	sort(mz_mins.begin(),mz_mins.end());
	sort(mz_means.begin(),mz_means.end());
	sort(mz_medians.begin(),mz_medians.end());
	
	sort(sn_maxs.begin(),sn_maxs.end());
	sort(sn_mins.begin(),sn_mins.end());
	sort(sn_means.begin(),sn_means.end());
	sort(sn_medians.begin(),sn_medians.end());

	sort(szs.begin(),szs.end());
	sort(xreas.begin(),xreas.end());
	sort(bldiffs.begin(),bldiffs.end());
	sort(smdiffs.begin(),smdiffs.end());
	
	// write R header for mapwise statistics
  out.open( outfile_mw.c_str() );
	out << "it_max it_min it_mean it_median ";
	out << "mz_max mz_min mz_mean mz_median ";
	out << "sn_max sn_min sn_mean sn_median ";
	out << "szs xreas bldiffs smdiffs ";
	out << "tic_min tic_max tic_skew tic_kurtosis" << endl;

	// write stats
  out << gsl_stats_median_from_sorted_data(&it_maxs[0], 1, it_maxs.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&it_mins[0], 1, it_mins.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&it_means[0], 1, it_means.size())  << " ";
  out << gsl_stats_median_from_sorted_data(&it_medians[0], 1, it_medians.size())  << " ";
    
  out << gsl_stats_median_from_sorted_data(&mz_maxs[0], 1, mz_maxs.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&mz_mins[0], 1, mz_mins.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&mz_means[0], 1, mz_means.size())  << " ";
  out << gsl_stats_median_from_sorted_data(&mz_medians[0], 1, mz_medians.size())  << " ";

  out << gsl_stats_median_from_sorted_data(&sn_maxs[0], 1, sn_maxs.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&sn_mins[0], 1, sn_mins.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&sn_means[0], 1, sn_means.size())  << " ";
  out << gsl_stats_median_from_sorted_data(&sn_medians[0], 1, sn_medians.size())  << " ";
  
  out << gsl_stats_median_from_sorted_data(&szs[0], 1, szs.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&xreas[0], 1, xreas.size()) << " ";
  out << gsl_stats_median_from_sorted_data(&bldiffs[0], 1, bldiffs.size())  << " ";
  out << gsl_stats_median_from_sorted_data(&smdiffs[0], 1, smdiffs.size())  << " ";

	out << gsl_stats_min(&tic[0],1,tic.size()) << " ";
	out << gsl_stats_max(&tic[0],1,tic.size()) << " ";
	out << gsl_stats_skew(&tic[0],1,tic.size()) << " ";
	out << gsl_stats_kurtosis(&tic[0],1,tic.size()) << endl;

  out.close();  
}
