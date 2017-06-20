////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                              ///
/// Evolutionary model ///
/// ------------------------------------------------------------------------------------------------------------ ///
///                                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <random>

#include <gsl/gsl_rng.h>

using namespace std;

// to build the program:
//   g++ SVModel.cpp -o SVModel -L/usr/local/lib/ -lgsl -I/usr/local/include
//
// to run the model:
//   ./SVModel Seed=1 Model=2

//////////////////////////////////////////////////////////
///                                                    ///
///   MAIN                                             ///
///                                                    ///
//////////////////////////////////////////////////////////
int main (int argc, char * const argv[]) {


  int InternalN;
  int RpNumber;
	
  ////////////////////////////////////////////////////////////
  ///Generating parameters and random seed from bash script///
  ////////////////////////////////////////////////////////////
	
  for (int CommandLineCounter=1; CommandLineCounter < argc; CommandLineCounter++) {
				
    for (int StringPosition=0; StringPosition < strlen(argv[CommandLineCounter]); StringPosition++) {
						
      int EqualPosition = -1;
			
      if (argv[CommandLineCounter][StringPosition]=='=') {

	EqualPosition = StringPosition;
				
	//std::cout << "= is at: " << StringPosition << std::endl;  				
				
      }
			
      
      if (EqualPosition > 1) {
	
	int VariableNameLength = EqualPosition;

	
	char VariableNameChar[VariableNameLength];
				
	for (int i=0; i<EqualPosition; i++) {
					
	  VariableNameChar[i] = argv[CommandLineCounter][i];

					
	}
				
	int VariableStringLength = strlen(argv[CommandLineCounter])-EqualPosition;

	char VariableChar[VariableStringLength];
				
	for (int i=0; i<VariableStringLength; i++) {
					
	  VariableChar[i] = argv[CommandLineCounter][i+1+EqualPosition];
					
	}
								
	double t1,t2; 
	char Bunk[2];
	Bunk[0] = VariableNameChar[0];
	Bunk[1] = VariableNameChar[1];
	
	if(Bunk[0]=='S'){
	  
	  t1 = strtod(VariableChar,NULL);
	  InternalN = int(t1);
	  
	  //std::cout << "Seed " << InternalN << std::endl;					
					
	}
	else if(Bunk[0]=='M'){
	  
	  RpNumber = strtod(VariableChar,NULL);
	  RpNumber = atoi(VariableChar);
	  //std::cout << "Model " << RpNumber << std::endl;
					
	}
									
				
      }
    }
  }

  // initialise rng
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  

  //fixed configuration parameters

  int NChr = 3;
  double lowB = 0.05;
  double chrlength = 1.0e+8;
  double Curvemax = 1.0;
  //int Nparam = 14;
  int Nparticles = 1;
  //int len_stat = 1000;

  
  int seed = 253564389 - 2000 - 6*InternalN - 100 - 400;	     
  srand(seed);	

  ////////////////////////////////////////////////////////////
  ///   Sample from the prior                              ///
  ////////////////////////////////////////////////////////////

  int Npop = (int) 1000; //2000 + (8000-2000)*gsl_rng_uniform (r);
  int ngen = (int) 1000; //1000 + (20000-1000)*gsl_rng_uniform (r);
  
  // pksv?
  double g_d = gsl_rng_uniform (r);
  // cell-wide mutation probability
  double mu_i = 1.0;
  // probability of translocation per generation
  double p_tran = gsl_rng_uniform (r);
  // mu_k?
  double mu_i_range = gsl_rng_uniform (r);
  
  // alpha_sv: gradient of the sv probability function with length
  double svgradient = 0.01 + (0.9-0.01)*gsl_rng_uniform (r);
  double svp2 = 5 + (7-5)*gsl_rng_uniform (r);
  
  // alpha_trans: gradient of the tranlocation probability function with length
  double svtransgrad = 0.01 + (0.9-0.01)*gsl_rng_uniform (r);
  double trp2 = 5 + (7-5)*gsl_rng_uniform (r);
  
  // clim: the gradient of the fitness functions
  double threshold_g = 200 + (300-200)*gsl_rng_uniform (r);
  double maxchr =  250e6;
  double minchr = 1e6;
  
  // max and minimum size of SVs
  double SV_min = 1000;  
  double SV_max = 0.5*250e6;
  

  // probability of genome doubling
  double gnmdou = gsl_rng_uniform (r);
    
  std::cout << "Sampled parameters:\n" << std::endl;
  std::cout << "\tNpop:" << "\t" << Npop << std::endl;
  std::cout << "\tngen:" << "\t" << ngen << std::endl;

  
 

  ////////////////////////////////////////////////////////////
  ///   Loop over the particles                            ///
  ////////////////////////////////////////////////////////////
    
  for (int bd = 0; bd < Nparticles; bd++) {

    vector<vector<vector<double> > > Mprev;
    vector<vector<vector<double> > > MCprev;
    vector<vector<vector<double> > > Cprev;
    vector<vector<vector<double> > > DivCprev;
    vector<vector<vector<double> > > CMixprev;	
	
    vector<vector<vector<double> > > NTprev;
    vector<vector<vector<double> > > NIprev;
    vector<vector<vector<double> > > NDprev;
	
    vector<vector<double> > rdelprev;
    vector<vector<double> > rinsprev;
    vector<vector<double> > rtransprev;
    vector<vector<double> > svweightprev;
    vector<int> GDprev;
    vector<double> GDprobprev;
    vector<double> CSize(3);

    vector<vector<vector<double> > > M;
    vector<vector<vector<double> > > MC;
    vector<vector<vector<double> > > C;
    vector<vector<vector<double> > > DivC;
    vector<vector<vector<double> > > CMix;

    vector<vector<vector<double> > > NT;
    vector<vector<vector<double> > > NI;
    vector<vector<vector<double> > > ND;

    vector<vector<double> > rdel;
    vector<vector<double> > rins;
    vector<vector<double> > rtrans;
    vector<vector<double> > svweight;
    vector<int> GD;
    vector<double> GDprob;
  
    GD.resize(Npop);
    GDprob.resize(Npop);
    GDprev.resize(Npop);
    GDprobprev.resize(Npop);
    M.resize(Npop);
    MC.resize(Npop);
    C.resize(Npop);
    DivC.resize(Npop);
    CMix.resize(Npop);
	
    NT.resize(Npop);
    NI.resize(Npop);
    ND.resize(Npop);
	
    rdel.resize(Npop);
    rins.resize(Npop);
    rtrans.resize(Npop);
    svweight.resize(Npop);
    

    //Initialize
    for (int i = 0; i < Npop; i++) {

      C[i].resize(NChr);
      DivC[i].resize(NChr);
      M[i].resize(NChr);
      MC[i].resize(NChr);
      CMix[i].resize(NChr);
	
      NT[i].resize(NChr);
      NI[i].resize(NChr);
      ND[i].resize(NChr);
	
      for(int j = 0; j < NChr; j++){
	C[i][j].resize(1);
	CMix[i][j].resize(1);
	DivC[i][j].resize(1);
	
	NT[i][j].resize(1);
	NI[i][j].resize(1);
	ND[i][j].resize(1);
	
	M[i][j].resize(NChr);
	MC[i][j].resize(NChr);
	if(j==0){
	  C[i][j][0]=chrlength;
	  CMix[i][j][0]=chrlength;
	  CSize[j]=chrlength;
	}
	else{
	  C[i][j][0]=chrlength/(j + j*1.5);
	  CMix[i][j][0]=chrlength/(j + j*1.5);
	  CSize[j]=chrlength/(j + j*1.5);
	}
	DivC[i][j][0]=0;
    
	NT[i][j][0]=0;
	NI[i][j][0]=0;
	ND[i][j][0]=0;
    
	for(int k = 0; k < NChr; k++){
	  if(j==k){
	    M[i][j][k]=1;
	    if(j==0){
	      MC[i][j][k]=chrlength;
	    }
	    else{
	      MC[i][j][k]=chrlength/(j + j*1.5);
	    }
	  }
	  else{
	    M[i][j][k]=0;
	    MC[i][j][k]=0;
	  }
	}
      }
      svweight[i].resize(NChr);
      rdel[i].resize(NChr);
      rins[i].resize(NChr);
      rtrans[i].resize(NChr);
    
    }
	
    for(int i=0; i<Npop; i++){
    
      GDprob[i] = gnmdou*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
      GDprobprev[i] = GDprob[i];
      GD[i] = 1;
      GDprev[i] = 1;
    
      double division_rate = mu_i*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
	
      for(int j = 0; j < NChr; j++){

	double p_ins = g_d*((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) ));
	double rintensityj = mu_i_range*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
	double p_del =1-p_ins;
	svweight[i][j] = svp2;
	  
	if(j==1){
	  rdel[i][j] = 0;
	  rins[i][j] = 0;
	}
	else{
  
	  rdel[i][j] = rintensityj*p_del;
	  rins[i][j] = rintensityj*p_ins;
	
	}
	rtrans[i][j]= p_tran*((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) ));
      }
      
    }
     	
    Mprev.resize(Npop);
    MCprev.resize(Npop);
    Cprev.resize(Npop);
    DivCprev.resize(Npop);
    CMixprev.resize(Npop);
		
    NTprev.resize(Npop);
    NIprev.resize(Npop);
    NDprev.resize(Npop);
	
    rdelprev.resize(Npop);
    rinsprev.resize(Npop);
    rtransprev.resize(Npop);
    svweightprev.resize(Npop);
	
    for (int i = 0; i < Npop; i++) {

      Cprev[i].resize(NChr);
      CMixprev[i].resize(NChr);
      DivCprev[i].resize(NChr);
      Mprev[i].resize(NChr);
      MCprev[i].resize(NChr);
      rdelprev[i].resize(NChr);
      rinsprev[i].resize(NChr);
      rtransprev[i].resize(NChr);
      svweightprev[i].resize(NChr);
	
      NTprev[i].resize(NChr);
      NIprev[i].resize(NChr);
      NDprev[i].resize(NChr);
	
      for(int j = 0; j < NChr; j++){
	
	rdelprev[i][j] = rdel[i][j];
	rinsprev[i][j] = rins[i][j];
	rtransprev[i][j] = rtrans[i][j];
	svweightprev[i][j] = svweight[i][j];
	
	Cprev[i][j].resize(1);
	CMixprev[i][j].resize(1);
	DivCprev[i][j].resize(1);
	Mprev[i][j].resize(NChr);
	MCprev[i][j].resize(NChr);
	Cprev[i][j][0]=C[i][j][0];
	CMixprev[i][j][0]=CMix[i][j][0];
	DivCprev[i][j][0]=DivC[i][j][0];
    
	NTprev[i][j].resize(1);
	NIprev[i][j].resize(1);
	NDprev[i][j].resize(1);
	
	NTprev[i][j][0]=NT[i][j][0];
	NIprev[i][j][0]=NI[i][j][0];
	NDprev[i][j][0]=ND[i][j][0];
    
	for(int k = 0; k < NChr; k++){    
	  Mprev[i][j][k]=M[i][j][k];
	  MCprev[i][j][k]=MC[i][j][k];
	}
      }
    }
	

    ////////////////////////////////////////////////////////////
    ///   Loop over the generations                          ///
    ////////////////////////////////////////////////////////////
    //std::cout << "Looping over generations" << std::endl;

    for(int gg=1; gg<ngen; gg++){
      
      //////////////////////////////////////////////////////
      /////////Mutate: insertions and deletions/////////////
      //////////////////////////////////////////////////////
	
      if(gg % 100 == 0) std::cout << "particle/ngen : " << bd << " / " << gg << std::endl;

      for(int i=0; i<Npop; i++){
	//std::cout << "GDprev[i]: " << GDprev[i] << std::endl;
	double r_gd = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			
	//
	if(r_gd<GDprob[i]){
			
	  //std::cout << "oh no" << std::endl;
	  //sample a genome in the current ploidy
	  int sample_chr = ceil(GDprev[i]*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) ))-1;
				
	  vector<vector<double> > MCins;
	  vector<vector<double> > Mins;
	  vector<vector<double> > CMixins;
				
	  MCins.resize(NChr);
	  Mins.resize(NChr);
	  CMixins.resize(NChr);
				
	  for(int j=0; j<NChr; j++){
	    MCins[j].resize(NChr);
	    Mins[j].resize(NChr);
	    CMixins[j].resize(1);
				
	    for(int k=0; k<NChr; k++){
				
	      MCins[j][k]=MCprev[i][j+ sample_chr*NChr][k];
	      Mins[j][k]=Mprev[i][j+ sample_chr*NChr][k];
	      CMixins[j][0]=CMixprev[i][j+ sample_chr*NChr][0];
				
	    }
	      
	  }
			
	  MCprev[i].reserve(MCprev[i].size() + MCins.size());
	  MCprev[i].insert(MCprev[i].end(), MCins.begin(), MCins.end());
	  int a=GDprev[i]-1;
	  
	  Mprev[i].reserve(Mprev[i].size() + Mins.size());
	  Mprev[i].insert(Mprev[i].end(), Mins.begin(), Mins.end());
	    
	  CMixprev[i].reserve(CMixprev[i].size() + CMixins.size());
	  CMixprev[i].insert(CMixprev[i].end(), CMixins.begin(), CMixins.end());
	  DivCprev[i][0][0] += Cprev[i][0][sample_chr];
	  DivCprev[i][1][0] += Cprev[i][1][sample_chr];
	  DivCprev[i][2][0] += Cprev[i][2][sample_chr];
				

	  for(int m=0; m<NChr; m++){
	    NIprev[i][m][0] = NIprev[i][m][0] + NIprev[i][m][0];
	  }

	  for(int j = 0; j< NChr; j++){
	    int a=sample_chr;
	    Cprev[i][j].push_back(Cprev[i][j][a]);
	    
	  }
	  GDprev[i] += 1;
	}
	  
      }

      for(int i=0; i<Npop; i++){
		
	for(int g_d=0; g_d<GDprev[i]; g_d++){
		
	  for(int j=0; j<NChr; j++){
			
	    for(int k=0; k<NChr; k++){
				
	      double proportional_prob = ceil(Mprev[i][j + g_d*NChr][k])*Curvemax/(1+exp(-1*svgradient*(log10((MCprev[i][j + g_d*NChr][k]))-svp2)));

				
	      double p_i = rinsprev[i][k]*proportional_prob;
	      double p_d = rdelprev[i][k]*proportional_prob;		


	      double r_i = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
	      double r_d = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			
	      if(r_i < p_i){
				
		double qq = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) ); 
		double sz = SV_min - log(qq)*SV_max; 
		DivCprev[i][k][0] += sz;
		Cprev[i][k][g_d] = Cprev[i][k][g_d] + sz;
		MCprev[i][j + g_d*NChr][k] = MCprev[i][j + g_d*NChr][k] + sz;
		NIprev[i][k][0] = NIprev[i][k][0] + 1;

		  
	      }
				
	      if(r_d < p_d){

		double qq = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) ); 
		double sz = SV_min - log(qq)*SV_max; 
		DivCprev[i][k][0] -= sz;
		Cprev[i][k][g_d] = Cprev[i][k][g_d] - sz;
		MCprev[i][j + g_d*NChr][k] = MCprev[i][j + g_d*NChr][k] - sz;
		NDprev[i][k][0] = NDprev[i][k][0] + 1;

	      }
		
	    }
			
			
	    
	  }
		
	}
      }
		
		
		
      /////////////////////////////
      /////////UpdateM/////////////
      /////////////////////////////
      //std::cout << "\tUpdating M" << std::endl;
      for(int i=0; i<Npop; i++){
		
	for(int g_d=0; g_d<GDprev[i]; g_d++){
		
	  for(int j=0; j<NChr; j++){
			
	    double sum_j = 0;
			
	    for(int k=0; k<NChr; k++){
	      
	      Mprev[i][j + g_d*NChr][k] = MCprev[i][j + g_d*NChr][k]/Cprev[i][k][g_d];
				
	      sum_j = Mprev[i][j + g_d*NChr][k]*Cprev[i][k][g_d] + sum_j;
							
	    }
			
	    CMixprev[i][j + g_d*NChr][0]=sum_j;

	  }
	  
	}
		
      }


      //////////////////////////////////////////////////////
      /////////Mutate: translocations///////////////////////
      //////////////////////////////////////////////////////
      //std::cout << "\tTranslocations" << std::endl;
      for(int i=0; i<Npop; i++){
	
	int jtick = 0;
	int jcounte = 0;
	for(int j=0; j<NChr*GDprev[i]; j++){
		
	  if(jtick==NChr){
	    jcounte = jcounte + 1;
	    jtick = 0;
	  }
			
	  for(int k=0; k<NChr; k++){
			
	    //translocate between C1 and C2//
	    double proportional_prob = Mprev[i][j][k]*Curvemax/(1+exp(-1*svtransgrad*(log10(MCprev[i][j][k])-trp2)));
	    double p1t = proportional_prob*rtransprev[i][k];
	    double p1tp = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			
			
	    if( p1tp < p1t ){
			
	      ///Mutate: second chomsome///
	      int j1tick = 0;
	      int j1counte = 0;
	      for(int j1=0; j1<NChr*GDprev[i]; j1++){
		  
		if(j1tick==NChr){
		  j1counte = j1counte + 1;
		  j1tick = 0;
		}
					
		if(j!=j1){
		    
		  for(int k1=0; k1<NChr; k1++){
					
		    if(k!=k1){
			
		      double proportional_probj = Mprev[i][j1][k1]*Curvemax/(1+exp(-1*svtransgrad*(log10(MCprev[i][j1][k1])-trp2)));
		      double p1tj = proportional_probj*rtransprev[i][k1];
		      double p1tpj = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
						
		      if( p1tpj < p1tj ){
			
			double sum_j0 = 0;
			double sum_j1 = 0;
			
			for(int ka=0; ka<NChr; ka++){

			  double p = lowB + (1-lowB)*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			  double a1 = Mprev[i][j][ka];
			  double b1 = Mprev[i][j1][ka];
			  Mprev[i][j][ka] = p*(a1+b1);
			  Mprev[i][j1][ka] = (1-p)*(a1+b1);
			  sum_j0 += Mprev[i][j][ka]*Cprev[i][ka][jcounte];
			  sum_j1 += Mprev[i][j1][ka]*Cprev[i][ka][j1counte];
				
			
			}
			
			CMixprev[i][j][0]=sum_j0;
			CMixprev[i][j1][0]=sum_j1;
			
			NTprev[i][k][0] = NTprev[i][k][0] + 1;
			NTprev[i][k1][0] = NTprev[i][k1][0] + 1;
			
		      }
			
		    }
			
			
		  }
					
		}
			
			
		j1tick = j1tick + 1;
			
	      }
			
	    }
			
	  }
			
	  jtick = jtick + 1;
			
	}
			

      }
		
	
      for(int i=0; i<Npop; i++){
	
	for(int g_d=0; g_d<GDprev[i]; g_d++){
		
	  for(int j=0; j<NChr; j++){
			
	    double sum_j = 0;
			
	    for(int k=0; k<NChr; k++){
				
	      MCprev[i][j + g_d*NChr][k] = Mprev[i][j + g_d*NChr][k]*Cprev[i][k][g_d];
				
	      sum_j = Mprev[i][j + g_d*NChr][k]*Cprev[i][k][g_d] + sum_j;
							
	    }
			
	    CMixprev[i][j + g_d*NChr][0]=sum_j;
	    
	  }
	  
	}
	
      }


      //////////////////////////////////////////////////////
      /////////Selection////////////////////////////////////
      //////////////////////////////////////////////////////
      //std::cout << "\tSelection" << std::endl;
      vector<int> keep;
	
      keep.resize(Npop);
	
      int n_remaining = 0;
	
      for(int i=0; i<Npop; i++){
	
	int check_len = 0;
		
	keep[i] = 0;
		
	for(int g_d=0; g_d<GDprev[i]; g_d++){
			
	  for(int j=0; j<NChr; j++){
			
	    double p1 = 1 - (Curvemax/(1+exp(-1*threshold_g*((CMixprev[i][j + g_d*NChr][0])-minchr))));
	    double p2 = (Curvemax/(1+exp(-1*threshold_g*((CMixprev[i][j + g_d*NChr][0])-maxchr))));
	      
	    double r1 =(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
	    double r2 =(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
	    
	    if(r1<p1){
			
	      check_len = 1;
	    }
	      
	    if(r2<p2){
	      check_len = 1;
	      
	    }
			
				
	  }
			
	}
	  
	if(check_len==0){
		
	  keep[i] = 1;
	  n_remaining += 1;
	}
	
      }
	
      vector<int> sample_vec;
      sample_vec.resize(n_remaining);
      int counter = 0;
      for(int q=0; q<Npop; q++){
	
	if(keep[q]==1){
	
	  sample_vec[counter] = q;
	  counter += 1; 
	
	}
	
      }

      //////////////////////////////////////////////////////
      /////////Resample/////////////////////////////////////
      //////////////////////////////////////////////////////
      //std::cout << "\tResample" << std::endl;
      for (int i = 0; i < Npop; i++) {
	int newgd = 0;
	
	if(keep[i]==1){
	
	  if(GD[i]!=GDprev[i]){

	    MC[i].resize(MCprev[i].size());
	    M[i].resize(Mprev[i].size());
	    CMix[i].resize(CMixprev[i].size());

	    for(int j=0; j<M[i].size(); j++){
	
	      MC[i][j].resize(NChr);
	      M[i][j].resize(NChr);
	      CMix[i][j].resize(1);
	
	    }
	
	    for(int j = 0; j< NChr; j++){
				
	      C[i][j].resize(Cprev[i][j].size());
	      
	    }
	  }
	
	  for(int g_d=0; g_d<GDprev[i]; g_d++){
	
	    for(int j = 0; j < NChr; j++){

	      C[i][j][g_d]=Cprev[i][j][g_d];
	      DivC[i][j][0]=DivCprev[i][j][0];
	      CMix[i][j + g_d*NChr][0]=CMixprev[i][j + g_d*NChr][0];
    
	      rdel[i][j] = rdelprev[i][j];
	      rins[i][j] = rinsprev[i][j];
	      rtrans[i][j] = rtransprev[i][j];
	      svweight[i][j] = svweightprev[i][j];
    
	      NT[i][j][0] = NTprev[i][j][0];
	      NI[i][j][0] = NIprev[i][j][0];
	      ND[i][j][0] = NDprev[i][j][0];
    
	      for(int k = 0; k < NChr; k++){   
     
		M[i][j + g_d*NChr][k]=Mprev[i][j + g_d*NChr][k];
		MC[i][j + g_d*NChr][k]=MCprev[i][j + g_d*NChr][k];
		  
	      }
	
	    }
	  }
	  newgd = GDprev[i];
	  GDprob[i] = GDprobprev[i];
	
	}
	else{

	  double rnew = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
	  int position = floor(rnew*n_remaining);
	
	  if(GD[i]!=GDprev[sample_vec[position]]){
	
	    MC[i].resize(MCprev[sample_vec[position]].size());
	    M[i].resize(Mprev[sample_vec[position]].size());
	    CMix[i].resize(CMixprev[sample_vec[position]].size());
	
	    for(int j=0; j<M[i].size(); j++){
	
	      MC[i][j].resize(NChr);
	      M[i][j].resize(NChr);
	      CMix[i][j].resize(1);
	
	    }
	
	    for(int j = 0; j< NChr; j++){
	
				
	      C[i][j].resize(Cprev[sample_vec[position]][j].size());
			
	    }

	
	  }
	
	  for(int g_d=0; g_d<GDprev[sample_vec[position]]; g_d++){
	
	    for(int j = 0; j < NChr; j++){

	      C[i][j][g_d]=Cprev[sample_vec[position]][j][g_d];
	      DivC[i][j][0]=DivCprev[sample_vec[position]][j][0];
	      CMix[i][j + g_d*NChr][0]=CMixprev[sample_vec[position]][j + g_d*NChr][0];
    
	      rdel[i][j] = rdelprev[sample_vec[position]][j];
	      rins[i][j] = rinsprev[sample_vec[position]][j];
	      rtrans[i][j] = rtransprev[sample_vec[position]][j];
	      svweight[i][j] = svweightprev[sample_vec[position]][j];
    
	      NT[i][j][0] = NTprev[sample_vec[position]][j][0];
	      NI[i][j][0] = NIprev[sample_vec[position]][j][0];
	      ND[i][j][0] = NDprev[sample_vec[position]][j][0];
    
	      for(int k = 0; k < NChr; k++){   
     
		M[i][j + g_d*NChr][k]=Mprev[sample_vec[position]][j + g_d*NChr][k];
		MC[i][j + g_d*NChr][k]=MCprev[sample_vec[position]][j + g_d*NChr][k];
		
	      }
	
	
	    }

	  }


	  newgd = GDprev[sample_vec[position]];
	  GDprob[i] = GDprobprev[sample_vec[position]];

	}
	
	GD[i] = newgd;
	
      }
	
      //////////////////////////////////////////////////////
      /////////Update vectors///////////////////////////////
      //////////////////////////////////////////////////////

      for (int i = 0; i < Npop; i++) {
	
	if(GDprev[i]!=GD[i]){

	  MCprev[i].resize(MC[i].size());
				
	  Mprev[i].resize(M[i].size());
			
	  CMixprev[i].resize(CMix[i].size());
	
	  for(int j=0; j<Mprev[i].size(); j++){
	
	    MCprev[i][j].resize(NChr);
	    Mprev[i][j].resize(NChr);
	    CMixprev[i][j].resize(1);
	
	  }
	
	  for(int j = 0; j< NChr; j++){
	
				
	    Cprev[i][j].resize(C[i][j].size());
				
	  }
			
			
	
	}
	
	for(int g_d=0; g_d<GD[i]; g_d++){

	  for(int j = 0; j < NChr; j++){	
	
	    Cprev[i][j][g_d]=C[i][j][g_d];
   
	    for(int k = 0; k < NChr; k++){   
     
	      Mprev[i][j + g_d*NChr][k]=M[i][j + g_d*NChr][k];
	      MCprev[i][j + g_d*NChr][k]=MC[i][j + g_d*NChr][k];
		
    
	    }
	    DivCprev[i][j][0]=DivC[i][j][0];
	    CMixprev[i][j + g_d*NChr][0]=CMix[i][j + g_d*NChr][0];
    
	    rdelprev[i][j] = rdel[i][j];
	    rinsprev[i][j] = rins[i][j];
	    rtransprev[i][j] = rtrans[i][j];
	    svweightprev[i][j] = svweight[i][j]; 
    
	    NTprev[i][j][0] = NT[i][j][0];
	    NIprev[i][j][0] = NI[i][j][0];
	    NDprev[i][j][0] = ND[i][j][0];
    
    
	
	
	  }
	
	}
	GDprev[i] = GD[i];
      }
	
    } // Loop over generations

      
    //////////////////////////////////////////////////////
    /////////Print out simulation results ////////////////
    //////////////////////////////////////////////////////
      
    char Datins[100];
    int Dnins;
    Dnins=sprintf(Datins,"results-S%d-P%d-M%d.dat",InternalN,bd,RpNumber);
    fstream myfileIns; 
    myfileIns.open(Datins,ios::out);
 	
    if(myfileIns.is_open()){
 	
      myfileIns << "Num" << ", " << "CNum" <<  ", " << "CSize" ", " << "PDiv" << ", " << "CDiv" << ", " << "RelLen" ", " << "NTrans" << ", " << "NIns" << ", " << "NDel" << ", " << "Rins" << ", " << "Rdel" << ", " << "Rtrans" << '\n';
      
      for (int i = 0; i < Npop; i++) {
	for(int j = 0; j < NChr; j++){

	  
	  double log_signed_div = 0;
	
	  if(DivCprev[i][j][0]==0){
	    log_signed_div = 0;
	  }
	  else if(DivCprev[i][j][0]<0){
	    log_signed_div = -1*log10(abs(DivCprev[i][j][0]));
	
	  }
	  else{
	    log_signed_div = log10(DivCprev[i][j][0]);
	  }
    
	  myfileIns << i << ", "  << j << ", " << CSize[j] << ", " << log_signed_div/CSize[j] << ", " << log_signed_div << ", " << DivCprev[i][j][0]/CSize[j] << ", " << NT[i][j][0] << ", " << NI[i][j][0] << ", " << ND[i][j][0] << ", " << rins[i][j] << ", " << rdel[i][j] << ", " << rtrans[i][j] << '\n';
	 
	}
      }
      myfileIns.close();

    }
     else{
       std::cerr << "error: open of results file unsuccessful: " << Datins << std::endl;
      return(1);
    }

    //////////////////////////////////////////////////////
    /////////Print out parameters       //////////////////
    //////////////////////////////////////////////////////
    char Dat[100];
    int Dn;
    Dn=sprintf(Dat,"parameters-S%d-P%d-M%d.dat",InternalN,bd,RpNumber);

    fstream myfileDatParam; 
    myfileDatParam.open(Dat,ios::out);

    if(myfileDatParam.is_open()){
	
      myfileDatParam << Npop << '\t';
      myfileDatParam << ngen << '\t';
      myfileDatParam << g_d << '\t';
      myfileDatParam << mu_i << '\t';
      myfileDatParam << p_tran << '\t';
      myfileDatParam << svp2 << '\t';
      myfileDatParam << trp2 << '\t';
      myfileDatParam << SV_min << '\t';
      myfileDatParam << svgradient << '\t';
      myfileDatParam << svtransgrad << '\t';
      myfileDatParam << threshold_g << '\t';
      myfileDatParam << SV_max << '\t';
      myfileDatParam << mu_i_range << '\t';
      myfileDatParam << gnmdou << '\t';
      myfileDatParam << maxchr << '\t';
      myfileDatParam << minchr << std::endl;
	
      myfileDatParam.close();
    }
    else{
      std::cerr << "error: open of parameter file unsuccessful: " << Dat << std::endl;
      return(1);
    }

  } // end loop over particles
  
  return(0); 
}


















