// g++ TimingFromScalerRF.cxx -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lTTigress -lTEmma -lTGRSIDetector -lTGenericDetector -o TimingFromScalerRF

// Takes in a fragment tree with RF scaler data and writes RF corrected timing values to a corresponding analysis tree
//
//
#include <utility>
#include <vector>
#include <cstdio>
#include <iostream>
#include <iomanip>

#include "TTree.h"
#include "TTreeIndex.h"
#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TApplication.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"

#include "TChannel.h"

#include "TScalerFrag.h"
#include "TTigress.h"
#include "TEmma.h"

#include <cstddef>

TH1D *periodHist, *phaseHist, *parHist;

std::vector<double> phaseBuffer, freqBuffer, parBuffer[4];
std::vector<ULong64_t> timestampBuffer;

char str[256];
char ftreename[30][256], atreename[30][256];

int BuildRFTimeTable(TTree* stree)
{

   TScalerFrag* currentFrag = nullptr;
   TBranch* branch = stree->GetBranch("TScalerFrag");
   branch->SetAddress(&currentFrag);

   int entries = stree->GetEntries();

   int FragsIn = 0;
   int RFFragsIn = 0;
   bool fe,p0e,p1e,p2e,de,te;//flags to check whether specific RF data exists in a fragment
   
   //reserve enough space for long runs
   freqBuffer.reserve(entries);
   parBuffer[0].reserve(entries);
   parBuffer[1].reserve(entries);
   parBuffer[2].reserve(entries);
   parBuffer[3].reserve(entries);
   timestampBuffer.reserve(entries);

   for(int i=0;i<entries;i++){

      stree->GetEntry(i);
      FragsIn++;

      fe=false;
      p0e=false;
      p1e=false;
      p2e=false;
      de=false;
      double val;

      timestampBuffer.push_back(currentFrag->GetTimeStamp());

      
      
      for(int j=0;j<currentFrag->fName.size();j++){
         val=currentFrag->GetData(j);
         //printf("%s: %f\n",currentFrag->fName[j].c_str(),val);
         if(j<currentFrag->fData.size()){
            if(strcmp(currentFrag->fName[j].c_str(),"RF Frequency")==0){ 
               freqBuffer.push_back(val);
               fe=true;
            }else if(strcmp(currentFrag->fName[j].c_str(),"RF par 0")==0){
               parBuffer[0].push_back(val);
               p0e=true;
            }else if(strcmp(currentFrag->fName[j].c_str(),"RF par 1")==0){
               parBuffer[1].push_back(val);
               p1e=true;
            }else if(strcmp(currentFrag->fName[j].c_str(),"RF par 2")==0){
               parBuffer[2].push_back(val);
               p2e=true;
            }else if(strcmp(currentFrag->fName[j].c_str(),"RF Determinant")==0){
               parBuffer[3].push_back(val);
               de=true;
            }
         }
      }


      if(fe && p0e && p1e && p2e && de){
         RFFragsIn++;
      }else{
         //make all vectors the same size and then remove the last element of each (corresponding to the bad fragment)
         if(!fe){
            freqBuffer.push_back(0.0);
         }else if(!p0e){
            parBuffer[0].push_back(0.0);
         }else if(!p1e){
            parBuffer[1].push_back(0.0);
         }else if(!p2e){
            parBuffer[2].push_back(0.0);
         }else if(!de){
            parBuffer[3].push_back(0.0);
         }else{
            printf("RF fragment error\n");
            exit(-1);
         }
         freqBuffer.pop_back();
         parBuffer[0].pop_back();
         parBuffer[1].pop_back();
         parBuffer[2].pop_back();
         parBuffer[3].pop_back();
         timestampBuffer.pop_back();
      }

   }

   printf("%i RF fragment(s) read in.\n", FragsIn);

   double t0,T,A,s,c;
   double par[3];
   
   for(int i=0;i<RFFragsIn;i++){
      
      if(parBuffer[3][i]!=0.){
         par[0]=parBuffer[0][i]/parBuffer[3][i];
         par[1]=parBuffer[1][i]/parBuffer[3][i];
         par[2]=parBuffer[2][i]/parBuffer[3][i];
         T = 1.34217728E9/freqBuffer[i]; //period in ns
         //printf("parBuffer[0][i]: %f, parBuffer[1][i]: %f, parBuffer[2][i]: %f, parBuffer[3][i]: %f\n",parBuffer[0][i],parBuffer[1][i],parBuffer[2][i],parBuffer[3][i]);
         //printf("par0: %f, par1: %f, par2: %f, T: %f\n",par[0],par[1],par[2],T);
         //getc(stdin);
         A = sqrt(par[1]*par[1] + par[0]*par[0]);
         s = -par[0] / A;
         c = par[1] / A;
         if(s >= 0) {
            t0 = acos(c) * T / (2 * TMath::Pi());
         } else {
            t0 = (1 - acos(c) / (2 * TMath::Pi())) * T;
         }
         phaseBuffer.push_back(t0); //in ns
         phaseHist->Fill(t0);
         periodHist->Fill(T); //period in ns
         parHist->Fill(A);
      }

      printf("ts=%llu, phase=%f, T=%f\n",timestampBuffer[i],t0,T);
      getc(stdin);
      
   }
   

   return RFFragsIn;
}

void PhaseInterpTest(int numRFFrags){

   TH1D *interpHist = new TH1D("interpolated - actual phase shift","interpolated - actual phase shift",20000,-100,100);
   interpHist->Reset();

   double phase, interpPhaseShift, interpPhase;
   double T, numTinterp;

   for(int i=1;i<numRFFrags;i++){
      phase=phaseBuffer[i];
      T = 1.34217728E9/freqBuffer[i]; //period in ns
      if(timestampBuffer[i] > timestampBuffer[i-1])
         numTinterp = (timestampBuffer[i]-timestampBuffer[i-1])/T;
      else //handle rollover of timestamps
         numTinterp = ((timestampBuffer[i]+43980465111040)-timestampBuffer[i-1])/T; //timestamp is a 42-bit unsigned int
      interpPhaseShift = (numTinterp - (int)numTinterp)*T;
      if(phaseBuffer[i-1] + interpPhaseShift < T)
         interpPhase = phaseBuffer[i-1] + interpPhaseShift;
      else
         interpPhase = interpPhaseShift - (T - phaseBuffer[i-1]); 
      //printf("phase: %f, numTinterp: %f, interpPhase: %f\n",phase,numTinterp,interpPhase);
      //getc(stdin);
      interpHist->Fill(interpPhase - phase);
   }

   TCanvas* c1 = new TCanvas("c1","c1",800,600);
   interpHist->GetXaxis()->SetTitle("interpolated - actual phase shift (ns)");
   interpHist->Draw();

}

double GetPhaseNsForTimestamp(ULong64_t timestamp, int numRFFrags){
   double T, numTinterp, interpPhaseShift;
   double interpPhase = -100;
   for(int i=0;i<numRFFrags;i++){
      if((i>0)&&(timestampBuffer[i]>timestamp)&&(timestampBuffer[i-1]<=timestamp)){ //ignoring rollover for now
         T = 1.34217728E9/freqBuffer[i-1]; //period in ns
         numTinterp = (timestamp-timestampBuffer[i-1])/T;
         interpPhaseShift = (numTinterp - (int)numTinterp)*T;
         if(phaseBuffer[i-1] + interpPhaseShift < T)
            interpPhase = phaseBuffer[i-1] + interpPhaseShift;
         else
            interpPhase = interpPhaseShift - (T - phaseBuffer[i-1]);
         break;
      }
   }
   return interpPhase;
}

void MapPhaseTest(TTree* atree, int numRFFrags, TH1D *interpHist, bool append){
   
   if(!append)
      interpHist->Reset();

   ULong64_t ts;

   //TEmma* currentFrag = nullptr;
   //TBranch* branch = atree->GetBranch("TEmma");
   TTigress* currentFrag = nullptr;
   TBranch* branch = atree->GetBranch("TTigress");
   branch->SetAddress(&currentFrag);

   int entries = atree->GetEntries();

   int FragsIn = 0;

   for(int i=0;i<entries;i++){

      if(i%100==0)
         printf("Entry %i / %i\r",i,entries);

      atree->GetEntry(i);
      FragsIn++;    

      //for(int j=0;j<currentFrag->GetICMultiplicity();j++){
      for(int j=0;j<currentFrag->GetMultiplicity();j++){
         //ts = currentFrag->GetICHit(j)->GetTimeStampNs();
         double e = currentFrag->GetTigressHit(j)->GetEnergy()/1000.0;
         //printf("energy: %f\n",e);
         if(e>0.){
            ts = currentFrag->GetTigressHit(j)->GetTimeStampNs();
            //printf("tigress ts: %lu\n",ts);
            interpHist->Fill(GetPhaseNsForTimestamp(ts,numRFFrags));
         }
         
         //printf("emma ts: %llu\n",ts);
      }
   }

}

void MapPhase2DTest(TTree* atree, int numRFFrags){

   TH1D *interpHist = new TH1D("TIGRESS-EMMA timing","TIGRESS-EMMA timing",4000,-10000,10000);
   interpHist->Reset();

   ULong64_t ts1,ts2;
   double modts1,modts2;

   TEmma* currentEmmaFrag = nullptr;
   TTigress* currentTigFrag = nullptr;
   TBranch* branchEmma = atree->GetBranch("TEmma");
   branchEmma->SetAddress(&currentEmmaFrag);
   TBranch* branchTig = atree->GetBranch("TTigress");
   branchTig->SetAddress(&currentTigFrag);

   int entries = atree->GetEntries();

   for(int i=0;i<entries;i++){

      if(i%100==0)
         printf("Entry %i / %i\r",i,entries);

      atree->GetEntry(i);

      for(int j=0;j<currentTigFrag->GetMultiplicity();j++){
         for(int k=0;k<currentEmmaFrag->GetICMultiplicity();k++){
            ts1 = currentTigFrag->GetTigressHit(j)->GetTimeStampNs();
            ts2 = currentEmmaFrag->GetICHit(k)->GetTimeStampNs();
            modts1 = static_cast<double>(ts1) + GetPhaseNsForTimestamp(ts1,numRFFrags);
            modts2 = static_cast<double>(ts2) + GetPhaseNsForTimestamp(ts2,numRFFrags);
            //modts1 = fmod(ts1 + GetPhaseNsForTimestamp(ts1,numRFFrags),84.841646);
            //modts2 = fmod(ts2 + GetPhaseNsForTimestamp(ts2,numRFFrags),84.841646);
            //printf("filling at: %f %f\n",modts2-modts1,modts2);
            //printf("diffs: %f %f\n",modts1-ts1,modts2-ts2);
            //printf("ts1: %lu, ts2: %lu, modts1: %f, modts2: %f, diffts: %f, diffmodts: %f\n",ts1,ts2,modts1,modts2,static_cast<double>(ts2)-static_cast<double>(ts1),modts2-modts1);
            interpHist->Fill(static_cast<double>(ts2)-static_cast<double>(ts1));
            //interpHist->Fill(modts2-modts1);
         }
      }
      
   }

   TCanvas* c1 = new TCanvas("c1","c1",800,600);
   interpHist->GetXaxis()->SetTitle("t_{emma} - t_{tigress} (ns)");
   interpHist->Draw();

}

void WriteCorrectedTimes(TTree* inptree, TTree* outtree)
{
   
   TBranch* branch = inptree->GetBranch("TTigress");
   if(branch == nullptr) {
      printf("Tigress branch not in analysis tree.\n");
   } else {
      printf("Tigress branch in analysis tree.\n");
      TTigress* currentDet = nullptr;
      branch->SetAddress(&currentDet);

      int entries = inptree->GetEntries();

      int FragsIn = 0;

      for(int i=0;i<entries;i++){

         inptree->GetEntry(i);
         FragsIn++;
         printf("Fragment %i\n",FragsIn);
         for(int j=0;j<currentDet->GetHitVector().size();j++){
            printf("ts=%li\n",currentDet->GetHitVector()[j]->GetTimeStampNs());
         }
         getc(stdin);

      }

   }
   
   return;
}

#ifndef __CINT__

int main(int argc, char** argv)
{
   TApplication *theApp;

   FILE *list, *list2;
   TFile* ffile;
   TFile* inpfile;

   if(argc != 3) {
      printf("try again (usage: %s <fragment tree list> <input analysis tree list>).\n", argv[0]);
      printf("\nThis code takes in a fragment tree with RF scaler data and writes RF corrected timing values to a corresponding analysis tree.\n");
      return 0;
   }

   TH1D *interpHist = new TH1D("TIGRESS mapped phase shift","TIGRESS mapped phase shift",20000,-101,100);
   interpHist->Reset();

   //read in tree list file
   if((list=fopen(argv[1],"r"))==NULL){
      printf("ERROR: Cannot open the fragment tree list file: %s\n",argv[1]);
      exit(-1);
   }
   int ind=0;
   while(fscanf(list,"%s",str)!=EOF){
      strcpy(ftreename[ind],str);
      ind++;
   }

   if((list2=fopen(argv[2],"r"))==NULL){
      printf("ERROR: Cannot open the analysis tree list file: %s\n",argv[1]);
      exit(-1);
   }
   ind=0;
   while(fscanf(list2,"%s",str)!=EOF){
      strcpy(atreename[ind],str);
      ind++;
   }

   
   //scan the list files for ROOT files
   for(int i=0;i<ind;i++){

      ffile = new TFile(ftreename[i]);
      if((ffile == nullptr) || (!ffile->IsOpen())) {
         printf("Failed to open file '%s'!\n", argv[1]);
         exit(-1);
      }

      inpfile = new TFile(atreename[i]);
      if((inpfile == nullptr) || (!inpfile->IsOpen())) {
         printf("Failed to open file '%s'!\n", argv[2]);
         exit(-1);
      }

      TTree* scalertree = dynamic_cast<TTree*>(ffile->Get("ScalerTree"));
      if(scalertree == nullptr) {
         printf("Failed to find scaler tree in file: %s\n", argv[1]);
         exit(-1);
      } else {
         std::cout<<scalertree->GetEntries()<<" scaler entries"<<std::endl;
      }

      TTree* inptree = dynamic_cast<TTree*>(inpfile->Get("AnalysisTree"));
      if(inptree == nullptr) {
         printf("Failed to find analysis tree in file: %s\n", argv[2]);
         exit(-1);
      } else {
         std::cout<<inptree->GetEntries()<<" entries"<<std::endl;
      }

      //initialize arrays
      timestampBuffer = {};
      phaseBuffer = {};
      freqBuffer = {};
      for(int i=0;i<4;i++)
         parBuffer[i] = {};

      periodHist = new TH1D("period","period",20000,84.83,84.85);
      periodHist->Reset();
      phaseHist = new TH1D("phase","phase",20000,0,200);
      phaseHist->Reset();
      parHist = new TH1D("parameter 0 1 sumsquare","parameter 0 1 sumsquare",20000,-20,20);
      parHist->Reset();
      

      int numFrags = BuildRFTimeTable(scalertree);
      //PhaseInterpTest(numFrags);

      

      MapPhaseTest(inptree,numFrags,interpHist,true);
      //MapPhase2DTest(inptree,numFrags);

   }
   
   //TCanvas* c1 = new TCanvas("c1","c1",800,600);
   //periodHist->Draw();
   //phaseHist->Draw();
   //parHist->Draw();
   theApp=new TApplication("App", &argc, argv);
   TCanvas* c1 = new TCanvas("c1","c1",800,600);
   interpHist->GetXaxis()->SetTitle("phase shift (ns)");
   interpHist->Draw();
   theApp->Run(kTRUE);
   

   //WriteCorrectedTimes(inptree,outtree);

   //auto* outfile = new TFile(argv[3], "recreate");
   //list->Write();

   ffile->Close();
   inpfile->Close();

   return 0;
}

#endif
