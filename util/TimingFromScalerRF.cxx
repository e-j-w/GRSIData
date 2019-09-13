// g++ TimingFromScalerRF.cxx -I$GRSISYS/include -I$GRSISYS/GRSIData/include -L$GRSISYS/libraries `root-config --cflags --libs` -lTreePlayer -o TimingFromScalerRF

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

#include "TChannel.h"

#include "TScalerFrag.h"
#include "TTigress.h"

std::vector<float> timestampBuffer, phaseBuffer;

void BuildRFTimeTable(TTree* stree)
{

   TScalerFrag* currentFrag = nullptr;
   TBranch* branch = stree->GetBranch("TScalerFrag");
   branch->SetAddress(&currentFrag);

   int entries = stree->GetEntries();

   int FragsIn = 0;

   for(int i=0;i<entries;i++){

      stree->GetEntry(i);
      FragsIn++;
      
      for(unsigned int j=0;j<currentFrag->fName.size();j++){
         if(j<currentFrag->fData.size()){
            if(strcmp(currentFrag->fName[j].c_str(),"RF par 1")==0){ //phase appears to be par 1 or 2
               printf("sn=%i, ts=%li, phase=%f\n",currentFrag->fDaqId,currentFrag->fDaqTimeStamp,currentFrag->fData[j]);
               phaseBuffer.push_back(currentFrag->fData[j]);
               timestampBuffer.push_back(currentFrag->fDaqTimeStamp);
               break;
            }
         }
      }

   }

   printf("%i RF fragment(s) read in.\n", FragsIn);

   return;
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
         

      }

   }
   
   return;
}

#ifndef __CINT__

int main(int argc, char** argv)
{

   if(argc != 4) {
      printf("try again (usage: %s <fragment tree file> <input analysis tree file> <output analysis tree file>).\n", argv[0]);
      printf("\nThis code takes in a fragment tree with RF scaler data and writes RF corrected timing values to a corresponding analysis tree.\n");
      return 0;
   }

   auto* ffile = new TFile(argv[1]);
   if((ffile == nullptr) || (!ffile->IsOpen())) {
      printf("Failed to open file '%s'!\n", argv[1]);
      return 1;
   }

   auto* inpfile = new TFile(argv[2]);
   if((inpfile == nullptr) || (!inpfile->IsOpen())) {
      printf("Failed to open file '%s'!\n", argv[2]);
      return 1;
   }

   auto* outfile = new TFile(argv[3],"RECREATE");
   if((outfile == nullptr) || (!outfile->IsOpen())) {
      printf("Failed to open file '%s'!\n", argv[3]);
      return 1;
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

   TTree* outtree = new TTree();

   //initialize arrays
   timestampBuffer = {};
   phaseBuffer = {};
   
   BuildRFTimeTable(scalertree);

   WriteCorrectedTimes(inptree,outtree);

   //auto* outfile = new TFile(argv[3], "recreate");
   //list->Write();

   ffile->Close();
   inpfile->Close();
   outfile->Close();

   return 0;
}

#endif
