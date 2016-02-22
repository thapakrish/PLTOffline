////////////////////////////////////////////////////////////////////
//
//
// Krishna Thapa <kthapa@cern.ch>
// Tue Feb 16 2015
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include "PLTAlignment.h"
#include "PLTEvent.h"
#include "PLTU.h"
#include "PLTTrack.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"

/*
  For all tracks that touch three planes of a telescope, keep track of Slopes, Residuals, and Beamspot position. Categorize the BX's by looking at the number of tracks/events. For now, it is set up arbitrarily as follows:

  1) Colliding: Any BX with at least 10% of the count of BX_max:Non_Colliding_Threshold(NCT)
  2) Non-Colliding: Any BX that is non-colliding, but above "empty", where empty = NCT / 2  
  3) Rest: < NCT / 2
  
  Non-colliding bunches' count would be bigger than the average of the "rest".  
 */

int BXTracks (std::string const DataFileName, std::string const GainCalFileName, std::string const AlignmentFileName)
{
  std::cout << "DataFileName:    " << DataFileName << std::endl;
  std::cout << "GainCalFileName:    " << GainCalFileName << std::endl;

  PLTU::SetStyle();
  TFile *f = new TFile("BXTracks.root","RECREATE");


  // Grab the plt event reader
  //****************************************************************************************
  PLTEvent Event(DataFileName, GainCalFileName, AlignmentFileName);
  Event.SetPlaneFiducialRegion(PLTPlane::kFiducialRegion_All);
  Event.SetPlaneClustering(PLTPlane::kClustering_AllTouching, PLTPlane::kFiducialRegion_All);

  gStyle->SetOptStat(1111);

  Int_t nBX;
  Float_t BsX, BsY, Sx, Sy, Rx, Ry;
  Float_t TxY, TxZ,TyX,TyZ;

  const int nbunches = 3565;

  TH1F *bunchTriggered = new TH1F ("bunchTriggered","bunchTriggered",nbunches,-0.5,3564.5);

  
  TH1F *nTracks = new TH1F ("nTracks","#Tracks",nbunches,-0.5,3564.5);
  TH1F *tCoinc = new TH1F ("tCoinc","#Triple Coincidences",nbunches,-0.5,3564.5);
  TH1F *aHits = new TH1F ("aHits","Hits for all Telescopes",nbunches,-0.5,3564.5);
  TH1F *tHits = new TH1F ("tHits","Hits @ Triple Coinc.",nbunches,-0.5,3564.5);          
  
  std::vector<std::vector<std::pair<Float_t,Float_t> > >  Slope_SxSy(nbunches);
  std::vector<std::vector<std::pair<Float_t,Float_t> > >  BeamSpot_BxBy(nbunches);

  std::vector<std::vector<std::pair<Float_t,Float_t> > >  Zpos_XzYz(nbunches);
  std::vector<std::vector<std::pair<Float_t,Float_t> > >  Residuals_RxRy(nbunches);    
  

  unsigned int begin = 10000;
  unsigned int end = 30000;  

  std::vector<int> cb(3564); // to get the maximum count among all BX's

  //start event looping
  for (int ientry = 0; Event.GetNextEvent() >= 0; ++ientry) {
    nBX = (Event.BX() - 24 ) % 3564;
    //    time = Event.Time();
    //    std::cout << "Start time: " << time << std::endl;
    
    if (ientry % 10000 == 0) {
      std::cout << "Processing entry: " << ientry << std::endl;
    }
    
    //    if (ientry >= end){break;}
    //    if (ientry < begin){continue;}    
    
    //    if (ientry >= begin){
      
      bunchTriggered->Fill(nBX);
      
      for (size_t it = 0; it != Event.NTelescopes(); ++it) {

        // THIS telescope is
        PLTTelescope* Telescope = Event.Telescope(it);

        // all hits
        for (size_t hits = 0; hits != Telescope->NHits(); ++hits) {
          aHits->Fill(nBX);       
        }
        //        check
        //        std::cout << Telescope->NHits() << " " << Telescope->Plane(0)->NHits()+ Telescope->Plane(1)->NHits()+Telescope->Plane(2)->NHits()<< std::endl;
        
        int phit = Telescope->HitPlaneBits();

            
        // If hits in all three planes        
        if(phit==0x7){
        
          tCoinc->Fill(nBX);

          for (size_t hits = 0; hits != Telescope->NHits(); ++hits) {
            tHits->Fill(nBX);     
          }
          
          cb[nBX]++;        
          for (uint trk = 0; trk != Telescope->NTracks(); ++trk) {

            // THIS tracks is
            PLTTrack* Track = Telescope->Track(trk);
            
            Sx = (Track->Cluster(2)->TX()-Track->Cluster(0)->TX())/7.54;
            Sy = (Track->Cluster(2)->TY()-Track->Cluster(0)->TY())/7.54;


            //@Z=0 plane
            BsX = Track->fPlaner[2][0]; //x
            BsY = Track->fPlaner[2][1]; //y            
            
            //@X=0 plane
            TxY = Track->fPlaner[0][1]; //y
            TxZ = Track->fPlaner[0][2]; //z

            //@Y=0 plane
            TyX = Track->fPlaner[1][0]; //x
            TyZ = Track->fPlaner[1][2]; //z

            
            Float_t signx = 0; Float_t signy = 0;
            if (Track->LResidualX(0) < 0 ){signx = -1;}
            else {signx=1;}
            if (Track->LResidualY(0) < 0 ){signy = -1;}
            else {signy=1;}

            Rx = signx*(pow(Track->LResidualX(0),2)+pow(Track->LResidualX(1),2)+pow(Track->LResidualX(2),2));
            Ry = signy*(pow(Track->LResidualY(0),2)+pow(Track->LResidualY(1),2)+pow(Track->LResidualY(2),2));


            nTracks->Fill(nBX);

            Slope_SxSy[nBX].push_back(std::make_pair(Sx, Sy));            
            BeamSpot_BxBy[nBX].push_back(std::make_pair(BsX, BsY));
            Zpos_XzYz[nBX].push_back(std::make_pair(TxZ, TyZ));            
            Residuals_RxRy[nBX].push_back(std::make_pair(sqrt(fabs(Rx)), sqrt(fabs(Ry))));

          }
        }
      }
      //    } // ientry > begin
  } // ientry
  

  TH2F *allS = new TH2F("allS","Sx vs Sy,All",100,-0.1,0.1,100,-0.1,0.12);  
  TH2F *cbS = new TH2F("cbS","Sx vs Sy,Colliding Bunches",100,-0.1,0.1,100,-0.1,0.12);
  TH2F *ncbS = new TH2F("ncbS","Sx vs Sy,Non Colliding Bunches",20,-0.1,0.1,20,-0.1,0.12);
  TH2F *restS = new TH2F("restS","Sx vs Sy, Rest",20,-0.1,0.1,20,-0.1,0.12);    


  TH2F *allB = new TH2F("allB","Bx vs By,All",100,-15,15,100,-15,15);  
  TH2F *cbB = new TH2F("cbB","Bx vs By,Colliding Bunches",100,-15,15,100,-15,15);
  TH2F *ncbB = new TH2F("ncbB","Bx vs By,Non Colliding Bunches",20,-15,15,20,-15,15);
  TH2F *restB = new TH2F("restB","Bx vs By, Rest",20,-15,15,20,-15,15);    

  TH2F *allT = new TH2F("allT","TxZ vs TyZ,All",300,-325,325,300,-325,325);
  TH2F *cbT = new TH2F("cbT","TxZ vs TyZ,Colliding Bunches",300,-325,325,300,-325,325);
  TH2F *ncbT = new TH2F("ncbT","TxZ vs TyZ,Non Colliding Bunches",50,-325,325,50,-325,325);
  TH2F *restT = new TH2F("restT","TxZ vs TyZ, Rest",50,-325,325,50,-325,325);


  TH2F *allR = new TH2F("allR","Rx vs Ry,All",100,-0.01,0.1,100,-0.01,0.1);
  TH2F *cbR = new TH2F("cbR","Rx vs Ry,Colliding Bunches",100,-0.01,0.1,100,-0.01,0.1);
  TH2F *ncbR = new TH2F("ncbR","Rx vs Ry,Non Colliding Bunches",20,-0.01,0.1,20,-0.01,0.1);
  TH2F *restR = new TH2F("restR","Rx vs Ry, Rest",20,-0.01,0.1,20,-0.01,0.1);

  

  int max = 0;
  for (int i=0; i< 3564; i++) {
    //    std::cout << cb[i] << std::endl;    
    if (max < cb[i] ){
      max = cb[i];
    }
  }
  
  unsigned  int ncbThresh = max / 10;
  unsigned  int emptyThresh = ncbThresh / 2; 
  std::cout << "Max value is " << max << std::endl;  
  std::cout << "NcbThreshold: " << ncbThresh << " EmptyThreshold: " << emptyThresh  << std::endl;  

  
  unsigned int vecSize = 0;  // Num of entries in a bin/BX

  unsigned int cbNum = 0;    // Num of colliding bunches with tracks
  unsigned int ncbNum = 0;   // Num of non-colliding bunches with tracks
  unsigned int ebNum = 0;    // Num of empty bunches with tracks

  // Loop through all the bunches
  for(std::vector<std::vector<std::pair<Float_t,Float_t> > >::size_type i = 0; i != Slope_SxSy.size(); i++) {
    
    std::vector<std::pair<Float_t,Float_t> > vecS = Slope_SxSy[i];
    std::vector<std::pair<Float_t,Float_t> > vecB = BeamSpot_BxBy[i];
    std::vector<std::pair<Float_t,Float_t> > vecR = Residuals_RxRy[i];
    std::vector<std::pair<Float_t,Float_t> > vecT = Zpos_XzYz[i];    
    
    vecSize = Slope_SxSy[i].size(); // no of entries with this BX
    // all of it
    for ( std::vector<std::pair<Float_t,Float_t> >::size_type j = 0; j !=vecSize; j++) {
      allS->Fill(vecS[j].first,vecS[j].second);
      allB->Fill(vecB[j].first,vecB[j].second);
      allR->Fill(vecR[j].first,vecR[j].second);
      allT->Fill(vecT[j].first,vecT[j].second);      
    }

    // If BX has some entries
    if (vecSize != 0) {
      // If the number of entries is above non-colliding bunch threshold
      if ( vecSize > ncbThresh ) {
        cbNum++;        
        for ( std::vector<std::pair<Float_t,Float_t> >::size_type j = 0; j !=vecSize; j++) {
          cbS->Fill(vecS[j].first,vecS[j].second);
          cbB->Fill(vecB[j].first,vecB[j].second);
          cbR->Fill(vecS[j].first,vecR[j].second);
          cbT->Fill(vecT[j].first,vecT[j].second);          
          //          std::cout << vecS[j].first << " " << vecS[j].second << std::endl;            
        }
      } else if ( vecSize < ncbThresh && vecSize > emptyThresh ) {
        // Else if the number of entries is less than non-colliding bunch threshold and greater than empty threshold. Empty = not from colliding, not from non-colliding
        // i.e. non-colliding bunches only
        ncbNum++;        
        for ( std::vector<std::pair<Float_t,Float_t> >::size_type j = 0; j !=vecSize; j++) {
          ncbS->Fill(vecS[j].first,vecS[j].second);
          ncbB->Fill(vecB[j].first,vecB[j].second);
          ncbR->Fill(vecR[j].first,vecR[j].second);
          ncbT->Fill(vecT[j].first,vecT[j].second);          
        }
      } else {
        // If the number of entries is less than empty threshold
        // =Rest of the entries. 
        ebNum++;        
        for ( std::vector<std::pair<Float_t,Float_t> >::size_type j = 0; j !=vecSize; j++) {
          restS->Fill(vecS[j].first,vecS[j].second);
          restB->Fill(vecB[j].first,vecB[j].second);
          restR->Fill(vecR[j].first,vecR[j].second);
          restT->Fill(vecT[j].first,vecT[j].second);          
        }
      }
    }
  }

  std::cout << "There were \n" << cbNum << " colliding bunches\n" << ncbNum<<" non-colliding bunches\n" << ebNum <<" empty filled bunches\n" << "between " << begin <<"  and " << end <<  " slink events for this fill " << std::endl;
  
  
  f->Write();
  f->Close();

  
  return 0;
}

int main (int argc, char* argv[])
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " [DataFile.dat] [GainCal.dat] [AlignmentFile.dat]" << std::endl;
    return 1;
  }

  std::string const DataFileName = argv[1];
  std::string const GainCalFileName = argv[2];
  std::string const AlignmentFileName = argv[3];

  BXTracks(DataFileName, GainCalFileName, AlignmentFileName);

  return 0;
}


