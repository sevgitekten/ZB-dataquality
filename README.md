# ZB-dataquality
#include <iostream>
#define ZBlowpu_cxx
#include "ZBlowpu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include <set>
#include "TRandom3.h"
#include <vector>
#include "basic.h"
#include <string>
#include <map>
#include <utility>
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;

void ZBlowpu::Loop()
{
    if (fChain == 0) return;
    TChain *chain=new TChain("ak4/ProcessedTree","analizpu"); // ak4'un icindeki Processed Tree deki treeleri chaine ekliyoruz
    chain->Add("ZB_*.root");      // -> 1GB
    ZBlowpu *t=new ZBlowpu(chain);
    
    //histogramımızı chainden cekelim
    TH1F *triggers = dynamic_cast<TH1F*>(fChain->GetCurrentFile()->Get("ak4/TriggerNames"));
    TAxis *xax = triggers->GetXaxis(); //histogramdan x axisindeki degerleri cektik
    vector <string> vectrigname; //trigger isimlerini koyacagımız bir vector olusturduk
    for (int trgidx = xax->GetFirst(); trgidx <= xax->GetLast(); trgidx++)
    {   
        string trgName = xax->GetBinLabel(trgidx); //triggerları bin bin alıyor
        if (trgName.compare("")==0) continue; // bos binleri atla
        vectrigname.push_back(trgName);
        cout<<"Trigger "<<trgidx<<" :"<<trgName<<endl;
        
    }
    
    Long64_t nentries = t->fChain->GetEntries();
    //nentries = 200000;
    cout<<"nentries"<<nentries<<endl;
    
    Long64_t nbytes = 0, nb = 0;
    // TFile myFile("lowpuHEGtrigturnon.root", "RECREATE");
    TFile myFile("5Nisan_alleta_daq_ZB.root", "RECREATE");
    
    //CondFormat'ın icinden cektigimiz text dosyaları,degistirdiklerim var!!///
    JetCorrectorParameters *pfchs_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L1FastJet_AK4PFchs.txt"); //burayı degistirdim//
    JetCorrectorParameters *pfchs_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L2Relative_AK4PFchs.txt");//burayı degistirdim//
    JetCorrectorParameters *pfchs_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L3Absolute_AK4PFchs.txt"); //burayı degistirdim//
    JetCorrectorParameters *pfchs_l2res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L2Residual_AK4PFchs.txt"); //burayı degistirdim//
    
    vector<JetCorrectorParameters> vParam_pfchs;
    vParam_pfchs.push_back(*pfchs_l1);
    vParam_pfchs.push_back(*pfchs_l2);
    vParam_pfchs.push_back(*pfchs_l3);
    vParam_pfchs.push_back(*pfchs_l2res);
    FactorizedJetCorrector *pfchs_jec = new FactorizedJetCorrector(vParam_pfchs);
    
    bool isdijet_p;
    Float_t pf_jtpt[100],pf_jt_eta[100],pf_jt_phi[100],ptdijet_pf;
    static const int netabins = 10;
    //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.7,4.7};
    static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,3.7,4.2,4.7};
    //static const double etabins[netabins+1] = {1.5,3.0};
    const double x[10][35]=
    {
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507},
        {1,5,6,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507}
    
    };
    
    //const int nx[1] = {34};
    const int nx[10] = {34,34,34,34,34,34,34,34,34,34};
	//const int nx[8] = {34,34,34,34,34,34,34,34};
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    std::vector<TProfile*> vectptt;
    std::vector<TProfile*> vectchf;
    std::vector<TProfile*> vectnhf;
    std::vector<TProfile*> vectnemf;
    std::vector<TProfile*> vectcemf;
    std::vector<TProfile*> vectmuf;
    //std::vector<TH1F*> vectSumetvsEvents;
    
    std::vector< std::vector< std::vector< TProfile*> > > VVdaq;
    std::vector< std::vector< TProfile*> > V2daq;
    
    //butun eventlerın histogramlarını buraya tanımlıyoruz//
    TH1F *hpt[netabins];
    
    for ( int k = 0; k<netabins; ++k)
    {
	
    TProfile *chf= new TProfile(Form("chf_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
	vectptt.push_back (chf);chf->Sumw2();
    TProfile *nhf= new TProfile(Form("nhf_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
	vectptt.push_back (nhf);nhf->Sumw2();
    TProfile *nemf= new TProfile(Form("nemf_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
	vectptt.push_back (nemf);nemf->Sumw2();
    TProfile *cemf= new TProfile(Form("cemf_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
	vectptt.push_back (cemf);cemf->Sumw2();
    TProfile *muf= new TProfile(Form("muf_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
	vectptt.push_back (muf);muf->Sumw2();
	
	V2daq.push_back(vectptt);vectptt.clear(); 
    }
    VVdaq.push_back(V2daq);V2daq.clear(); 
    
    TLorentzVector p4pf;
    std::vector<double>pfchs_ptcorr;
        int pf_njt = t->PFJetsCHS__;

    //---------------------------events dongusu icine girdik---------
  for (int i=0; i<nentries; i++)  
    {
        
        cout<<"\r"<<"event number: "<<i<<"/ "<<nentries<<flush;
        chain->GetEntry(i);
           
        if (t->TriggerDecision_.empty()) continue;
        int pf_njt =t->PFJetsCHS__;	//corr
        pfchs_ptcorr.resize(pf_njt);	//corr
        float pf_jtpt[100],pf_jt_eta[100];	//corr
        int itrig=0;
		
        for(int j=0; j<t->PFJetsCHS__; j++)
        {
            
            p4pf.SetPxPyPzE(t->PFJetsCHS__P4__fCoordinates_fX[j],t->PFJetsCHS__P4__fCoordinates_fY[j],
                            t->PFJetsCHS__P4__fCoordinates_fZ[j],t->PFJetsCHS__P4__fCoordinates_fT[j]);
             
            ////////////correction icin eklediklerimmm///
            pf_jtpt[j]   = p4pf.Pt();
            pf_jt_eta[j] = p4pf.Eta();
            
            float rho=t->EvtHdr__mPFRho;
            double area=t->PFJetsCHS__area_[j];
            //double ptraw = pf_jtpt[j]/PFJetsCHS__cor_[j]+rho*area;
            double ptraw = pf_jtpt[j]+rho*area;
            //cout<< "   rho: "<<rho<<"   area:  "<<area<<"   raw:  "<<ptraw<<endl;
            pfchs_jec->setJetEta(pf_jt_eta[j]);
            pfchs_jec->setJetPt(ptraw);
            pfchs_jec->setRho(rho);
            pfchs_jec->setJetA(area);
            
            pfchs_ptcorr[j]= ptraw*pfchs_jec->getCorrection();
            //cout<<"pfchs_ptcorr :"<< pfchs_ptcorr[j]<<endl;
            
            double pf_pt=pfchs_ptcorr[j];
            // double pf_pt=pf_jtpt[j];
            double pf_jteta=pf_jt_eta[j];
            
            ///////correction bolgesi bitti////////
            
            if ( t->PFJetsCHS__tightID_[j] && (t->PFMet__et_<0.3*t->PFMet__sumEt_) && (t->PFJetsCHS__nhf_[j]<0.9) && (t->PFJetsCHS__nemf_[j]<0.9)&& (t->PFJetsCHS__muf_[j]<0.9))
				{
                
                for ( int k=0;k<netabins;++k){
                    
                    if ( fabs(pf_jteta)>=etabins[k] && fabs(pf_jteta)<etabins[k+1]){
                                            
                        for(int trnameindex=0; trnameindex<1; trnameindex++)
                        {
							for(int trdecindex=0; trdecindex<t->TriggerDecision_.size(); trdecindex++)
                            {
								if(t->TriggerDecision_[trdecindex]==trnameindex)
                                {
                        VVdaq[trnameindex][k][0]->Fill(pf_pt,t->PFJetsCHS__chf_[j]);
                        VVdaq[trnameindex][k][1]->Fill(pf_pt,t->PFJetsCHS__nhf_[j]);
                        VVdaq[trnameindex][k][2]->Fill(pf_pt,t->PFJetsCHS__nemf_[j]);
                        VVdaq[trnameindex][k][3]->Fill(pf_pt,t->PFJetsCHS__cemf_[j]);
                        VVdaq[trnameindex][k][4]->Fill(pf_pt,t->PFJetsCHS__muf_[j]);
                        
                    }
                } //trnameindex
               } // 
            } /// fabs ifi bitiyor
          } ///etabin bitiyor
        } ///tightid li olan if cutı burda bitiyor
        } //jetler bitiyor  
    } ///event bitiyor
    
    myFile.cd();
    myFile.Write();
    myFile.Close();
}


