#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );


void hadd() {
   // Prepare the files to me merged
   if(gSystem->AccessPathName("hsimple1.root")) {
     gSystem->CopyFile("hsimple.root", "hsimple1.root");
     gSystem->CopyFile("hsimple.root", "hsimple2.root");
   }

   // in an interactive ROOT session, edit the file names
   // Target and FileList, then
   // root > .L hadd.C
   // root > hadd()

   Target = TFile::Open( "full_result2.root", "RECREATE" );

   FileList = new TList();
   FileList->Add( TFile::Open("histos_17903_18826.root") );
  // FileList->Add( TFile::Open("histos_18827_19803.root") ); //fission only
   FileList->Add( TFile::Open("histos_19804_21670.root") );
  //FileList->Add( TFile::Open("histos_21671_22430.root") ); //fission only
   FileList->Add( TFile::Open("histos_22431_26854.root") );
  // FileList->Add( TFile::Open("histos_26855_27822.root") ); //fission only
   FileList->Add( TFile::Open("histos_27823_38350.root") );
   FileList->Add( TFile::Open("histos_38352_49763.root") );
   FileList->Add( TFile::Open("histos_49773_51749.root") );
 // FileList->Add( TFile::Open("histos_51750_54730.root") ); //fission only
   FileList->Add( TFile::Open("histos_54731_73543.root") );  
   FileList->Add( TFile::Open("histos_73544_73678.root") );
   FileList->Add( TFile::Open("histos_73679_73857.root") );
   FileList->Add( TFile::Open("histos_73858_75801.root") );
   FileList->Add( TFile::Open("histos_75802_78633.root") );
   FileList->Add( TFile::Open("histos_78634_95328.root") );
   FileList->Add( TFile::Open("histos_95329_99518.root") );
   FileList->Add( TFile::Open("histos_99519_99999.root") );
   FileList->Add( TFile::Open("histos_100000_100273.root") );
//   FileList->Add( TFile::Open("histos_100274_101025.root") ); //fission only 
   FileList->Add( TFile::Open("histos_101026_102141.root") ); //needs shift for stilbene
  // FileList->Add( TFile::Open("histos_102142_102295.root") ); //fission only
//  FileList->Add( TFile::Open("histos_102296_102507.root") ); //fission only
   FileList->Add( TFile::Open("histos_102508_103391.root") ); //needs shift for stilbene
 //  FileList->Add( TFile::Open("histos_103392_103601.root") ); //fission only
  // FileList->Add( TFile::Open("histos_104794_104869.root") );
//  FileList->Add( TFile::Open("histos_104870_105155.root") ); //fission only
 //FileList->Add( TFile::Open("histos_105156_105635.root") ); //fission only
 //FileList->Add( TFile::Open("histos_105636_105728.root") ); //fission only
  //  FileList->Add( TFile::Open("histos_105729_106535.root") );
// FileList->Add( TFile::Open("histos_106536_107307.root") ); //fission only

   MergeRootfile( Target, FileList );

}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {

   //  cout << "Target path: " << target->GetPath() << endl;
   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );

   TFile *first_source = (TFile*)sourcelist->First();
   first_source->cd( path );
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   // loop over all keys in this directory
   TChain *globChain = 0;
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {

      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

      // read object from first source file
      first_source->cd( path );
      TObject *obj = key->ReadObj();

      if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
         // descendant of TH1 -> merge it

         //      cout << "Merging histogram " << obj->GetName() << endl;
         TH1 *h1 = (TH1*)obj;

         // loop over all source files and add the content of the
         // correspondant histogram to the one pointed to by "h1"
         TFile *nextsource = (TFile*)sourcelist->After( first_source );
         while ( nextsource ) {

            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
            if (key2) {
               TH1 *h2 = (TH1*)key2->ReadObj();
               h1->Add( h2 );
               delete h2;
            }

            nextsource = (TFile*)sourcelist->After( nextsource );
         }
      }
      else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

         // loop over all source files create a chain of Trees "globChain"
         const char* obj_name= obj->GetName();

         globChain = new TChain(obj_name);
         globChain->Add(first_source->GetName());
         TFile *nextsource = (TFile*)sourcelist->After( first_source );
         //      const char* file_name = nextsource->GetName();
         // cout << "file name  " << file_name << endl;
         while ( nextsource ) {

            globChain->Add(nextsource->GetName());
            nextsource = (TFile*)sourcelist->After( nextsource );
         }

      } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
         // it's a subdirectory

         cout << "Found subdirectory " << obj->GetName() << endl;

         // create a new subdir of same name and title in the target file
         target->cd();
         TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

         // newdir is now the starting point of another round of merging
         // newdir still knows its depth within the target file via
         // GetPath(), so we can still figure out where we are in the recursion
         MergeRootfile( newdir, sourcelist );

      } else {

         // object is of no type that we know or can handle
         cout << "Unknown object type, name: "
         << obj->GetName() << " title: " << obj->GetTitle() << endl;
      }

      // now write the merged histogram (which is "in" obj) to the target file
      // note that this will just store obj in the current directory level,
      // which is not persistent until the complete directory itself is stored
      // by "target->Write()" below
      if ( obj ) {
         target->cd();

         //!!if the object is a tree, it is stored in globChain...
         if(obj->IsA()->InheritsFrom( TTree::Class() ))
            globChain->Merge(target->GetFile(),0,"keep");
         else
            obj->Write( key->GetName() );
      }

   } // while ( ( TKey *key = (TKey*)nextkey() ) )

   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
}
