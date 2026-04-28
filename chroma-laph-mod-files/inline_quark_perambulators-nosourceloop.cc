#include "inline_quark_perambulators.h"
#include "chroma.h"
using namespace std;

#ifdef TESTING
#include "tests.h"
#endif


namespace Chroma {
using namespace LaphEnv;

#if (QDP_ND == 4)

  using namespace LaphEnv;
  namespace InlineQuarkPerambulatorsEnv {

    //  The crucial create measurement routine. Must be in the *.cc
    //  so that it is local to this file.  Dynamically allocates
    //  and instantiates an object of our class "QuarkPerambulatorsInlineMeas".

AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
                                        const std::string& path) 
{
 return new QuarkPerambulatorsInlineMeas(xml_in, path);
}


const std::string name = "QUARK_PERAMBULATORS";

    // Registration boolean hidden in anonymous namespace.
namespace {
   bool registered = false;
}

    // Register all the factories.  This function may be called many
    // times by other measurements, so we only want to register this
    // inline measurement once.  Hence, the use of the "registered"
    // boolean above (which must be hidden in an anonymous namespace).

bool registerAll() 
{
 bool success = true; 
 if (!registered){
    success &= TheInlineMeasurementFactory::Instance().registerObject(
                      name, createMeasurement);
    registered = true;}
 return success;
}

// *********************************************************************
	

bool QuarkPerambulatorsInlineMeas::setComputationSet(int nLaphEigvecs,
                GaugeConfigurationInfo &gaugeinfo)
{
 XmlReader xml_rdr(xml_rd);
 if (xml_tag_count(xml_rdr,"ComputationSet")!=1) return false;
 
 XmlReader xmlrd(xml_rdr,"./descendant-or-self::ComputationSet");
 xmlread(xmlrd,"SourceTime",comp.SourceTime,"QUARK_PERAMBULATORS");
 if (xml_tag_count(xmlrd,"SourceLaphEigvecIndices")==1){
    multi1d<int> srcev_inds;
    xmlread(xmlrd,"SourceLaphEigvecIndices",srcev_inds,"QUARK_PERAMBULATORS");
    for (int k=0;k<srcev_inds.size();++k){
       comp.SourceLaphEigvecIndices.insert(srcev_inds[k]);}}
 else if (  (xml_tag_count(xmlrd,"SourceLaphEigvecIndexMin")==1)
          &&(xml_tag_count(xmlrd,"SourceLaphEigvecIndexMax")==1)){
    int sevmin=-1, sevmax=-1;
    xmlread(xmlrd,"SourceLaphEigvecIndexMin",sevmin,"QUARK_PERAMBULATORS");
    xmlread(xmlrd,"SourceLaphEigvecIndexMax",sevmax,"QUARK_PERAMBULATORS");
    if (sevmin<0) sevmin=0;
    if (sevmax>=nLaphEigvecs) sevmax=nLaphEigvecs-1;
    if (sevmax<sevmin) sevmax=sevmin;
    for (int ev=sevmin;ev<=sevmax;++ev){
       comp.SourceLaphEigvecIndices.insert(ev);}}
 else{
    for (int ev=0;ev<nLaphEigvecs;++ev){
       comp.SourceLaphEigvecIndices.insert(ev);}}
 multi1d<int> sink_times;
 xmlread(xmlrd,"SinkTimes",sink_times,"QUARK_PERAMBULATORS");
 for (int k=0;k<sink_times.size();++k){
    comp.SinkTimes.insert(sink_times[k]);}
 QDPIO::cout << endl << "QUARK_PERAMBULATORS computation set:"<<endl;
 QDPIO::cout <<endl<< "Source Time: "<<comp.SourceTime<<endl;
 QDPIO::cout <<"Source Laph Eigvec Indices: ";
 for (set<int>::const_iterator st=comp.SourceLaphEigvecIndices.begin();
                st!=comp.SourceLaphEigvecIndices.end();++st){
    QDPIO::cout <<" "<<*st;}
 QDPIO::cout <<endl<<"Sink Times: ";
 for (set<int>::const_iterator st=comp.SinkTimes.begin();
                st!=comp.SinkTimes.end();++st){
    QDPIO::cout <<" "<<*st;}
 QDPIO::cout <<endl<<endl;
    // reset gaugeinfo min and max time so the perambulator handler
    // does not look for laph eigenvectors at times that are not needed
 /*uint mintime=comp.SourceTime;
 uint maxtime=comp.SourceTime;
 for (set<int>::const_iterator it=comp.SinkTimes.begin();it!=comp.SinkTimes.end();++it){
    if (*it>maxtime) maxtime=*it;
    if (*it<mintime) mintime=*it;}
 gaugeinfo.resetMinTime(mintime);
 gaugeinfo.resetMaxTime(maxtime);*/
 return true;
}


bool QuarkPerambulatorsInlineMeas::setMergeInfo(GaugeConfigurationInfo& gaugeinfo)
{
 XmlReader xml_rdr(xml_rd);
 if (xml_tag_count(xml_rdr,"Merge")!=1) return false;
 XmlReader xmlrd(xml_rdr,"./descendant-or-self::Merge");
 mrgfiles=new FileListInfo(xmlrd);
 return true;
}


// *********************************************************************
	
     // Subroutine which does all of the work!!  Input parameters
     // must be as shown (specified by Chroma).  Actual input to
     // this routine is through the private data member
     //     XMLReader xlm_rdr


void QuarkPerambulatorsInlineMeas::operator()(unsigned long update_no,
                                              XMLWriter& xmlout) 
{

 XmlReader xml_rdr(xml_rd);

    // create the handler and set up the info from the
    // XML <QuarkPerambulatorInfo> tag

 if (xml_tag_count(xml_rdr,"QuarkPerambulatorInfo")!=1){
    QDPIO::cerr << "Must have one <QuarkPerambulatorInfo> tag"<<endl;
    QDP_abort(1);}
 XmlReader xmlr(xml_rdr,"./QuarkPerambulatorInfo");
 string gauge_xml;
 GaugeConfigurationInfo gaugeinfo(xmlr);
 //GaugeConfigurationInfo gaugeinfo(xmlr,gauge_xml);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 string smeared_quark_filestub;   
 xmlread(xmlr,"SmearedQuarkFileStub",smeared_quark_filestub,
         "QUARK_PERAMBULATORS");
 QuarkActionInfo quark(xmlr);
 //QuarkActionInfo quark(xmlr,gaugeinfo);
 FileListInfo files(xmlr);
 InverterInfo invinfo(xmlr);
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>=1){
    upper_spin_only=true;}
 bool verbose=false;
 string verbosity;
 xmlreadif(xmlr,"Verbosity",verbosity,"QUARK_PERAMBULATORS");
 if (tidyString(verbosity)=="full") verbose=true;

 QDPIO::cout << endl << endl;
 QDPIO::cout << " ***********************************************************"<<endl;
 QDPIO::cout << " *                                                         *"<<endl;
 QDPIO::cout << " *   Laph Task 2: Compute the quark perambulators          *"<<endl;
 QDPIO::cout << " *                and write to file                        *"<<endl;
 QDPIO::cout << " *                                                         *"<<endl;
 QDPIO::cout << " ***********************************************************"<<endl;
 QDPIO::cout << endl;
 //QDPIO::cout <<endl<<gaugeinfo.output()<<endl;
 QDPIO::cout << "XML header in the gauge configuration:"<<endl;
 QDPIO::cout << gauge_xml<<endl<<endl;
 QDPIO::cout <<endl<<endl<<"Gluon Smearing:"<<endl<<gSmear.output()<<endl<<endl;
 QDPIO::cout <<endl<<endl<<"Quark Smearing:"<<endl<<qSmear.output()<<endl<<endl;
 if (upper_spin_only){
    QDPIO::cout <<endl<< "Only upper spin components used"<<endl;}
 else{
    QDPIO::cout <<endl<< "All spin components used"<<endl;}
 QDPIO::cout <<"SmearedQuarkFileStub: "<<smeared_quark_filestub<<endl;
 QDPIO::cout << endl<<"QuarkAction:"<<endl<< quark.output()<<endl<<endl;
 QDPIO::cout << endl<<"Inverter Info:"<<endl<<invinfo.output()<<endl;

 XmlBufferWriter xml_out;
 push(xml_out,"QUARK_PERAMBULATORS");
 gaugeinfo.output(xml_out);
 gSmear.output(xml_out);
 qSmear.output(xml_out);
 quark.output(xml_out);
 invinfo.output(xml_out);
 pop(xml_out);
 xmlout << xml_out.str();

 //int mintime_orig=gaugeinfo.getMinTime();
 //int maxtime_orig=gaugeinfo.getMaxTime();
 bool mergemode=false;

    // read the source time, source laph eigvec indices, sink times
 bool cflag=setComputationSet(qSmear.getNumberOfLaplacianEigenvectors(),gaugeinfo);
    // if not computation, then get merge info
 if (xml_tag_count(xml_rdr,"PerambulatorInsertionInfo")==0){
    if (!cflag){
    setMergeInfo(gaugeinfo);
    mergemode=true;}
 
  // create handler
 PerambulatorHandler Q(gaugeinfo,gSmear,qSmear,quark,files,
                   smeared_quark_filestub,upper_spin_only,mergemode);

  // if not cflag, then doing a merge
 if (!cflag){
        Q.mergeData(*mrgfiles);
        delete mrgfiles;
        return;}

  // set the inverter info
 Q.setInverter(invinfo);
 QDPIO::cout << endl <<"Info initialized in QuarkHandler"<<endl;
 Q.outputSuffixMap();
  //gaugeinfo.resetMinTime(mintime_orig);
  //gaugeinfo.resetMaxTime(maxtime_orig);

 START_CODE();
 StopWatch outer;
 outer.start();
 QDPIO::cout <<endl<<endl;
 QDPIO::cout <<" *************************************************************"<<endl;
 QDPIO::cout <<" *"<<endl;
 QDPIO::cout <<" *  Now starting computation:"<<endl;
 QDPIO::cout <<" *"<<endl;
 QDPIO::cout <<" *************************************************************"<<endl;

 Q.computePerambulators(comp.SourceTime,comp.SourceLaphEigvecIndices,
                        comp.SinkTimes,verbose);

 outer.stop();
 QDPIO::cout << name << ": total time = " << outer.getTimeInSeconds() 
          << " secs" << endl;
 QDPIO::cout << name << ": ran successfully" << endl;
 END_CODE();
 } else if (xml_tag_count(xml_rdr,"PerambulatorInsertionInfo")==1){
    QDPIO::cerr << " <PerambulatorInsertionInfo> tag detected."<<endl;
    QDPIO::cerr << "Perambulator Insertion Computation is initiated."<<endl;

    XmlReader xmlri(xml_rdr,"./PerambulatorInsertionInfo");
    FileListInfo insertionfiles(xmlri);
    PerambulatorInsertionInfo insertinfo(xmlri);

    QDPIO::cout << endl<<"Insertion Info:"<<endl<<insertinfo.output()<<endl;
    
      // create insertion handler 
    PerambulatorInsertionHandler QI(gaugeinfo,gSmear,qSmear,quark,insertionfiles,
                       smeared_quark_filestub,insertinfo,upper_spin_only,mergemode);
    
      // if not cflag, then doing a merge
    if (!cflag){
        QI.mergeData(*mrgfiles);
        delete mrgfiles;              //<<< merge has to be written for insertion
        return;}

      // set the inverter info
    QI.setInverter(invinfo);
    QDPIO::cout << endl <<"Info initialized in QuarkHandler"<<endl;
    QI.outputSuffixMap();
    //gaugeinfo.resetMinTime(mintime_orig);
    //gaugeinfo.resetMaxTime(maxtime_orig);

    START_CODE();
    StopWatch outer;
    outer.start();
    QDPIO::cout <<endl<<endl;
    QDPIO::cout <<" *************************************************************"<<endl;
    QDPIO::cout <<" *"<<endl;
    QDPIO::cout <<" *  Now starting computation:"<<endl;
    QDPIO::cout <<" *"<<endl;
    QDPIO::cout <<" *************************************************************"<<endl;

    QI.computePerambulatorsInsertion(comp.SourceTime,comp.SourceLaphEigvecIndices,
                        comp.SinkTimes,verbose);

    outer.stop();
    QDPIO::cout << name << ": total time = " << outer.getTimeInSeconds() 
             << " secs" << endl;
    QDPIO::cout << name << ": ran successfully" << endl;
    END_CODE();   
 } else {
    QDPIO::cerr << "There are multiple <PerambulatorInsertionInfo> tags detected."<<endl;
    QDPIO::cerr << "Perambulator Computation is finished."<<endl;

 }

} 

// ******************************************************************
  }
#endif
} // namespace Chroma
