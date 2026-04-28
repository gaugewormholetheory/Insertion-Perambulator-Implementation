//#include "inline_variational.h"
//#include "inline_correlators.h"
#include "inline_convert_perambulators.h"
#include "inline_convert_perambulators_insertion.h"
//include "inline_convert_mode_doublets.h"
//include "inline_convert_mode_triplets.h"
//include "inline_compute_two_point_corrs_dist.h"
//include "inline_compute_two_point_curr_corrs_dist.h"
//#include "inline_display.h"
//#include "inline_vevs.h"
#include "xml_handler.h"
#include "named_obj_map.h"
#include "stopwatch.h"
#include <cstdio>
#include <ctime>
#include <map>
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
using namespace LaphEnv;


void output_datetime()
{
 time_t rawtime;
 struct tm *timeinfo;
 time(&rawtime);
 timeinfo=localtime(&rawtime);
 cout << "  Current date/time: "<<asctime(timeinfo);
}


  //  Constructor sets up the known tasks.  Call
  //  member function "do_task" to perform the task.
  //  "TaskMap" is a map that associates a string
  //  containing a task name to a function pointer.

class Tasker
{
   typedef void (*task_ptr)(XMLHandler&);
   map<string, task_ptr> TaskMap;

 public:
    
   Tasker();
   ~Tasker(){}
   void do_task(XMLHandler& xml_in, bool echo);

};

     // set up the known tasks

Tasker::Tasker()
{
// TaskMap["LAPH_CORRELATORS"]=&doCorrelators;
 //TaskMap["LAPH_CORR_REFORMAT"]=&doCorrReformat;
 TaskMap["LAPH_PERAMB_CONVERT"]=&doPerambulatorConvert;
 TaskMap["LAPH_PERAMB_CONVERT_CL"]=&doPerambulatorConvertCL;
 TaskMap["LAPH_PERAMB_CONVERT_SPARSE_GRID"]=&doPerambulatorConvertSG;
 TaskMap["LAPH_PERAMB_INSERT_CONVERT_CL"]=&doPerambulatorInsertionConvertCL;
 //TaskMap["LAPH_MODE_DOUBLET_CONVERT"]=&doModeDoubletConvert;
 //TaskMap["LAPH_MODE_TRIPLET_CONVERT"]=&doModeTripletConvert;
 //TaskMap["LAPH_COMPUTE_TWO_POINT_CORRS_DIST"]=&doComputeTwoPointCorrsDist;
 //TaskMap["LAPH_COMPUTE_TWO_POINT_CURR_CORRS_DIST"]=&doComputeTwoPointCurrCorrsDist;
// TaskMap["LAPH_SINGLE_MESON_CORRELATOR"]=&doSingleMesonCorrelator;
 //TaskMap["LAPH_SINGLE_BARYON_CORRELATOR"]=&doSingleMesonCorrelator;
// TaskMap["LAPH_VEVS"]=&doVEVs;
 //TaskMap["LAPH_CORR_MERGE"]=&doCorrMerge;
// TaskMap["LAPH_VEV_MERGE"]=&doVEVMerge;
 //TaskMap["LAPH_VAR_IMPROVE_MESON"]=&doVariationalImproveMesons;
 //TaskMap["LAPH_VAR_IMPROVE_BARYON"]=&doVariationalImproveBaryons;
 //TaskMap["LAPH_DISPLAY_CORR"]=&doCorrDisplay;
 //TaskMap["LAPH_CHECK_CORR"]=&doCorrCheck;
// TaskMap["LAPH_DISPLAY_VEV"]=&doVEVDisplay;
// TaskMap["LAPH_CHECK_VEV"]=&doVEVCheck;
};


void Tasker::do_task(XMLHandler& xml_task, bool echo)
{
 if (echo){
    cout << "Input XML for this task:"<<endl
         <<xml_task.output()<<endl;}
 XMLHandler xmlt(xml_task);
 xml_child_assert(xmlt,"Name","do_task");
 if (!xmlt.is_simple_element()) 
    throw(std::runtime_error("Name tag is not simple XML element"));
 string task_name=xmlt.get_text_content();
 cout << "  Task name = "<<task_name<<endl;

 map<string,task_ptr >::iterator taskit=TaskMap.find(task_name);
 if (taskit!=TaskMap.end()){
    (*(taskit->second))(xml_task);}  // do the task!!
 else{
    throw(std::runtime_error("Unknown task name"));}   // unknown task?
}


// ****************************************************************
// *                                                              *
// *          Main driver program to run all tasks                *
// *                                                              *
// *   Program takes a single argument that is the name of the    *
// *   input file, and standard output is used.  Use redirection  *
// *   for output to a log file.  Input file must contain a       *
// *   single XML document with root tag named "LastLaph".        *
// *   Inside the root tag should be one or more <Task> tags,     *
// *   and inside each <Task> element should be a <Name> tag      *
// *   whose content is the name of the task.  The name must      *
// *   be one of the allowed names specified in the "do_task"     *
// *   subroutine.  If a tag <EchoXML/> is present as a child of  *
// *   the root tag, then the XML input is echoed to standard     *
// *   output.                                                    *
// *                                                              *
// *   Sample input XML:                                          *
// *                                                              *
// *    <LastLaph>                                                *
// *       <EchoXML/>                                             *
// *       <Task>                                                 *
// *          <Name>Task 1</Name>                                 *
// *       </Task>                                                *
// *       <Task>                                                 *
// *          <Name>Task 2</Name>                                 *
// *       </Task>                                                *
// *    </LastLaph>                                               *
// *                                                              *
// ****************************************************************

std::map<std::string, NamedObjBase*> NamedObjMap::the_map;


int main(int argc, const char* argv[])
{
 cout << endl << "Starting LAST_LAPH"<<endl;
 output_datetime();
  
 if (argc!=2){
    cerr << "  LastLaph requires one argument: "
         << " the name of an input XML file"<<endl;
    cerr << "   ... exiting..."<<endl;
    exit(1);}

 {StopWatch rolex,swatch;
 rolex.reset();
 rolex.start();
 string input_file = string( argv[1] );
 cout << "  Name of input XML file: "<<input_file<<endl<<endl;

 XMLHandler xml_in;
 xml_in.set_from_file(input_file);
 if (xml_in.fail()){
    cerr << "  Unable to read/parse XML content in input file"<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}
 if (xml_in.get_node_name()!="LastLaph"){
    cerr << "  Root tag of input XML must be named \"LastLaph\""<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

 bool echo=false;
 int ntasks = 0;
 try{
    xml_in.seek_first_child();
    while (xml_in.good()){
       if (xml_in.get_node_name()=="Task") ntasks++;
       else if (xml_in.get_node_name()=="EchoXML") echo=true;
       else throw(std::runtime_error("Invalid XML input"));
       xml_in.seek_next_sibling();}}
 catch(const std::exception& err){
    cerr << "  Error: "<<err.what()<<endl;
    cerr << "  ... exiting..."<<endl;
    exit(1);}

 cout << "  Number of tasks is "<<ntasks<<endl;
#if defined(_OPENMP)
 cout << "  Maximum number of threads is "<<omp_get_max_threads()<<endl;
#endif
 xml_in.seek_root();
 xml_in.seek_first_child();
 Tasker T;

 for (int task=1;task<=ntasks;++task){
    swatch.reset();
    swatch.start();
    try{
       cout << endl<<endl<<"Starting Task "<<task<<endl<<endl;
       if (xml_in.fail()) throw(std::runtime_error("XML input error"));
       if (xml_in.get_node_name()=="EchoXML")
          xml_in.seek_next_sibling();

       XMLHandler xml_task(xml_in);

           // the main task
       T.do_task(xml_task,echo);

       }
    catch(const std::exception& err){
       cerr << "Error on Task "<<task<<":  "<<err.what()<<endl;}
    catch(...){
       cerr << "Error on Task "<<task<<":  generic"<<endl;}
    swatch.stop();
    cout << "Task "<<task<<" done using time = " << swatch.getTimeInSeconds() 
         << " secs" << endl;
    xml_in.seek_next_sibling();
    }

 rolex.stop();
 cout <<endl<<endl;
 cout << "LAST_LAPH: total time = "<< rolex.getTimeInSeconds() 
      << " secs" << endl;}
 cout << "LAST_LAPH: completion" << endl;
 output_datetime();

 return 0;
}

