#include "quark_action_info.h"
//#include "gauge_configuration_handler.h"
#include "perambulator_insertion.h"
#include <sstream>

using namespace std;

//namespace Chroma {
  namespace LaphEnv {


// ************************************************************

PerambulatorInsertionInfo::PerambulatorInsertionInfo(const XMLHandler& xml_in)
{ 
 //xml_tag_assert(xml_in,"InsertionInfo");
 XMLHandler xmlr(xml_in, "InsertionInfo");
 set_info(xmlr);
}

void PerambulatorInsertionInfo::set_info(XMLHandler& xmlr)
{
 try {
    xmlread(xmlr, "Name", insertion_id,"InsertionInfo");}
 catch(...){
    //xml_cerr(xmlr,"Invalid read of <Name> in InsertionInfo");
    xmlreadfail(xmlr,"InsertionInfo","Unsupported name in InsertionInfo");}
 insertion_id=tidyString(insertion_id);
 if(insertion_id=="Twist"){
   setTwistVector(xmlr, twist_vector);
 }
}

void PerambulatorInsertionInfo::setTwistVector(XMLHandler& xmlr, vector<int>& twist_vector)
{
 string vec_str;
 vector<int> temp_vec;
 int val;

 if (xml_tag_count(xmlr,"twist_vector")>0){
    try{
       xmlread(xmlr, "twist_vector", vec_str,"InsertionInfo");
       istringstream iss(vec_str);
       int val;

       while (iss >> val) {
         temp_vec.push_back(val);}
      
       vector<int> my_vec(temp_vec.size());
       for (int i = 0; i < temp_vec.size(); ++i) {
         my_vec[i] = temp_vec[i];}
      
       twist_vector = my_vec;}

    catch(...){
       //xml_cerr(xmlr,"Invalid read of <twist_vector> in InsertionInfo");
       xmlreadfail(xmlr,"InsertionInfo","Invalid read of <twist_vector> in InsertionInfo");}}
 else{
    //xml_cerr(xmlr,"could not set twist vector in InsertionInfo");
    xmlreadfail(xmlr,"InsertionInfo","could not set twist vector in InsertionInfo");}
}

  // ************************************************************

     // copy constructor

PerambulatorInsertionInfo::PerambulatorInsertionInfo(const PerambulatorInsertionInfo& rhs)
                : insertion_id(rhs.insertion_id), twist_vector(rhs.twist_vector) {}

						
PerambulatorInsertionInfo& PerambulatorInsertionInfo::operator=(const PerambulatorInsertionInfo& rhs)
{
 insertion_id = rhs.insertion_id;
 if(rhs.insertion_id=="Twist"){
   twist_vector = rhs.twist_vector;
 }
 return *this;
}


void PerambulatorInsertionInfo::checkEqual(const PerambulatorInsertionInfo& rhs) const
{
 if ((insertion_id != rhs.insertion_id)||(twist_vector!=rhs.twist_vector)){
    std::cerr << "InsertionInfo checkEqual failed"<<endl;
    std::cerr << "LHS:"<<endl<<output()<<endl<<"RHS:"<<endl<<rhs.output()<<endl;
    throw(std::invalid_argument("InsertionInfo xml contents do not match"));}
}



string PerambulatorInsertionInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<InsertionInfo>"<<endl;
 oss << pad << "   <Name>" << insertion_id << "</Name>"<<endl;
 oss << pad << "   <twist_vector>" << twist_vector[0] << " " << twist_vector[1] << " " << twist_vector[2] << " " << twist_vector[3]  << "</twist_vector>"<<endl;
 oss << pad << "</InsertionInfo>"<<endl;
 return oss.str();
}

void PerambulatorInsertionInfo::output(XMLHandler& xmlout) const
{
 std::ostringstream oss;

 for (int i = 0; i < twist_vector.size(); ++i) {
    oss << twist_vector[i];
    if (i != twist_vector.size() - 1)
        oss << " ";  // add space between numbers
 }
   
 /*push(xmlout,"InsertionInfo");
 write(xmlout,"Name",insertion_id);
 write(xmlout,"twist_vector",oss.str());
 pop(xmlout);*/
 xmlout.set_root("InsertionInfo");
 xmlout.put_child("Name",insertion_id);
 xmlout.put_child("twist_vector",oss.str());
}



// ***********************************************************
  }
//}
