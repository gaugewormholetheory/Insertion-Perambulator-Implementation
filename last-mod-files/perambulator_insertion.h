#ifndef PERAMBULATOR_INSERTION_H
#define PERAMBULATOR_INSERTION_H


#include "xml_handler.h"
//#include "gauge_configuration_info.h"

//namespace Chroma {
  namespace LaphEnv {

// ********************************************************************
// *                                                                  *
// *  Class "PerambulatorInsertionInfo" holds information about a     *
// *  insertion of a perambulator.  It also checks that its config    *
// *  info is the same when compared against another object of class  *
// *  PerambulatorInsertionInfo.                                      *
// *                                                                  *
// *  The constructor takes either an XmlReader or a string.  The     *
// *  XmlReader version is called as follows:                         *
// *                                                                  *
// *    XmlReader xml_in(...);                                        *
// *    GaugeConfigurationInfo U(xml_in);                             *
// *                                                                  *
// *  There are three ways of constructing a GaugeConfigurationInfo   *
// *  depending on whether a "<gauge_id>" or a "<getFromFile>" is     *
// *  given as a tag in the XmlReader:                                *
// *                                                                  *
// *    (1) If "<gauge_id>" is given as a tag:                        *
// *                                                                  *
// *      --> This constructor expects XML of the form                *
// *               <GaugeConfigurationInfo>                           *
// *                  <gauge_id>....</gauge_id>                       *
// *                  <ConfigType> .... <ConfigType>                  *
// *               </GaugeConfigurationInfo>                          *
// *          then information about the configuration is obtained    *
// *          from TheNamedObjMap.                                    *
// *                                                                  *
// *    (2) If "<getFromFile>" is given as a tag:                     *
// *                                                                  *
// *      --> This constructor expects XML of the form                *
// *               <GaugeConfigurationInfo>                           *
// *                  <getFromFile> gaugefilename </getFromFile>      *
// *               <ConfigType> .... <ConfigType>                     *
// *               </GaugeConfigurationInfo>                          *
// *          then information about the configuration is obtained    *
// *          from xml header in the gauge configuration file named   *
// *          in the tag.                                             *
// *                                                                  *
// *    (3) If neither "<gauge_id>" nor "<getFromFile>" is given as   *
// *        a tag:                                                    *
// *                                                                  *
// *      --> This version of the constructor expects all of the      *
// *          information to be given in the XML itself.              *
// *            <GaugeConfigurationInfo>                              *
// *               <HMCTrajectoryNumber>...<HMCTrajectoryNumber>      *
// *               <TimeDir> .... <TimeDir>                           *
// *               <ConfigType> .... <ConfigType>                     *
// *               <FileName> .... <FileName>                         *
// *               <TimeExtent> .... <TimeExtent>                     *
// *               <NumberOfDir> .... <NumberOfDir>                   *
// *               <LatticeExtents> .... <LatticeExtents>             *
// *            </GaugeConfigurationInfo>                             *
// *                                                                  *
// *                                                                  *
// *    string headerInfo(...);                                       *
// *    GaugeConfigurationInfo U(headerInfo);                         *
// *                                                                  *
// *                                                                  *
// *  Example usage:                                                  *
// *                                                                  *
// *    GaugeConfigurationInfo u2(xml_in);                            *
// *    u.checkEqual(u2); --> checks that u2 and u have same          *
// *                           content; throws exception if not       *
// *                                                                  *
// *    string out = u.output();   // xml output                      *
// *    string out = u.output(2);  // indented xml output             *
// *    int j = u.getTrajNum();    // returns  number of RHMC         *
// *                               // trajectory for this config      *
// *    int nt = u.getTimeExtent(); // time extent of lattice         *
// *    int tdir = u.getTimeDir();  // index of time                  *
// *                                                                  * 
// *  With open temporal boundary conditions, we only want to perform *
// *  calculations on the central half of the timeslices. The         *
// *  members "minTime" and "maxTime" restrict the calculations and   *
// *  storage output to this time range.                              *
// *                                                                  *
// ********************************************************************


class PerambulatorInsertionInfo
{

  private:

    std::string insertion_id;
    std::vector<int> twist_vector;

  public:
  
    PerambulatorInsertionInfo() = default;
 
    PerambulatorInsertionInfo(const XMLHandler& xml_in);

    PerambulatorInsertionInfo(const PerambulatorInsertionInfo& rhs);

    PerambulatorInsertionInfo& operator=(const PerambulatorInsertionInfo& rhs);

    ~PerambulatorInsertionInfo(){}

    void checkEqual(const PerambulatorInsertionInfo& rhs) const;


    //Added to do make current_handler.cc functional

    std::string output(int indent = 0) const;

    void output(XMLHandler& xmlout) const;

    std::vector<int> getTwistVector() const {return twist_vector;}
    
    std::string getInsertionName() const {return insertion_id;}


  private:

    void setTwistVector(XMLHandler& xmlr, std::vector<int>& twist_vector);
    void set_info(XMLHandler& xmlr);
};


// ****************************************************************
  }
//}
#endif

