#ifndef __INLINE_QUARK_PERAMBULATORS_H__
#define __INLINE_QUARK_PERAMBULATORS_H__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "perambulator_handler.h"
#include "perambulator_insertion_handler.h"
#include "xml_help.h"


// ******************************************************************
// *                                                                *
// * Driver inline measurement that constructs the quark line ends  *
// * (sinks only).  Assumes that the eigenvectors of the            *
// * gauge-smeared covariant Laplacian are available.  Must be run  *
// * in 4D chroma_laph.  XML input must have the form:              *
// *                                                                *
// *  <chroma>                                                      *
// *   <RNG><Seed> ... </Seed></RNG>                                *
// *   <Cfg> ... </Cfg>                                             *
// *   <Param>                                                      *
// *    <nrow>12 12 12 96</nrow>   # lattice size Nx Ny Nz Nt       *
// *    <InlineMeasurements>                                        *
// *    <elem>                                                      *
// *                                                                *
// *     <Name> QUARK_PERAMBULATORS </Name>                         *
// *     <QuarkPerambulatorInfo>                                    *
// *       <GaugeConfigurationInfo> ... </GaugeConfigurationInfo>   *
// *       <GluonStoutSmearingInfo> ... </GluonStoutSmearingInfo>   *
// *       <QuarkLaphSmearingInfo> ... </QuarkLaphSmearingInfo>     *
// *       <SmearedQuarkFileStub> ... </SmearedQuarkFileStub>       *
// *       <QuarkActionInfo> ... </QuarkActionInfo>                 *
// *       <FileListInfo> ... </FileListInfo>                       *
// *       <UpperSpinComponentsOnly/>  (optional)                   *
// *       <InvertParam> ... </InvertParam>                         *
// *       <Verbosity> ... </Verbosity>  (optional)                 *
// *     </QuarkPerambulatorInfo>                                   *
// *                                                                *
// *     <ComputationSet>                                           *
// *       <SourceTime>23</SourceTime>                              *
// *       <SourceLaphEigvecIndices>3 5 9</SourceLaphEigvecIndices> *
// *             or                                                 *
// *       <SourceLaphEigvecIndexMin>3</SourceLaphEigvecIndexMin>   *
// *       <SourceLaphEigvecIndexMax>32</SourceLaphEigvecIndexMax>  *
// *       <SinkTimes>24 25</SinkTimes>                             *
// *     </ComputationSet>                                          *
// *                     OR                                         *
// *     <Merge>                                                    *
// *       <FileListInfo> ... </FileListInfo>                       *
// *     </Merge>                                                   *
// *                                                                *
// *    </elem>                                                     *
// *    </InlineMeasurements>                                       *
// *   </Param>                                                     *
// *  </chroma>                                                     *
// *                                                                *
// * If the tag "Verbosity" is included with value "full", then the *
// * quark sink solutions will be echoed to standard output.        *
// *                                                                *
// * Usually all of the perambulators for all source ev indices     *
// * cannot be done in a single run, and so for safety, the         *
// * data for different sets of source ev indices will get saved    *
// * into separate files.  These then need to be combined into      *
// * a single file for subsequent reading to compute the pion       *
// * and nucleon correlator.  The <Merge> tag above is meant        *
// * to accomplish this.                                            *
// *                                                                *
// ******************************************************************


namespace Chroma {
                 
#if (QDP_ND == 4) 

  namespace InlineQuarkPerambulatorsEnv {

 // **************************************************************


extern const std::string name;
bool registerAll();


class QuarkPerambulatorsInlineMeas : public AbsInlineMeasurement 
{

   XMLReader xml_rd;   // holds the XML input for this inline
                       // measurement, for use by the operator()
                       // member below

   struct ComputationSet {
      std::set<int> SourceTimes;
      std::set<int> SourceLaphEigvecIndices;
      std::set<int> SinkTimes;
   };
   
   ComputationSet comp;
   LaphEnv::FileListInfo* mrgfiles;

 public:

   QuarkPerambulatorsInlineMeas(XMLReader& xml_in, const std::string& path) 
                              : xml_rd(xml_in, path) {}

   ~QuarkPerambulatorsInlineMeas() {}
      
   bool setComputationSet(int nLaphEigvecs, LaphEnv::GaugeConfigurationInfo& gaugeinfo);

   bool setMergeInfo(LaphEnv::GaugeConfigurationInfo& gaugeinfo);

      //! Do the measurement
   void operator()(const unsigned long update_no, XMLWriter& xmlout); 

   unsigned long getFrequency() const {return 1;}
   
};
	

// ***********************************************************
  }
#endif
}

#endif
