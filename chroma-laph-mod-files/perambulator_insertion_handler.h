#ifndef PERAMBULATOR_INSERTION_HANDLER_H
#define PERAMBULATOR_INSERTION_HANDLER_H

#include "qdp.h"
#include "chromabase.h"
#include <memory>
#include "gauge_configuration_info.h"
#include "gauge_configuration_handler.h"
#include "xml_help.h"
#include "field_smearing_info.h"
#include "field_smearing_handler.h"
#include "quark_action_info.h"
#include "inverter_info.h"
#include "filelist_info.h"
#include "data_io_handler.h"
#include "perambulator_insertion.h"

namespace Chroma {
  namespace LaphEnv {


// *****************************************************************
// *                                                               *
// *  "PerambulatorHandler" handles computation of and subsequent  *
// *  access to the quark perambulators.                           *
// *                                                               *
// *  One of these handlers deals with quark perambulators         *
// *  for **one** set of info parameters, given by                 *
// *                                                               *
// *         GaugeConfigurationInfo                                *
// *         GluonSmearingInfo                                     *
// *         QuarkSmearingInfo                                     *
// *         Nspin  (4 or 2)                                       *
// *         QuarkActionInfo                                       *
// *         FileListInfo                                          *
// *                                                               *
// *  but one handler deals with different                         *
// *                                                               *
// *         src_time, snk_time, spin and Laph eigvec indices      *
// *                                                               *
// *  File structure and contents:                                 *
// *                                                               *
// *   - Results manipulated by one handler are contained in       *
// *     several files.  Each file has the same stub, but          *
// *     different positive integer suffices.                      *
// *             stub.0                                            *
// *             stub.1                                            *
// *             ...                                               *
// *             stub.N                                            *
// *     The files included are specified in a FileListInfo.       *
// *                                                               *
// *   - The header info in each file has a common part and a      *
// *     part that is specific to that one file:                   *
// *        <PerambulatorHandlerDataFile>                          *
// *           - common part                                       *
// *           - specific part                                     *
// *        </PerambulatorHandlerDataFile>                         *
// *                                                               *
// *   - The common header info includes                           *
// *         GaugeConfigurationInfo                                *
// *         GluonSmearingInfo                                     *
// *         QuarkSmearingInfo                                     *
// *         QuarkActionInfo                                       *
// *         Nspin  (4 or 2)                                       *
// *                                                               *
// *   - The specific header info includes (FileKey)               *
// *         int snk_time, src_time                                *
// *                                                               *
// *   - Each file contains several records whose access is given  *
// *     by a RecordKey.  The RecordKey contains an integer        *
// *     that specifies source spin and source eigenvec index.     *
// *     The data in each record (DataType) is a multi1d<Complex>  *
// *     containing "nev" * Nspin complex numbers, where "nev" is  *
// *     the number of Laph eigenvectors.   The Dirac-Pauli spin   *
// *     convention is used.                                       *
// *                                                               *
// *   - In 3d, the "get" function takes source and sink times     *
// *     and spins and returns a multi2d<Complex> which is a       *
// *     square matrix in terms of the LapH eigenvector indices.   *
// *     A smaller "neigsize" by "neigsize" matrix can also be     *
// *     returned.                                                 *
// *                                                               *
// *  All Laph Handlers follow the member naming convention:       *
// *                                                               *
// *    compute....()  to do original computation                  *
// *                                                               *
// *    get...()       provides access to results                  *
// *                                                               *
// *  When using "set" and "get", individual spin components and   *
// *  time slices are returned.  At this time, no quark            *
// *  covariant displacements are accommodated.                    *
// *                                                               *
// *  The usual use of a PerambulatorHandler is as follows:        *
// *                                                               *
// *   - to compute the quark sources/sinks:                       *
// *                                                               *
// *       PerambulatorHandler Q;  // declare                      *
// *       Q.setInfo(...);            // input common info         *
// *       Q.setInverter(...);        // input inverter info       *
// *       Q.computePerambulators(...);                            *
// *                                                               *
// *   - to subsequently use results to compute                    *
// *     hadron source/sinks:                                      *
// *                                                               *
// *       PerambulatorHandler Q;                                  *
// *       Q.setInfo(...);                                         *
// *       Q.getData(...);                                         *
// *                                                               *
// *  Usually all of the perambulators for all source ev indices   *
// *  cannot be done in a single run, and so for safety, the       *
// *  data for different sets of source ev indices will get saved  *
// *  into separate files.  These then need to be combined into    *
// *  a single file for subsequent reading to compute the pion     *
// *  and nucleon correlator.  The method                          *
// *                                                               *
// *     mergeData(...)                                            *
// *                                                               *
// *  is provided to accomplish this.                              *
// *                                                               *
// *****************************************************************


class PerambulatorInsertionHandler
{

 public:

   struct FileKey
   {
      int snk_time;
      int src_time;

      FileKey() : snk_time(0), src_time(0) {}
      FileKey(int in_snktime, int in_srctime);
      FileKey(XmlReader& xmlr);
      FileKey(const FileKey& rhs);
      FileKey& operator=(const FileKey& rhs);
      ~FileKey(){}
      void output(XmlWriter& xmlw) const;
      bool operator<(const FileKey& rhs) const;
      bool operator==(const FileKey& rhs) const;
      bool operator!=(const FileKey& rhs) const;
   };

  typedef multi1d<Complex> DataType;

 private:

       // pointers to internal infos (managed by this handler
       // with new and delete)

   const GaugeConfigurationInfo *uPtr;
   const GluonSmearingInfo *gSmearPtr;
   const QuarkSmearingInfo *qSmearPtr;
   const QuarkActionInfo *qactionPtr;
   const FileListInfo *fPtr;
   const InverterInfo *invertPtr;
   uint Nspin;
   bool merge_mode;

   const PerambulatorInsertionInfo *insertPtr;

       // sub-handler pointers

   static std::unique_ptr<QuarkSmearingHandler> qSmearHandler;
   static std::unique_ptr<GaugeConfigurationHandler> gaugeHandler;

   static int qSmearCounter;
   static int gaugeCounter;

       // Prevent copying ... handler might contain large
       // amounts of data

   PerambulatorInsertionHandler(const PerambulatorInsertionHandler&);
   PerambulatorInsertionHandler& operator=(const PerambulatorInsertionHandler&);


       // data I/O handler pointers

   DataPutHandlerMF<PerambulatorInsertionHandler,FileKey,UIntKey,DataType> *DHputPtr;
   DataGetHandlerMF<PerambulatorInsertionHandler,FileKey,UIntKey,DataType> *DHgetPtr;


 public:


   PerambulatorInsertionHandler();

      // in 4-d, "gauge_str" should be a gauge_id,
      // in 3-d, "gauge_str" should be the smeared gauge file name

   PerambulatorInsertionHandler(const GaugeConfigurationInfo& gaugeinfo,
                       const GluonSmearingInfo& gluonsmear,
                       const QuarkSmearingInfo& quarksmear,
                       const QuarkActionInfo& quark,
                       const FileListInfo& flist,
                       const std::string& smeared_quark_filestub,
                       const PerambulatorInsertionInfo& insertinfo,
                       bool upper_spin_components_only=false,
                       bool mergemode=false,
                       const std::string& gauge_str="default_gauge_field");
   
   void setInfo(const GaugeConfigurationInfo& gaugeinfo,
                const GluonSmearingInfo& gluonsmear,
                const QuarkSmearingInfo& quarksmear,
                const QuarkActionInfo& quark,
                const FileListInfo& flist,
                const std::string& smeared_quark_filestub,
                const PerambulatorInsertionInfo& insertinfo,
                bool upper_spin_components_only=false,
                bool mergemode=false,
                const std::string& gauge_str="default_gauge_field");

   ~PerambulatorInsertionHandler();

   void clear();

   bool isInfoSet() const;

   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

   const GluonSmearingInfo& getGluonSmearingInfo() const;

   const QuarkSmearingInfo& getQuarkSmearingInfo() const;

   const QuarkActionInfo& getQuarkActionInfo() const;

   const FileListInfo& getFileListInfo() const;

   const PerambulatorInsertionInfo& getPerambulatorInsertionInfo() const;

   uint getNumberOfLaplacianEigenvectors() const;

   int getTimeExtent() const;

   void getHeader(XmlWriter& xmlout) const;


#if (QDP_ND == 4)

   void getFileMap(XmlWriter& xmlout) const;

   void outputSuffixMap();

   void outputSuffixMap(TextFileWriter& fout);

   void setInverter(const InverterInfo& invinfo);

   const InverterInfo& getInverterInfo() const;

        // compute quark perambulators (exact distillation); useful for smearing studies

   void computePerambulatorsInsertion(int src_time, const std::set<int>& src_lapheigvec_indices,
                             const std::set<int>& snk_times, bool verbose=false);

   void mergeData(const FileListInfo& input_files);



#elif (QDP_ND == 3)

   void getSourceTimes(std::set<int>& source_times) const;

   multi2d<Complex> getData(int snk_time, int snk_spin, int src_time, int src_spin, 
                            int nEigsUse);

   bool queryData(int src_time, int src_spin, int snk_time, int snk_spin, int nEigsUse);

   QuarkSmearingHandler& getQuarkSmearingHandler()
    {return *qSmearHandler;}

   void clearGaugeData();
   
#endif

   

 private:


   void set_info(const GaugeConfigurationInfo& gaugeinfo,
                 const GluonSmearingInfo& gluonsmear,
                 const QuarkSmearingInfo& quarksmear,
                 const QuarkActionInfo& quark,
                 const FileListInfo& flist,
                 const std::string& smeared_quark_filestub,
		           const PerambulatorInsertionInfo& insertinfo,
                 bool upper_spin_components_only,
                 const std::string& gauge_str, bool mergemode);


   bool checkHeader(XmlReader& xmlr, int suffix);
   void writeHeader(XmlWriter& xmlout, const FileKey& fkey,
                    int suffix);

   void check_info_set(const std::string& name) const;
   void check_info_set_Insert(const std::string& name) const;

   	 // insertion multiplications

   void D1Insert(LatticeFermion& chi, const LatticeFermion& psi);
   void D1InsertTwist(LatticeFermion& chi, const LatticeFermion& psi);
   void D1InsertDeltaM(LatticeFermion& chi, const LatticeFermion& psi);


         //  sub-handler connections

   void connectGaugeConfigurationHandler(const std::string& gauge_id="default_gauge_field");
   void connectQuarkSmearingHandler(const std::string& smeared_quark_filestub);

   void disconnectGaugeConfigurationHandler();
   void disconnectQuarkSmearingHandler();


   friend class DataPutHandlerMF<PerambulatorInsertionHandler,FileKey,UIntKey,DataType>;
   friend class DataGetHandlerMF<PerambulatorInsertionHandler,FileKey,UIntKey,DataType>;

};



// **************************************************************************************
  }
}
#endif  
