#include "perambulator_insertion_handler.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "util/ferm/diractodr.h"
#include "time_slices.h"
#include "meas/smear/displace.h"
#include "multi_compare.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "openqcd_handler.h"
#include <vector>
using namespace std;

#ifdef TESTING
#include "tests.h"
#endif

#if (QDP_ND == 4)
#if QDP_USE_LEXICO_LAYOUT == 1
#define USE_CURRENT
#endif
#endif


namespace Chroma {
  namespace LaphEnv {
   
// *************************************************************************


Handle< SystemSolver<LatticeFermion> > RR;  //renamed from PP to avoid overlap with perambulator_handler

PerambulatorInsertionHandler::FileKey::FileKey(int in_snktime, int in_srctime)
    : snk_time(in_snktime),src_time(in_srctime) {}

PerambulatorInsertionHandler::FileKey::FileKey(XmlReader& xmlr)
{
 try{
    XmlReader xmlf(xmlr,"./descendant-or-self::FileKey");
    xmlread(xmlf,"SinkTime",snk_time,"PerambulatorInsertionHandler::FileKey");
    xmlread(xmlf,"SourceTime",src_time,"PerambulatorInsertionHandler::FileKey");}
 catch(...){
    QDPIO::cerr << "Could not read PerambulatorInsertionHandler::FileKey"<<endl;
    QDP_abort(1);}
}

PerambulatorInsertionHandler::FileKey::FileKey(const FileKey& rhs)
    : snk_time(rhs.snk_time),src_time(rhs.src_time) {}

PerambulatorInsertionHandler::FileKey& PerambulatorInsertionHandler::FileKey::operator=(const FileKey& rhs)
{
 snk_time=rhs.snk_time;
 src_time=rhs.src_time;
 return *this;
}

void PerambulatorInsertionHandler::FileKey::output(XmlWriter& xmlw) const
{
 push(xmlw,"FileKey");
 write(xmlw,"SinkTime",snk_time);
 write(xmlw,"SourceTime",src_time);
 pop(xmlw);
}

bool PerambulatorInsertionHandler::FileKey::operator<(const FileKey& rhs) const
{
 return multiLessThan(src_time,rhs.src_time,snk_time,rhs.snk_time);
}

bool PerambulatorInsertionHandler::FileKey::operator==(const FileKey& rhs) const
{
 return multiEqual(src_time,rhs.src_time,snk_time,rhs.snk_time);
}

bool PerambulatorInsertionHandler::FileKey::operator!=(const FileKey& rhs) const
{
 return multiNotEqual(src_time,rhs.src_time,snk_time,rhs.snk_time);
}


// *************************************************************************




PerambulatorInsertionHandler::PerambulatorInsertionHandler()
          : uPtr(0), gSmearPtr(0), qSmearPtr(0), qactionPtr(0),
            fPtr(0), invertPtr(0), DHputPtr(0), 
            DHgetPtr(0), Nspin(4), merge_mode(false), insertPtr(0)  {}


PerambulatorInsertionHandler::PerambulatorInsertionHandler(const GaugeConfigurationInfo& gaugeinfo,
                                         const GluonSmearingInfo& gluonsmear,
                                         const QuarkSmearingInfo& quarksmear,
                                         const QuarkActionInfo& quark,
                                         const FileListInfo& flist,
                                         const string& smeared_quark_filestub,
										 const PerambulatorInsertionInfo& insertinfo,
                                         bool upper_spin_components_only,
                                         bool mergemode, const string& gauge_str)
          : invertPtr(0),  DHgetPtr(0)
{
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,flist,
          smeared_quark_filestub,insertinfo,upper_spin_components_only,gauge_str,mergemode);
}

void PerambulatorInsertionHandler::setInfo(const GaugeConfigurationInfo& gaugeinfo,
                                  const GluonSmearingInfo& gluonsmear,
                                  const QuarkSmearingInfo& quarksmear,
                                  const QuarkActionInfo& quark,
                                  const FileListInfo& flist,
                                  const string& smeared_quark_filestub,
								  const PerambulatorInsertionInfo& insertinfo,
                                  bool upper_spin_components_only,
                                  bool mergemode, const string& gauge_str)
{
 clear();
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,flist,
          smeared_quark_filestub,insertinfo,upper_spin_components_only,gauge_str,mergemode);
}


void PerambulatorInsertionHandler::set_info(const GaugeConfigurationInfo& gaugeinfo,
                                   const GluonSmearingInfo& gluonsmear,
                                   const QuarkSmearingInfo& quarksmear,
                                   const QuarkActionInfo& quark,
                                   const FileListInfo& flist,
                                   const string& smeared_quark_filestub,
								   const PerambulatorInsertionInfo& insertinfo,
                                   bool upper_spin_components_only,
                                   const string& gauge_str, bool mergemode)
{
 try{
    uPtr = new GaugeConfigurationInfo(gaugeinfo);
    gSmearPtr = new GluonSmearingInfo(gluonsmear);
    qSmearPtr = new QuarkSmearingInfo(quarksmear);
    qactionPtr = new QuarkActionInfo(quark);
    fPtr = new FileListInfo(flist);
    Nspin = (upper_spin_components_only) ? 2 : 4;
    merge_mode = mergemode;

	insertPtr = new PerambulatorInsertionInfo(insertinfo);

#if (QDP_ND == 4)
    DHputPtr=new DataPutHandlerMF<PerambulatorInsertionHandler,FileKey,UIntKey,DataType>(
                    *this,*fPtr,"Laph--QuarkPeramb","PerambulatorInsertionHandlerDataFile");
#elif (QDP_ND == 3)
    DHgetPtr=new DataGetHandlerMF<PerambulatorInsertionHandler,FileKey,UIntKey,DataType>(
                    *this,*fPtr,"Laph--QuarkPeramb","PerambulatorInsertionHandlerDataFile");
#endif
    }
 catch(...){
    QDPIO::cerr << "allocation problem in PerambulatorInsertionHandler"<<endl;
    QDP_abort(1);}

 if (!mergemode){
#if (QDP_ND == 4)
    connectGaugeConfigurationHandler(gauge_str);
#endif
    connectQuarkSmearingHandler(smeared_quark_filestub);}
}


PerambulatorInsertionHandler::~PerambulatorInsertionHandler()
{
 clear();
}


void PerambulatorInsertionHandler::clear()
{
 try{
    delete uPtr;
    delete gSmearPtr;
    delete qSmearPtr;
    delete qactionPtr;
    delete fPtr;
    delete invertPtr;

	delete insertPtr;}
 catch(...){ QDP_abort(1);}
 uPtr=0;
 gSmearPtr=0;
 qSmearPtr=0;
 qactionPtr=0;
 fPtr=0;
 invertPtr=0;

 insertPtr=0;

 if (!merge_mode){
#if (QDP_ND == 4)
    disconnectGaugeConfigurationHandler();
#endif
    disconnectQuarkSmearingHandler();}

#if (QDP_ND == 3)
 delete DHgetPtr; DHgetPtr=0;
#elif (QDP_ND == 4)
 delete DHputPtr; DHputPtr=0;
#endif
}

  // ********************************
  // *
  // *    insertion functions
  // *
  // ********************************

void PerambulatorInsertionHandler::D1Insert(LatticeFermion& chi, const LatticeFermion& psi)
{
     if(insertPtr!=0){
        std::string insert_info = insertPtr->getInsertionName();
        if(insert_info=="Twist"){
            D1InsertTwist(chi, psi);
        }else if(insert_info=="DeltaM"){
            D1InsertDeltaM(chi, psi);
        }else{
            QDP_abort(1);
        }
     }
}

void PerambulatorInsertionHandler::D1InsertTwist(LatticeFermion& chi, const LatticeFermion& psi)
{
 // int Tdir = uPtr->getTimeDir();
 
 int Tdir = 3;

 multi1d<LatticeColorMatrix> u = gaugeHandler->getData();

 multi1d<int> tw_vector = insertPtr->getTwistVector();

 chi = zero;

 if(0!=Tdir){
	chi += ( spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi), FORWARD, 0))
	- spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi), BACKWARD, 0)) )*tw_vector[0];
 }
 if(1!=Tdir){
        chi += ( spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi), FORWARD, 1))
        - spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi), BACKWARD, 1)) )*tw_vector[1];
 }
 if(2!=Tdir){
        chi += ( spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi), FORWARD, 2))
        - spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi), BACKWARD, 2)) )*tw_vector[2];
 }
 if(3!=Tdir){
        chi += ( spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi), FORWARD, 3))
        - spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi), BACKWARD, 3)) )*tw_vector[3];
 }
 QDP::Complex b = cmplx(Real(0.0), Real(-0.5));

 chi = b*chi;
}

void PerambulatorInsertionHandler::D1InsertDeltaM(LatticeFermion& chi, const LatticeFermion& psi)
{
 chi = psi;
}



  // ********************************
  // *
  // *    sub-handler connections  (private)
  // *
  // ********************************

void PerambulatorInsertionHandler::connectGaugeConfigurationHandler(const string& gauge_id)
{
 if ((gaugeCounter==0)&&(gaugeHandler.get()==0)){
    try{
       gaugeHandler.reset(new GaugeConfigurationHandler(*uPtr,gauge_id));
       gaugeCounter=1;}
    catch(...){
       QDPIO::cerr << "allocation problem in PerambulatorInsertionHandler::connectGaugeConfigurationHandler"<<endl;
       QDP_abort(1);}}
 else{
    try{
       if (gaugeHandler.get()==0) throw(std::invalid_argument("error"));
       uPtr->checkEqual(gaugeHandler->getGaugeConfigurationInfo());
       gaugeCounter++;}
    catch(...){
       QDPIO::cerr << "inconsistent PerambulatorInsertionHandler::connectGaugeConfigurationHandler"<<endl;
       QDP_abort(1);}}
}

void PerambulatorInsertionHandler::disconnectGaugeConfigurationHandler()
{
 gaugeCounter--;
 if (gaugeCounter==0){
    try{ gaugeHandler.reset();}
    catch(...){
       QDPIO::cerr << "delete problem in PerambulatorInsertionHandler::disconnectGaugeConfigurationHandler"<<endl;
       QDP_abort(1);}}
}



void PerambulatorInsertionHandler::connectQuarkSmearingHandler(const string& smeared_quark_filestub)
{
 if ((qSmearCounter==0)&&(qSmearHandler.get()==0)){
    try{
       qSmearHandler.reset(new QuarkSmearingHandler(*gSmearPtr,*uPtr,*qSmearPtr,
                                                    smeared_quark_filestub));
       qSmearCounter=1;}
    catch(...){
       QDPIO::cerr << "allocation problem in PerambulatorInsertionHandler::connectQuarkSmearingHandler"<<endl;
       QDP_abort(1);}}
 else{
    try{
       if (qSmearHandler.get()==0) throw(std::invalid_argument("error"));
       uPtr->checkEqual(qSmearHandler->getGaugeConfigurationInfo());
       qSmearHandler->updateSmearing(*qSmearPtr);  // increase eigvecs if needed
       qSmearCounter++;}
    catch(...){
       QDPIO::cerr << "inconsistent PerambulatorInsertionHandler::connectQuarkSmearingHandler"<<endl;
       QDP_abort(1);}}
}

void PerambulatorInsertionHandler::disconnectQuarkSmearingHandler()
{
 qSmearCounter--;
 if (qSmearCounter==0){
    try{ qSmearHandler.reset();}
    catch(...){
       QDPIO::cerr << "delete problem in PerambulatorInsertionHandler::disconnectQuarkSmearingHandler"<<endl;
       QDP_abort(1);}}
}



#if (QDP_ND == 4)

void PerambulatorInsertionHandler::setInverter(const InverterInfo& invinfo)
{
 try{
    delete invertPtr;
    invertPtr = new InverterInfo(invinfo);}
 catch(...){
    QDPIO::cerr << "allocation error in PerambulatorInsertionHandler::setInverter"<<endl;
    QDP_abort(1);}
}

const InverterInfo& PerambulatorInsertionHandler::getInverterInfo() const 
{
 if (invertPtr!=0){
    QDPIO::cerr << "error in PerambulatorInsertionHandler:"<<endl;
    QDPIO::cerr << "  must setInverter before calling getInverterInfo"<<endl;
    QDP_abort(1);}
 return *invertPtr;
}

void PerambulatorInsertionHandler::getFileMap(XmlWriter& xmlout) const
{
 if (isInfoSet()) DHputPtr->getFileMap(xmlout);
}

void PerambulatorInsertionHandler::outputSuffixMap()
{
 check_info_set("getSuffixMap");
 map<int,PerambulatorInsertionHandler::FileKey> suffixmap=DHputPtr->getSuffixMap();
 QDPIO::cout <<endl<<"Suffix map:"<<endl;
 for (map<int,PerambulatorInsertionHandler::FileKey>::const_iterator it=suffixmap.begin();
      it!=suffixmap.end();++it){
    QDPIO::cout << "suffix "<<it->first<<":  source time "
      << it->second.src_time<<"  sink time = "
      << it->second.snk_time << endl;}
 QDPIO::cout << endl;
}

#endif




bool PerambulatorInsertionHandler::isInfoSet() const
{
 return ((uPtr!=0)&&(gSmearPtr!=0)&&(qSmearPtr!=0)&&(fPtr!=0)
        &&(qactionPtr!=0)&&(insertPtr!=0));
}


void PerambulatorInsertionHandler::check_info_set(const string& name) const
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in PerambulatorInsertionHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling "<<name<<endl;
    QDP_abort(1);}
}


const GaugeConfigurationInfo& PerambulatorInsertionHandler::getGaugeConfigurationInfo() const 
{
 check_info_set("getGaugeConfigurationInfo");
 return *uPtr;
}

const GluonSmearingInfo& PerambulatorInsertionHandler::getGluonSmearingInfo() const
{
 check_info_set("getGluonSmearingInfo");
 return *gSmearPtr;
}

const QuarkSmearingInfo& PerambulatorInsertionHandler::getQuarkSmearingInfo() const
{
 check_info_set("getQuarkSmearingInfo");
 return *qSmearPtr;
}

const QuarkActionInfo& PerambulatorInsertionHandler::getQuarkActionInfo() const 
{
 check_info_set("getQuarkActionInfo");
 return *qactionPtr;
}

const FileListInfo& PerambulatorInsertionHandler::getFileListInfo() const 
{
 check_info_set("getFileListInfo");
 return *fPtr;
}

uint PerambulatorInsertionHandler::getNumberOfLaplacianEigenvectors() const
{
 return qSmearPtr->getNumberOfLaplacianEigenvectors();
}

int PerambulatorInsertionHandler::getTimeExtent() const 
{
 check_info_set("getTimeExtent");
 return uPtr->getTimeExtent();
}


void PerambulatorInsertionHandler::getHeader(XmlWriter& xmlout) const
{
 if (isInfoSet()){
    push(xmlout,"PerambulatorInsertionHandlerDataFile");
    uPtr->output(xmlout);
    gSmearPtr->output(xmlout);
    qSmearPtr->output(xmlout);
    qactionPtr->output(xmlout); 
    write(xmlout,"NumSpinComponents",Nspin);
	insertPtr->output(xmlout);
    pop(xmlout);}
}

bool PerambulatorInsertionHandler::checkHeader(XmlReader& xml_in, int suffix)
{
 cout << "checking header!"<<endl;
 if (xml_tag_count(xml_in,"PerambulatorInsertionHandlerDataFile")!=1) return false;
 XmlReader xmlr(xml_in,"./descendant-or-self::PerambulatorInsertionHandlerDataFile");
 GaugeConfigurationInfo gauge_check(xmlr);
 GluonSmearingInfo gsmear_check(xmlr);
 QuarkSmearingInfo qsmear_check(xmlr);
 uint numspin;
 xmlread(xmlr,"NumSpinComponents", numspin, "PerambulatorInsertionHandler");
 QuarkActionInfo qaction_check(xmlr);

 PerambulatorInsertionInfo insert_check(xmlr);

 cout << "initialized check Info classes"<<endl;
 try {
    uPtr->checkEqual(gauge_check);
 		cout << "gauge check passed"<<endl;
    gSmearPtr->checkEqual(gsmear_check);
 		cout << "gauge smear check passed"<<endl;
    qSmearPtr->checkEqual(qsmear_check); 
 		cout << "quark smear check passed"<<endl;
    if (numspin!=Nspin){
       throw(std::invalid_argument("PerambulatorInsertion checkEqual failed...NumSpinComponents mismatch"));}
 		cout << "nspin check passed"<<endl;
    qactionPtr->checkEqual(qaction_check); 
 		cout << "quark action check passed"<<endl;
	insertPtr->checkEqual(insert_check);
		cout << "insertion check passed"<<endl;
 }
 catch(const exception& xp){QDPIO::cout << xp.what()<<endl; return false;}
 return true;
}

void PerambulatorInsertionHandler::writeHeader(XmlWriter& xmlout, 
                                     const PerambulatorInsertionHandler::FileKey& fkey,
                                     int suffix)
{
 push(xmlout, "PerambulatorInsertionHandlerDataFile");
 uPtr->output(xmlout);
 gSmearPtr->output(xmlout);
 qSmearPtr->output(xmlout);
 qactionPtr->output(xmlout);
 write(xmlout,"NumSpinComponents",Nspin);
 insertPtr->output(xmlout);
 fkey.output(xmlout);
 pop(xmlout);
}


#if (QDP_ND == 4)


 // ************************************************************************************


    //  do "perambulator" inversions for all spin indices, one time source, a set of source 
    //  laph eigenvector numbers and a set of sink times.  All sink spin indices and
    //  laph_eigvec indices at the sink are computed. This routine is useful
    //  for smearing studies.
 

/*
This is just to remember which form it had.
LatticeSpinVector sv = zero;
pokeSpin(sv,temp,srcspin-1);
*/

//Sparse grid at the source 
#if 0
void PerambulatorInsertionHandler::computeSparsePerambulatorsTest() {
	if ((!isInfoSet())||(invertPtr==0)||(merge_mode)){
		QDPIO::cerr << "cannot computeSink in PerambulatorHandler until"
			<< " info and inverter set and not in merge mode"<<endl;
		QDP_abort(1);}

	int Textent = uPtr->getTimeExtent();
	int minTime = 0;
	int maxTime = Textent-1;
	int Tdir = 3;
	int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();

	QDPIO::cout <<endl<< "Quark perambulators computation for one source time,"
		<< " several source Laph eigenvector indices and several sink times beginning"<<endl;
  int src_time=1;	

	START_CODE();
	StopWatch totaltime,bulova;
	totaltime.start();
	bulova.start();
	double inittime=0.0, srctime=0.0, snktime=0.0, invtime=0.0;

	Set timeslices;                         // needed for time slice masks
	timeslices.make(TimeSlice(Tdir)); 

	SpinMatrix SrcRotate = Gamma(8) * DiracToDRMat();   // 08.05  multiply by gamma_4
	SpinMatrix SnkRotate = adj(DiracToDRMat());    // rotate back to Dirac-Pauli

#if defined(BUILD_OPENQCD12) || defined(BUILD_OPENQCD14) || defined(BUILD_OPENQCD16)
	if (invertPtr->getId().find("OPENQCD")!=std::string::npos) {

		// if using openQCD for the first time, initialize deflation subspace etc.
		if (!TheNamedObjMap::Instance().check("OpenQCDHandler")) {
			QDPIO::cout << "Setting up OpenQCD Handler" << endl;
			TheNamedObjMap::Instance().create<OpenQCDHandler>("OpenQCDHandler");
			OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<
				OpenQCDHandler>("OpenQCDHandler");
			oqcdHandler.init(uPtr->getFileName(),*qactionPtr, *invertPtr);}
		// reset the mass to make sure we are computing the correct inversions
		else {
			OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<OpenQCDHandler>("OpenQCDHandler"
					);
			oqcdHandler.resetMass(*qactionPtr);}}
	else
#else
#error "OpenQCD must be used for the inversion" 
#endif
		bulova.stop();
	inittime+=bulova.getTimeInSeconds();

	// loop over source spin

		
		int x=0; int y=2; int z=3;
		multi1d<int> coords(4); 
		coords[0]=x; coords[1]=y; coords[2]=z; coords[3]=src_time;

		for (int srccolor=0;srccolor<3;srccolor++)
			for (int srcspin=1;srcspin<=Nspin;srcspin++){

				bulova.reset();bulova.start();
			QDPIO::cout << "source point = ("<<x<<","<<y<<","<<z<<","<<src_time<<"), source color = "<<srccolor<<" and source spin = "<<srcspin<<endl<<endl;

			//  set the source field
			Real one=1.0;
			Real zzero=0.0;
			Complex ctmp = cmplx(one,zzero);  
			LatticeComplex temp = zero;
			pokeSite(temp,ctmp,coords);	
			LatticeSpinVector sv = zero;
			pokeSpin(sv,temp,srcspin-1);
			LatticeFermion latfermA = zero;
			pokeColor(latfermA,sv,srccolor);
			LatticeFermion latfermB = SrcRotate * latfermA;  // 08.05  rotate to DeGrand-Rossi, mult by gamma_4

			QDPIO::cout << "Norm of source vector = "<<sqrt(norm2(latfermB))<<endl;
			latfermA = zero;
			SystemSolverResults_t res;

#if defined(BUILD_OPENQCD12) || defined(BUILD_OPENQCD14) || defined(BUILD_OPENQCD16)
			if (invertPtr->getId().find("OPENQCD")!=std::string::npos) {
				OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<
					OpenQCDHandler>("OpenQCDHandler");
				QDPIO::cout << "source created...starting inversion with openQCD"
					<< endl;
				bulova.stop();
				srctime += bulova.getTimeInSeconds();

				// now do the inversion
				bulova.reset(); bulova.start();
				res = oqcdHandler.solve(latfermA, latfermB); // solution in latfermA
				bulova.stop();
				invtime += bulova.getTimeInSeconds();
			}
#endif
			if ((res.n_count>0)&&(res.n_count<=invertPtr->getMaxIterations())){

				latfermB = SnkRotate * latfermA;         // rotate back to Dirac-Pauli
				multi1d<DComplex> tres(Textent);

				// sink = Vs^dagger * phi
				for (int sInd=0;sInd<Nspin;sInd++){
					bulova.reset();bulova.start(); 
					LatticeColorVector phi_s = peekSpin(latfermB, sInd);

					// TODO this projection should be batched up for several solutions
					for (int nInd=0;nInd<nEigs;nInd++){
						//QDPIO::cout << "projecting onto LapH eigenvector "<<n<<endl;
						const LatticeColorVector& Ws= qSmearHandler->getLaphEigenvector(nInd);
						qSmearHandler->closeLaphLevelFiles();  // prevent too many files being open
						LatticeComplex tmp;
						tmp = localInnerProduct( Ws, phi_s );  // color contraction
						tres = sumMulti(tmp,timeslices);  // spatial sums on each timeslice
						for (int tInd=0 ; tInd<Textent;tInd++){
							QDPIO::cout << "snk_spin = "<<sInd+1<<", snk_eig = "<<nInd<<", t = "<<tInd<<", P = "<<tres[tInd]<<endl;}}
				}
			}
			else{

				QDPIO::cout << endl<<endl;
				QDPIO::cout << "Inversion FAILED to converge before max iteration reached"<<endl;
				QDPIO::cout << "Solution NOT WRITTEN to file"<<endl;
				QDPIO::cout << endl<<endl;}

		}

	totaltime.stop();
	QDPIO::cout << endl<<endl;
	QDPIO::cout << "computeQuarkPerambulators (test) ran successfully" << endl;
	QDPIO::cout << endl<<endl;
	END_CODE();
	QMP_barrier();
}
#endif 

void PerambulatorInsertionHandler::computePerambulatorsInsertion(int src_time, const set<int>& src_lapheigvec_indices,
                                               const set<int>& snk_times, bool verbose)
{
	if ((!isInfoSet())||(invertPtr==0)||(merge_mode)){
		QDPIO::cerr << "cannot computeSink in PerambulatorInsertionHandler until"
			<< " info and inverter set and not in merge mode"<<endl;
		QDP_abort(1);}

	int Textent = uPtr->getTimeExtent();
	int minTime = 0;
	int maxTime = Textent-1;
	int Tdir = 3;
	int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();

	if ((src_time<minTime)||(src_time>maxTime)){
		QDPIO::cerr << "invalid source time given to computePerambulatorsInsertion"<<endl;
		QDP_abort(1);}
	for (set<int>::const_iterator st=src_lapheigvec_indices.begin();st!=src_lapheigvec_indices.end();++st){
		if ((*st<0)||(*st>=nEigs)){
			QDPIO::cerr << "invalid source LaphEigenvector indices given to computePerambulatorsInsertion"<<endl;
			QDP_abort(1);}}
	for (set<int>::const_iterator st=snk_times.begin();st!=snk_times.end();++st){
		if ((*st<minTime)||(*st>maxTime)){
			QDPIO::cerr << "invalid sink times given to computePerambulatorsInsertion"<<endl;
			QDP_abort(1);}}

	// adjust minTime and maxTime
	minTime=src_time;
	maxTime=src_time;
	for (set<int>::const_iterator st=snk_times.begin();st!=snk_times.end();++st){
		if (*st<minTime) minTime=*st;
		if (*st>maxTime) maxTime=*st;}
	int nTimes = maxTime - minTime + 1;

	QDPIO::cout <<endl<< "Quark perambulators insertion computation for one source time,"
		<< " several source Laph eigenvector indices and several sink times beginning"<<endl;
	QDPIO::cout << " Source time = "<<src_time<<endl;

	/*
		 if (qactionPtr->getTimeBoundaryConditions()==QuarkActionInfo::ZeroDirichlet){
		 if ((src_time==0)||(src_time==(Textent-1))){
		 QDPIO::cerr << "Source time cannot be 0 to Nt-1 for Zero Dirichlet boundary conditions"<<endl;
		 return;}}
		 */
	START_CODE();
	StopWatch totaltime,bulova;
	totaltime.start();
	bulova.start();
	double inittime=0.0, srctime=0.0, snktime=0.0, invtime=0.0;

	Set timeslices;                         // needed for time slice masks
	timeslices.make(TimeSlice(Tdir)); 

	Set centLat;
	centLat.make(TimeInterval(minTime, nTimes, Tdir));

	SpinMatrix SrcRotate = Gamma(8) * DiracToDRMat();   // 08.05  multiply by gamma_4
	SpinMatrix SnkRotate = adj(DiracToDRMat());    // rotate back to Dirac-Pauli
	/*
		 string fermact_xml = qactionPtr->getDescription();
		 string fermact_id = qactionPtr->getActionName();

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	GroupXML_t solverInfo;
	solverInfo.xml =  invertPtr->output();
	solverInfo.id = invertPtr->getId();
	solverInfo.path = "//InvertParam";

	// Initialize fermion action

	istringstream xml_s(fermact_xml);
	XMLReader fermacttop0(xml_s);
	XMLReader fermacttop(fermacttop0,"./descendant-or-self::Description");   // due to XmlReader bug

	Handle< FermionAction<T,P,Q> > S_f; 
	Handle< FermState<T,P,Q> > state;
	Handle< SystemSolver<LatticeFermion> > RR;
	*/
#if defined(BUILD_OPENQCD12) || defined(BUILD_OPENQCD14) || defined(BUILD_OPENQCD16)
	if (invertPtr->getId().find("OPENQCD")!=std::string::npos) {

		// if using openQCD for the first time, initialize deflation subspace etc.
		if (!TheNamedObjMap::Instance().check("OpenQCDHandler")) {
			QDPIO::cout << "Setting up OpenQCD Handler" << endl;
			TheNamedObjMap::Instance().create<OpenQCDHandler>("OpenQCDHandler");
			OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<
				OpenQCDHandler>("OpenQCDHandler");
			oqcdHandler.init(uPtr->getFileName(),*qactionPtr, *invertPtr);}
		// reset the mass to make sure we are computing the correct inversions
		else {
			OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<OpenQCDHandler>("OpenQCDHandler"
					);
			oqcdHandler.resetMass(*qactionPtr);}}
	else
#else
#error "OpenQCD must be used for the inversion" 
#endif
		bulova.stop();
	inittime+=bulova.getTimeInSeconds();

	// loop over source eigenvector number and source spin

	int nsrc=src_lapheigvec_indices.size()*Nspin;
	int srccount=0;
	for (set<int>::const_iterator st=src_lapheigvec_indices.begin();st!=src_lapheigvec_indices.end();++st)
		for (int srcspin=1;srcspin<=Nspin;srcspin++){
			int srcev_ind=*st;

			bulova.reset();bulova.start();
			QDPIO::cout <<endl<< "Starting src number "<<srccount<<" (with "<<nsrc-1<<" as last)"<<endl;
			QDPIO::cout << "Source lapheigvec index "<<srcev_ind<<" and source spin = "<<srcspin<<endl<<endl;
			srccount++;

			bool doneflag=true;
			for (set<int>::const_iterator it=snk_times.begin();it!=snk_times.end();++it){
				FileKey fkey(*it,src_time);
				DHputPtr->open(fkey);
				if (!DHputPtr->queryData(UIntKey((srcspin-1)*nEigs+srcev_ind))){
					doneflag=false; DHputPtr->close(); break;}
				DHputPtr->close();}
			if (doneflag){
				QDPIO::cout << "warning: these quark perambulators insertions already computed..."
					<< "skip re-computing"<<endl;
				continue;}

			//  set the source field
			Real one=1.0;
			Real zzero=0.0;
			LatticeComplex temp = zero;
			temp[timeslices[src_time]] = cmplx(one,zzero);  // temp has value 1 on src_time time slice, zero everywhere else
			LatticeSpinVector sv = zero;
			pokeSpin(sv,temp,srcspin-1);
			const LatticeColorVector& Vs = qSmearHandler->getLaphEigenvector(srcev_ind);
			LatticeFermion latfermA = zero;
			latfermA[timeslices[src_time]] = sv * Vs;
			LatticeFermion latfermB = SrcRotate * latfermA;  // 08.05  rotate to DeGrand-Rossi, mult by gamma_4

			QDPIO::cout << "Norm of source vector = "<<sqrt(norm2(latfermB))<<endl;
			latfermA = zero;
			SystemSolverResults_t res;
			SystemSolverResults_t res2;

#if defined(BUILD_OPENQCD12) || defined(BUILD_OPENQCD14) || defined(BUILD_OPENQCD16)
			if (invertPtr->getId().find("OPENQCD")!=std::string::npos) {
				OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<
					OpenQCDHandler>("OpenQCDHandler");
				QDPIO::cout << "source created...starting inversion with openQCD"
					<< endl;
				bulova.stop();
				srctime += bulova.getTimeInSeconds();

				// now do the inversion
				bulova.reset(); bulova.start();
				res = oqcdHandler.solve(latfermA, latfermB); // solution in latfermA
				latfermB = zero;
       			D1Insert(latfermB, latfermA);
       			latfermA = zero;
       			res2 = oqcdHandler.solve(latfermA, latfermB); // solution in latfermA
				bulova.stop();
				invtime += bulova.getTimeInSeconds();
			}
#endif
			if ((res.n_count>0)&&(res.n_count<=invertPtr->getMaxIterations())&&(res2.n_count>0)&&(res2.n_count<=invertPtr->getMaxIterations())){

				latfermB = SnkRotate * latfermA;         // rotate back to Dirac-Pauli
				 multi1d<int> coords(4);
				 coords[0] = 4;
				 coords[1] = 0;
				 coords[2] = 0;
				 coords[3] = 5;
				 Fermion ftmp = peekSite(latfermB, coords);
				 for (int cInd=0;cInd<3;cInd++)  
					 for (int sInd=0;sInd<4;sInd++) { 
					 ColorVector ctemp = peekSpin(ftmp, sInd);
					 Complex tmp = peekColor(ctemp, cInd);
				 QDPIO::cout << "S_________________P"<<endl;
				 QDPIO::cout << "SrcT = " << src_time<<", SrcSpin = " << srcspin <<", Srcev_ind = "<<srcev_ind<<
					", cInd = "<<cInd<<", SnkSpin = "<<sInd+1<<", P = "<< tmp <<endl;
				 QDPIO::cout << "S________________P"<<endl;       
					 }	
				 multi2d<Complex> quark_sink(nTimes,nEigs*Nspin);
				 multi1d<DComplex> tres(Textent);

				// sink = Vs^dagger * phi

				for (int s=0;s<Nspin;s++){
					bulova.reset();bulova.start(); 
					LatticeColorVector phi_s = peekSpin(latfermB, s);
					/*
						 int Lx,Ly,Lz,Lt; 
						 Lx=Layout::lattSize()[0]; 
						 Ly=Layout::lattSize()[1]; 
						 Lz=Layout::lattSize()[2];
						 Lt=Layout::lattSize()[3];
						 for (int t=0; t<Lt; t++) { 
						 QDPIO::cout << "src_ev = "<<srcev_ind<<", spin = "<<s+1<<", time = "<<t<<endl;
						 for (int z=0; z<Lz; z++)  
						 for (int y=0; y<Ly; y++) 
						 for (int x=0; x<Lx; x++) { 
						 multi1d<int> coords(4);
						 coords[0]=x;coords[1]=y;coords[2]=z;coords[3]=t;
						 ColorVector ctemp = peekSite(phi_s, coords);
						 for (int c=0;c<3;c++){
						 Complex tmp2 = peekColor(ctemp,c);
						 QDPIO::cout << "(x,y,z) = ("<<x<<", "<<y<<", "<<z<<"), color = "<<c<<", res = ("<<real(tmp2)<<", "<<imag(tmp2)<<")"<<endl;
						 }
						 }
						 }
						 */

					// TODO this projection should be batched up for several solutions
					for (int n=0;n<nEigs;n++){
						//QDPIO::cout << "projecting onto LapH eigenvector "<<n<<endl;
						const LatticeColorVector& Ws= qSmearHandler->getLaphEigenvector(n);
						qSmearHandler->closeLaphLevelFiles();  // prevent too many files being open
						LatticeComplex tmp;
						if (nTimes==Textent)
							tmp = localInnerProduct( Ws, phi_s );  // color contraction
						else
							tmp[centLat[0]] = localInnerProduct(Ws, phi_s); // color contraction
						tres = sumMulti(tmp,timeslices);  // spatial sums on each timeslice
						for (set<int>::const_iterator it=snk_times.begin();it!=snk_times.end();++it){
							quark_sink(*it-minTime,n+s*nEigs)=tres[*it];}}
					bulova.stop(); snktime+=bulova.getTimeInSeconds();
				}

				// output to file
				for (set<int>::const_iterator it=snk_times.begin();it!=snk_times.end();++it){
					FileKey fkey(*it,src_time);
					DHputPtr->open(fkey);
					UIntKey rkey((srcspin-1)*nEigs+srcev_ind);
					if (!DHputPtr->queryData(rkey)){
						DHputPtr->putData(rkey,quark_sink[*it-minTime]);
						DHputPtr->flush();}
					if (verbose){
						QDPIO::cout << "Compare Standard Perambulators Insertion! " << endl;
						QDPIO::cout << "Coefficients for source time "<<src_time<<", src_spin "<<
							srcspin<<", src_ev "<<srcev_ind<<", and sink time "<<*it<<endl;
						for (int s=0;s<Nspin;s++){
							for (int n=0;n<nEigs;n++){
								QDPIO::cout << "coef for sink spin = "<<s+1<<" and eigenlevel "<<n<<" = "
									<< quark_sink(*it-minTime,s*nEigs+n)<<endl;}}}}
			}

			else{

				QDPIO::cout << endl<<endl;
				QDPIO::cout << "Inversion FAILED to converge before max iteration reached"<<endl;
				QDPIO::cout << "Solution NOT WRITTEN to file"<<endl;
				QDPIO::cout << endl<<endl;}

		}

	totaltime.stop();
	QDPIO::cout << endl<<endl;
	QDPIO::cout << "computeQuarkPerambulatorsInsertion ran successfully" << endl;
	QDPIO::cout << "                 Total time = "<<totaltime.getTimeInSeconds() << " seconds" << endl;
	QDPIO::cout << "        Initialization time = "<<inittime<<" seconds"<<endl;
	QDPIO::cout << "   Total source set up time = "<<srctime<<" seconds"<<endl;
	QDPIO::cout << "     Total time in inverter = "<<invtime<<" seconds"<<endl;
	QDPIO::cout << " Total sink completion time = "<<snktime<<" seconds"<<endl;
	QDPIO::cout << endl<<endl;
	END_CODE();
	QMP_barrier();
}

 // ************************************************************************************

#if 0
void PerambulatorHandler::obtainPhi(LatticeFermion& latfermA, const LatticeFermion& latfermB){
QDPIO::cout << "Norm of source vector = "<<sqrt(norm2(latfermB))<<endl;
SystemSolverResults_t res;
#if defined(BUILD_OPENQCD12) || defined(BUILD_OPENQCD14) || defined(BUILD_OPENQCD16)
if (invertPtr->getId().find("OPENQCD")!=std::string::npos) {
   OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<
           OpenQCDHandler>("OpenQCDHandler");
   QDPIO::cout << "source created...starting inversion with openQCD"
           << endl;
   res = oqcdHandler.solve(latfermA, latfermB);
   }
#else 
  #error "OpenQCD must be used for the inversion" 
#endif
}
void PerambulatorHandler::computeGPerambulators(const set< pair<pair<int, int>, set<int>> >& Snk_Src_Ins_times, const set<string>& current_components, const set<int>& src_spin, const set<int>& snk_spin, bool verbose) {
	if ((!isInfoSet())||(invertPtr==0)||(merge_mode)){
		QDPIO::cerr << "cannot computeSink in PerambulatorHandler until"<< " info and inverter set and not in merge mode"<<endl;
		QDP_abort(1);}
	int Textent = uPtr->getTimeExtent();
	int minTime = 0;
	int maxTime = Textent-1;
	int Tdir = 3;
	int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();

	size_t crnts = current_components.size();
	size_t GPnr = Snk_Src_Ins_times.size() *  crnts * nEigs * nEigs * 4 * 4;

	START_CODE();
	StopWatch totaltime,bulova;
	totaltime.start();
	bulova.start();
	double inittime=0.0, srctime=0.0, snktime=0.0, invtime=0.0, triplet_timing;
	Set timeslices;                         // needed for time slice masks
	timeslices.make(TimeSlice(Tdir));

	SpinMatrix Id = 1;
	SpinMatrix SrcRotate = Gamma(8) * DiracToDRMat();   //  multiply by gamma_4
	SpinMatrix SnkRotate = adj(DiracToDRMat());    // rotate back to Dirac-Pauli
	SpinMatrix  Src2Rotate = Gamma(15) * DiracToDRMat();
	SpinMatrix Snk2Rotate = Gamma(15) * Id;

#if defined(BUILD_OPENQCD12) || defined(BUILD_OPENQCD14) || defined(BUILD_OPENQCD16)
	if (invertPtr->getId().find("OPENQCD")!=std::string::npos) {
		// if using openQCD for the first time, initialize deflation subspace etc.
		if (!TheNamedObjMap::Instance().check("OpenQCDHandler")) {
			QDPIO::cout << "Setting up OpenQCD Handler" << endl;
			TheNamedObjMap::Instance().create<OpenQCDHandler>("OpenQCDHandler");
			OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<
				OpenQCDHandler>("OpenQCDHandler");
			oqcdHandler.init(uPtr->getFileName(),*qactionPtr, *invertPtr);}
		// reset the mass to make sure we are computing the correct inversions
		else {
			OpenQCDHandler& oqcdHandler = TheNamedObjMap::Instance().getData<OpenQCDHandler>("OpenQCDHandler");
			oqcdHandler.resetMass(*qactionPtr);}}
	else
#else   
#error "OpenQCD must be used for the inversion"
#endif

		bulova.stop();
	inittime+=bulova.getTimeInSeconds();

	int srccount=0;
	QDPIO::cout<< " "<<endl;
	QDPIO::cout<<"The calculations are performed and structured as follows:"<<endl;
	QDPIO::cout<<"****************** Loop Through Source Spin ******************"<<endl;
	QDPIO::cout<<"*************** Loop Through Sink Spin ***************"<<endl;
	QDPIO::cout<<"*********** Loop Throguh Sink Laphev ************"<<endl;
	QDPIO::cout<<" "<<endl;
	QDPIO::cout<<"Solve Dirac Equation for"<<endl;
	QDPIO::cout<<"given sink time, spin"<<endl;
	QDPIO::cout<<"and Laphev"<<endl;
	QDPIO::cout<<" "<<endl;
	QDPIO::cout<<"********* Loop Through TimeTriplet *********"<<endl;
	QDPIO::cout<<"***** Loop Throguh Source Laphev ******"<<endl;
	QDPIO::cout<<" "<<endl;
	QDPIO::cout<<"Solve Dirac Equation for"<<endl;
	QDPIO::cout<<"given source time, spin"<<endl;
	QDPIO::cout<<"and Laphev"<<endl;
	QDPIO::cout<<" "<<endl;
	QDPIO::cout<<"**** Loop Through Currents ****"<<endl;
	QDPIO::cout<<"Perform localInnerProduct"<<endl;
	QDPIO::cout<<"and save the results"<<endl;
	QDPIO::cout<<"All Loops end here"<<endl;
	QDPIO::cout<<" "<<endl;
	QDPIO::cout<<"** In the given computation set there are  "<< GPnr <<" GPerambulators"<<endl;
	QDPIO::cout<<"**  computeGPerambulators starts running now  **"<<endl;
	QDPIO::cout<<" "<<endl;
	SpinMatrix GMatrix = Id*Gamma(7);
	SpinMatrix NGMatrix =  SrcRotate * GMatrix;
	int GPCounter = 0;

	for (set<int>::const_iterator it1=src_spin.begin();it1!=src_spin.end();++it1){
		// Loop Through Source Spin
		int delta0=*it1;
		for (set<int>::const_iterator it2=snk_spin.begin();it2!=snk_spin.end();++it2){
			// Loop Through Sink Spin
			int delta=*it2 ;
			//Construct the spin vector S^{snk} = (Gamma(delta, 0), Gamma(delta, 1), Gamma(delta, 2), Gamma(delta, 3))

			for (int snkev_ind=0;snkev_ind<nEigs;snkev_ind++){
				// Loop Through Sink Eigenvectors	
				for (const auto& triplet : Snk_Src_Ins_times) {
					// Loop Through Set ( Snk_time, Src_time, (Ins_times) )
					const auto& snk_src_pair = triplet.first;
					const auto& ins_times = triplet.second;
					int snk_time = snk_src_pair.first;
					int src_time = snk_src_pair.second;
					LatticeSpinVector s_snk = zero;
					Real one_t=1.0;
					Real zzero_t=0.0;
					LatticeComplex temp1_t = zero;
					temp1_t[timeslices[snk_time]] = cmplx(one_t,zzero_t);
					pokeSpin(s_snk,temp1_t, delta);
					s_snk = DiracToDRMat() * s_snk;
					//s_snk = Src2Rotate * s_snk;
					if ((snk_time<minTime)||(snk_time>maxTime)){
						QDPIO::cerr << "invalid sink times given to computePerambulators"<<endl;
						QDP_abort(1);} 
					triplet_timing=0.0;
					bulova.reset();
					bulova.start();
					if ((src_time<minTime)||(src_time>maxTime)){
						QDPIO::cout << minTime << " " << maxTime << endl;
						QDPIO::cerr << "invalid source time given to computePerambulators"<<endl;
						QDP_abort(1);}

					int nTimes = maxTime - minTime + 1;
					//1.Solve
					LatticeFermion Phi_tilde=zero;
					LatticeFermion Phi_tilde_temp;
					const LatticeColorVector& Vsnk = qSmearHandler->getLaphEigenvector(snkev_ind);
				  Phi_tilde[timeslices[snk_time]] = s_snk * Vsnk;
					PerambulatorHandler::obtainPhi(Phi_tilde_temp, Phi_tilde);
					Phi_tilde = Snk2Rotate *  Phi_tilde_temp; 
					for (int srcev_ind=0;srcev_ind<nEigs;srcev_ind++){
						// Loop Through Source Eigenvectors
						//.Solve
						// **************************** Correct ****************************
						LatticeSpinVector s_src = zero;
						Real one=1.0;
						Real zzero=0.0;
						LatticeComplex temp1 = zero;
						temp1[timeslices[src_time]] = cmplx(one,zzero);
						QDPIO::cout << "delta0 = " << delta0<<endl;
						QDPIO::cout << "delta = " << delta<<endl;
						pokeSpin(s_src,temp1, delta0);
						s_src = SrcRotate * s_src;
						const LatticeColorVector& Vsrc = qSmearHandler->getLaphEigenvector(srcev_ind);
						LatticeFermion Phi_src_k0_temp;
						LatticeFermion Phi_src_k0=zero;
						Phi_src_k0[timeslices[src_time]] = s_src * Vsrc;
						PerambulatorHandler::obtainPhi(Phi_src_k0_temp, Phi_src_k0);
						PerambulatorHandler::obtainPhi(Phi_src_k0, Phi_src_k0_temp);
						//Phi_src_k0 = SnkRotate * Phi_src_k0_temp;
						int Current_Counter = 0;
						int Jsize = current_components.size();
						// ****************************         ****************************
						//Phi_src_k0 = Phi_src_k0_temp;
						/*
						 * I get the same resul for the standard perambulator!
						// **************************** Add to Correct the suggestions by John **************************** 
						Phi_src_k0_temp =  SrcRotate * Phi_src_k0;
						PerambulatorHandler::obtainPhi(Phi_src_k0_temp, src_time, s_src, srcev_ind);
						Phi_src_k0 = SnkRotate * Phi_src_k0_temp;
						*/ 


						for (set<string>::const_iterator it3=current_components.begin();it3!=current_components.end();++it3){
							// Loop Through Current Components
							string currentMatrix =* it3;
							SpinMatrix jM;
							if(currentMatrix == "V1"){
								jM= Id*Gamma(1);}
							else if(currentMatrix == "V0"){
								jM = Id;}
							else if (currentMatrix == "V2"){
								jM = Id*Gamma(2);}
							else if (currentMatrix == "V3"){
								jM = Id*Gamma(4);}
							else if (currentMatrix == "V4"){
								jM = Id*Gamma(8);}
							else if (currentMatrix == "PV1"){
								jM = Id*Gamma(14);}
							else if (currentMatrix == "PV2"){
								jM = Id*Gamma(13);}
							else if (currentMatrix == "PV3"){
								jM = Id*Gamma(11);}
							else if (currentMatrix == "PV4"){
								jM = Id*Gamma(7);}
							else if (currentMatrix == "PS5"){
								jM = Id*Gamma(15);}
							else{QDPIO::cerr << "error in the GeneralizedPerambulator"<<endl;
								QDPIO::cerr << "The current argument must be in {V0, V1, V2, V3, V4, PV1, PV2,PV3, PV4, PS}"<<endl;
								QDP_abort(1);}
							// SpinMatrix GJ = GMatrix * jM; 
							SpinMatrix GJ = jM;
							multi1d<DComplex> tres(Textent);
							QDPIO::cout<<  "----------------------	"
								<<"Do local innerproduct to obtain the Gperambulator number  "
								<<  GPCounter<< "----------------------  "
								<< endl;
							Complex res_i;
							//This part, from here                               
							const LatticeColorVector& Vs = qSmearHandler->getLaphEigenvector(snkev_ind);
							LatticeFermion latfermB = zero;
							latfermB[timeslices[snk_time]] = s_snk * Vs;
							res_i = innerProduct(latfermB, Phi_src_k0);
	
   
							//res_i = innerProduct( Phi_tilde, Phi_src_k0 );
							QDPIO::cout << "For: srcTime"<<src_time<<"_snkTime"<<snk_time<<endl;
							QDPIO::cout <<  "delta0: "<< delta0 <<endl;
							QDPIO::cout <<  "delta: "<< delta <<endl;
							QDPIO::cout <<  "snkev_ind: "<< snkev_ind <<endl;
							QDPIO::cout <<  "srcev_ind: " << srcev_ind <<endl;   
							QDPIO::cout << "The summed GP is: " <<endl;
							QDPIO::cout << res_i<<endl;
							// until here, is extra! 
							LatticeComplex tmp = zero;
							multi1d<LatticeColorVector> Color_V_R(4);
							multi1d<LatticeColorVector> Color_V_L(4);
							for(int s = 0; s<4; s++){
								Color_V_R[s] = peekSpin(GJ * Phi_src_k0, s);
								Color_V_L[s] = peekSpin(Phi_tilde, s);
							}
							for(int s = 0; s<4; s++){
								tmp += localInnerProduct( Color_V_L[s],  Color_V_R[s]);
							} 


							//	tmp = localInnerProduct( Phi_tilde, GJ * Phi_src_k0 );
							//16.05		const LatticeColorVector& Vs = qSmearHandler->getLaphEigenvector(snkev_ind);
							//16.05               int Tdir = 3;                   
							//16.05                LatticeFermion latfermB = zero;
							//16.05                Set timeslices;                 
							//16.05                timeslices.make(TimeSlice(Tdir));
							//16.05                latfermB[timeslices[snk_time]] = s_snk * Vs;
							/* 
							 * succeeded
							 multi1d<int> coords(4);
							 coords[0] = 0;
							 coords[1] = 0;
							 coords[2] = 0;
							 coords[3] = 0; 
							 Fermion ftmp = peekSite(Phi_src_k0, coords);
							 ColorVector ctemp = peekSpin(ftmp, 0);
							 Complex tmp_1 = peekColor(ctemp, 0);
							 QDPIO::cout << "__________________"<<endl;
							 QDPIO::cout << "At SrcT: " << src_time<<endl;
							 QDPIO::cout << "and SrcSpin: " << delta0 <<endl;
							 QDPIO::cout << "and Srcev_ind:  "<< srcev_ind<<endl;
							 QDPIO::cout << "tmp = "<< tmp_1 <<endl;
							 QDPIO::cout << "__________________"<<endl;
							 */
							//16.05		LatticeColorVector stemp = peekSpin(GJ * Phi_src_k0, delta);	
							//16.05		tmp = localInnerProduct( Vs,  stemp);
							// Hier         latfermB is a LatticeFermion at snk_time
							// While in SP  latfermB is a LatticeColorVector&
							tres = sumMulti(tmp, timeslices);//This does the spatial sum on all time slices	
							multi1d<Complex> innerproduct_on_instime(ins_times.size());
							for(set<int>::const_iterator instime = ins_times.begin(); instime != ins_times.end(); ++instime){
								if (*instime >= 0 && *instime < tres.size()) {
									innerproduct_on_instime[*instime-*ins_times.begin()] = tres[*instime];}
								else {
									QDPIO::cerr << "ERROR: instime " << *instime << " is out of bounds!" << endl;}
							}
							QDPIO::cout <<"Inner product is performed. Now save the results"<<endl;

							FileKey fkey(snk_time,src_time);
							DHputPtr->open(fkey);

							UIntKey rkey(delta*(Nspin * nEigs * nEigs * Jsize)
									+ delta0 *( nEigs * nEigs * Jsize)
									+ srcev_ind *(nEigs * Jsize) 
									+ snkev_ind *( Jsize ) +
									Current_Counter);


							std::cout <<  "delta0"<< delta0 <<std::endl;
							std::cout <<  "delta"<< delta <<std::endl;
							std::cout <<  "Current_Counter"<< Current_Counter <<std::endl;
							std::cout <<  "Laph_Sink_Counter"<< snkev_ind <<std::endl;
							std::cout <<  "Laph_Source_Counter" << srcev_ind <<std::endl;	
							std::cout << "current file key is (" << fkey.snk_time << "," << fkey.src_time << ")" << std::endl; 
							std::cout << "current record key value is " << rkey.getValue() << std::endl;

							if (!DHputPtr->queryData(rkey)){
								DHputPtr->putData(rkey,innerproduct_on_instime);
								DHputPtr->flush();}
							GPCounter++;
							QDPIO::cout <<" Saving done!"<<endl;
							QDPIO::cout << "********************************************************" <<endl;

							//Nspin, nEigs, Jsize
							Current_Counter += 1;}// End Of Loop Through Current Components
					}// End Of Loop Through Source Eigenvectors

					triplet_timing+=bulova.getTimeInSeconds();

				}// End Of Loop Through Set ( Snk_time, Src_time, (Ins_times) )
			}// End Of Loop Through Sink Eigenvectors
		}// End Of Loop Through Sink Spin
	}// End Of Loop Through Source Spin
	totaltime.stop();
	QDPIO::cout << endl<<endl;
	QDPIO::cout << "computeQuarkPerambulators ran successfully" << endl;
	QDPIO::cout << "                 Total time = "<<totaltime.getTimeInSeconds() << " seconds" << endl;
	QDPIO::cout << "        Initialization time = "<<inittime<<" seconds"<<endl;
	QDPIO::cout << "   Total source set up time = "<<srctime<<" seconds"<<endl;
	QDPIO::cout << "     Total time in inverter = "<<invtime<<" seconds"<<endl;
	QDPIO::cout << " Total sink completion time = "<<snktime<<" seconds"<<endl;
	QDPIO::cout << endl<<endl;
	END_CODE();
	QMP_barrier();
}//End of the method: PerambulatorHandler::computeGPerambulators
#endif

// ***************************************************************************


void PerambulatorInsertionHandler::mergeData(const FileListInfo& input_files)
{
	if (!merge_mode) return;
	QDPIO::cout <<endl<<endl<<"Performing merge of perambulator insertion data files"<<endl;
	QDPIO::cout << "NumSpinComponents = "<<Nspin<<endl;

	for (int suffix=input_files.getMinFileNumber(); suffix<=input_files.getMaxFileNumber(); ++suffix){

		string infilename(input_files.getFileName(suffix));
		if (!fileExists(infilename)) continue;
		QDPIO::cout <<endl<< "Processing file "<<infilename<<endl;
		string header;
		IOMap<UIntKey,multi1d<Complex>> in_iom;
		try{
			in_iom.openReadOnly(infilename,"Laph--QuarkPeramb",header);}
		catch(...){continue;}
		XmlReader xmlh(header);
		bool flag=false;
		if (xml_tag_count(xmlh,"PerambulatorInsertionHandlerDataFile")==1){
			XmlReader xmlr(xmlh,"./descendant-or-self::PerambulatorInsertionHandlerDataFile");
			GaugeConfigurationInfo gauge_check(xmlr);
			GluonSmearingInfo gsmear_check(xmlr);
			QuarkSmearingInfo qsmear_check(xmlr);
			uint numspin;
			xmlread(xmlr,"NumSpinComponents", numspin, "PerambulatorInsertionHandler");
			QuarkActionInfo qaction_check(xmlr);
			PerambulatorInsertionInfo insert_check(xmlr);
			try {
				uPtr->checkEqual(gauge_check);
				gSmearPtr->checkEqual(gsmear_check);
				qSmearPtr->checkEqual(qsmear_check); 
				if (numspin!=Nspin){
					throw(std::invalid_argument("Perambulator insertion checkEqual failed...NumSpinComponents mismatch"));}
				qactionPtr->checkEqual(qaction_check); 
				insertPtr->checkEqual(insert_check);
				flag=true;}
			catch(const std::exception& xp){QDPIO::cout << xp.what()<<endl;}}
		if (!flag){
			QDPIO::cout << "   Header did not pass so skipping this file"<<endl;
			continue;}

		XmlReader xmlr(xmlh,"./descendant-or-self::PerambulatorInsertionHandlerDataFile");
		FileKey fkey(xmlr);
		std::set<UIntKey> keys;
		in_iom.getKeys(keys);
		uint count=0;
		DHputPtr->open(fkey);
		for (std::set<UIntKey>::iterator it=keys.begin();it!=keys.end();++it){
			multi1d<Complex> buf; 
			in_iom.get(*it,buf);
			if (!DHputPtr->queryData(*it)){   // new data so output
				DHputPtr->putData(*it,buf); ++count;}
			else{
				QDPIO::cout << "Encountered duplicate data for record key "<<it->getValue()
					<<" .. skipped non-first occurrence"<<endl;}}
		DHputPtr->flush();
		in_iom.close();
		QDPIO::cout << "Finished: processed "<<count<<" record keys"<<endl;}
}


// *******************************************************************

#elif (QDP_ND == 3)


void PerambulatorInsertionHandler::getSourceTimes(set<int>& source_times) const
{
	source_times.clear();
	check_info_set("getSourceTimes");
	set<FileKey> fkeys=DHgetPtr->getFileKeys();
	for (set<FileKey>::const_iterator it=fkeys.begin();it!=fkeys.end();++it){
		source_times.insert(it->src_time);}
}



multi2d<Complex> PerambulatorInsertionHandler::getData(int snk_time, int snk_spin, 
		int src_time, int src_spin, int nEigsUse)
{
	uint nEigs=getNumberOfLaplacianEigenvectors();
	int minTime = 0;
	int maxTime = uPtr->getTimeExtent()-1;
	if ((src_time<minTime)||(src_time>maxTime)||(snk_time<minTime)||(snk_time>maxTime)
			||(src_spin<1)||(src_spin>Nspin)||(snk_spin<1)||(snk_spin>Nspin)
			||(nEigsUse<1)||(nEigsUse>nEigs)){
		QDPIO::cout << "Invalid parameters passed to getData in PerambulatorHandler"<<endl;
		QDP_abort(1);}
	multi2d<Complex> result(nEigsUse,nEigsUse);
	try{
		FileKey fkey(snk_time,src_time);
		multi1d<Complex> buff;
		for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
			UIntKey rkey((src_spin-1)*nEigs+srcev_ind);
			buff=DHgetPtr->getData(fkey,rkey);
			for (int n=0;n<nEigsUse;n++){
				result(n,srcev_ind)=buff[(snk_spin-1)*nEigs+n];}}}
	catch(...){ QDP_abort(1);}
	return result;
}


// *****************************************************


bool PerambulatorInsertionHandler::queryData(int src_time, int src_spin, int snk_time, 
		int snk_spin, int nEigsUse)
{
	uint nEigs=getNumberOfLaplacianEigenvectors();
	int minTime = 0;
	int maxTime = uPtr->getTimeExtent()-1;
	if ((src_time<minTime)||(src_time>maxTime)||(snk_time<minTime)||(snk_time>maxTime)
			||(src_spin<1)||(src_spin>Nspin)||(snk_spin<1)||(snk_spin>Nspin)){
		return false;}
	FileKey fkey(src_time,snk_time);
	for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
		UIntKey rkey((src_spin-1)*nEigs+srcev_ind);
		if (!DHgetPtr->queryData(fkey,rkey)) return false;}
	return true;
}


void PerambulatorInsertionHandler::clearGaugeData()
{
	if (qSmearHandler.get()) qSmearHandler->clearLaphEigenvectors();
}


#endif


// ***************************************************************

//  static pointers (set to null in default constructor)

unique_ptr<QuarkSmearingHandler> PerambulatorInsertionHandler::qSmearHandler;
unique_ptr<GaugeConfigurationHandler> PerambulatorInsertionHandler::gaugeHandler;

int PerambulatorInsertionHandler::qSmearCounter=0;
int PerambulatorInsertionHandler::gaugeCounter=0;


// ***************************************************************
}
}

