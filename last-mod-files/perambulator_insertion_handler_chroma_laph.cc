#include "perambulator_insertion_handler_chroma_laph.h"
#include "multi_compare.h"
#include<array>

using namespace std;

namespace LaphEnv {
   
PerambulatorInsertionHandlerCL::FileKey::FileKey(XMLHandler& xmlr)
{
    XMLHandler xmlf(xmlr,"FileKey");
    xmlread(xmlf,"SourceTime",src_time,"PerambulatorInsertionHandler::FileKey");
    xmlread(xmlf,"SinkTime",snk_time,"PerambulatorInsertionHandler::FileKey");
}

PerambulatorInsertionHandlerCL::FileKey& PerambulatorInsertionHandlerCL::FileKey::operator=(const FileKey& rhs)
{
 src_time=rhs.src_time;
 snk_time=rhs.snk_time;
 return *this;
}

void PerambulatorInsertionHandlerCL::FileKey::output(XMLHandler& xmlw) const
{
 xmlw.set_root("FileKey");
 xmlw.put_child("SourceTime",make_string(src_time));
 xmlw.put_child("SinkTime",make_string(snk_time));
}

bool PerambulatorInsertionHandlerCL::FileKey::operator<(const FileKey& rhs) const
{
 return multiLessThan(src_time,rhs.src_time,snk_time,rhs.snk_time);
}

bool PerambulatorInsertionHandlerCL::FileKey::operator==(const FileKey& rhs) const
{
 return multiEqual(src_time,rhs.src_time,snk_time,rhs.snk_time);
}

bool PerambulatorInsertionHandlerCL::FileKey::operator!=(const FileKey& rhs) const
{
 return multiNotEqual(src_time,rhs.src_time,snk_time,rhs.snk_time);
}


// *************************************************************************

PerambulatorInsertionHandlerCL::PerambulatorInsertionHandlerCL()
          : uPtr(0), gSmearPtr(0), qSmearPtr(0), qActionPtr(0), fPtr(0), DHgetPtr(0), insertPtr(0)  {}


PerambulatorInsertionHandlerCL::PerambulatorInsertionHandlerCL(const GaugeConfigurationInfo& gaugeinfo,
                                         	 const GluonSmearingInfo& gluonsmear,
											 const QuarkSmearingInfo& quarksmear,
                       						 const QuarkActionInfo& quark,
											 const FileListInfo& flist,
											 const PerambulatorInsertionInfo& insertinfo,
											 bool upper_spin_only
																				 )
          : DHgetPtr(0)
{
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,flist,insertinfo,upper_spin_only);
}

void PerambulatorInsertionHandlerCL::setInfo(const GaugeConfigurationInfo& gaugeinfo,
                                  const GluonSmearingInfo& gluonsmear,
                                  const QuarkSmearingInfo& quarksmear,
                                  const QuarkActionInfo& quark,
								  const FileListInfo& flist,
								  const PerambulatorInsertionInfo& insertinfo,
								  bool upper_spin_only
																	)
{
 clear();
 set_info(gaugeinfo,gluonsmear,quarksmear,quark,flist,insertinfo,upper_spin_only);
}


void PerambulatorInsertionHandlerCL::set_info(const GaugeConfigurationInfo& gaugeinfo,
                                   const GluonSmearingInfo& gluonsmear,
                                   const QuarkSmearingInfo& quarksmear,
                				   const QuarkActionInfo& quark,
 								   const FileListInfo& flist,
								   const PerambulatorInsertionInfo& insertinfo,
								   bool upper_spin_only
																	 )
{
// try{
    uPtr = new GaugeConfigurationInfo(gaugeinfo);
	gSmearPtr = new GluonSmearingInfo(gluonsmear);
    qSmearPtr = new QuarkSmearingInfo(quarksmear);
    qActionPtr = new QuarkActionInfo(quark);
	fPtr = new FileListInfo(flist);
	Nspin = (upper_spin_only) ? 2 : 4;

	insertPtr = new PerambulatorInsertionInfo(insertinfo);

	DHgetPtr=new DataGetHandlerMF<PerambulatorInsertionHandlerCL,FileKey,UIntKey,DataType>(
                    *this,flist,"Laph--QuarkPeramb","PerambulatorInsertionHandlerDataFile");
// }
// catch(...){
 //   cerr << "allocation problem in PerambulatorPerambulatorInsertionHandler"<<endl;
  //  exit(1);}
}


PerambulatorInsertionHandlerCL::~PerambulatorInsertionHandlerCL()
{
 clear();
}


void PerambulatorInsertionHandlerCL::clear()
{
 try{
    delete uPtr;
    delete gSmearPtr;
    delete qSmearPtr;
    delete qActionPtr;
    delete fPtr;

	delete insertPtr;
   }
 catch(...){ exit(1);}
 uPtr=0;
 gSmearPtr=0;
 qSmearPtr=0;
 qActionPtr=0;
 fPtr=0;

 insertPtr=0;

 delete DHgetPtr; DHgetPtr=0;
}


bool PerambulatorInsertionHandlerCL::isInfoSet() const
{
 return ((uPtr!=0)&&(gSmearPtr!=0)&&(qSmearPtr!=0)&&(qActionPtr!=0)
		 &&(fPtr!=0)&&(insertPtr!=0));
}

void PerambulatorInsertionHandlerCL::check_info_set(const string& name) const
{
 if (!isInfoSet()){
    cerr << "error in PerambulatorInsertionHandlerCL:"<<endl;
    cerr << "  must setInfo before calling "<<name<<endl;
    exit(1);}
}


const GaugeConfigurationInfo& PerambulatorInsertionHandlerCL::getGaugeConfigurationInfo() const 
{
 check_info_set("getGaugeConfigurationInfo");
 return *uPtr;
}

const GluonSmearingInfo& PerambulatorInsertionHandlerCL::getGluonSmearingInfo() const
{
 check_info_set("getGluonSmearingInfo");
 return *gSmearPtr;
}

const QuarkSmearingInfo& PerambulatorInsertionHandlerCL::getQuarkSmearingInfo() const
{
 check_info_set("getQuarkSmearingInfo");
 return *qSmearPtr;
}

const QuarkActionInfo& PerambulatorInsertionHandlerCL::getQuarkActionInfo() const 
{
 check_info_set("getQuarkActionInfo");
 return *qActionPtr;
}

const FileListInfo& PerambulatorInsertionHandlerCL::getFileListInfo() const 
{
 check_info_set("getFileListInfo");
 return *fPtr;
}

uint PerambulatorInsertionHandlerCL::getNumberOfLaplacianEigenvectors() const
{
 return qSmearPtr->getNumberOfLaplacianEigenvectors();
}

int PerambulatorInsertionHandlerCL::getTimeExtent() const 
{
 check_info_set("getTimeExtent");
 return uPtr->getTimeExtent();
}

void PerambulatorInsertionHandlerCL::getHeader(XMLHandler& xmlout) const
{
 if (isInfoSet()){
    xmlout.set_root("PerambulatorInsertionHandlerDataFile");
    XMLHandler gxml,gsxml,qsxml,qxml,qxmml; 
		uPtr->output(gxml);
    	gSmearPtr->output(gsxml);
    	qSmearPtr->output(qsxml);
    	qActionPtr->output(qxml);
		insertPtr->output(qxmml); 
		xmlout.put_child(gxml);
		xmlout.put_child(gsxml);
		xmlout.put_child(qsxml);
		xmlout.put_child(qxml);
		xmlout.put_child(qxmml);
    xmlout.put_child("NumSpinComponents",make_string(Nspin));
 }
}

bool PerambulatorInsertionHandlerCL::checkHeader(XMLHandler& xml_in, int suffix)
{
 if (xml_tag_count(xml_in,"PerambulatorInsertionHandlerDataFile")!=1) return false;
 XMLHandler xmlr(xml_in,"PerambulatorInsertionHandlerDataFile");
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
			 cerr <<"numspin = "<<numspin<<", Nspin = "<<Nspin<<endl;
       throw(std::invalid_argument("PerambulatorInsertion checkEqual failed...NumSpinComponents mismatch"));}
    qActionPtr->checkEqual(qaction_check);
    insertPtr->checkEqual(insert_check); }
 catch(const exception& xp){cout << xp.what()<<endl; return false;}
 return true;
}

void PerambulatorInsertionHandlerCL::writeHeader(XMLHandler& xmlout, 
                                     const PerambulatorInsertionHandlerCL::FileKey& fkey,
                                     int suffix)
{
 xmlout.set_root("PerambulatorHandlerDataFile");
 XMLHandler gxml,gsxml,qsxml,qxml,fxml;
 uPtr->output(gxml);
 gSmearPtr->output(gsxml);
 qSmearPtr->output(qsxml);
 qActionPtr->output(qxml);
 fkey.output(fxml);
 xmlout.put_child(gxml);
 xmlout.put_child(gsxml);
 xmlout.put_child(qsxml);
 xmlout.put_child(qxml);
 xmlout.put_child(fxml);
 xmlout.put_child("NumSpinComponents",make_string(Nspin));
 insertPtr->output(xmlout);
}


void PerambulatorInsertionHandlerCL::getSourceSinkTimes(set<array<int,2>>& source_sink_times) const
{
 source_sink_times.clear();
 check_info_set("getSourceSinkTimes");
 set<FileKey> fkeys=DHgetPtr->getFileKeys();
 for (const auto& fkey : fkeys) {
		 array<int,2> tarr={fkey.src_time,fkey.snk_time};
		 source_sink_times.insert(tarr);
 }
}


Array<cmplx> PerambulatorInsertionHandlerCL::getData(int snk_time, int snk_spin, 
                                              int src_time, int src_spin, int nEigsUse)
{
 uint nEigs=getNumberOfLaplacianEigenvectors();
 if ((src_spin<1)||(src_spin>Nspin)||(snk_spin<1)||(snk_spin>Nspin)
    ||(nEigsUse<1)||(nEigsUse>nEigs)){
    cerr << "Invalid parameters passed to getData in PerambulatorInsertionHandler"<<endl;
    exit(1);}
 Array<cmplx> result(nEigsUse,nEigsUse);
 FileKey fkey(snk_time,src_time);
 bool g5hermMode=false; 
 if (!(DHgetPtr->queryFile(fkey))) {
	 g5hermMode=true; 
	 fkey.src_time=snk_time; 
	 fkey.snk_time=src_time; 
	 if (!(DHgetPtr->queryFile(fkey)))
		 throw runtime_error("perambulator indices not found"); 
 }

 int snk_spin_tmp=snk_spin;
 int src_spin_tmp=src_spin;
 double fact=1.0;
 if (g5hermMode) {
	 //cout << "using G5 mode! "<<endl;
	 Array<vector<int>> inds(4,4); 
	 inds(0,0) = {2,2,-1};
	 inds(0,1) = {2,3,-1};
	 inds(0,2) = {2,0,1};
	 inds(0,3) = {2,1,1};
 
	 inds(1,0) = {3,2,-1};
	 inds(1,1) = {3,3,-1};
	 inds(1,2) = {3,0,1};
	 inds(1,3) = {3,1,1};
 
	 inds(2,0) = {0,2,1};
	 inds(2,1) = {0,3,1}; 
	 inds(2,2) = {0,0,-1};
	 inds(2,3) = {0,1,-1};

	 inds(3,0) = {1,2,1}; 
	 inds(3,1) = {1,3,1}; 
	 inds(3,2) = {1,0,-1}; 
	 inds(3,3) = {1,1,-1}; 
 	
	 src_spin_tmp = inds(snk_spin-1,src_spin-1)[0]+1;
	 snk_spin_tmp = inds(snk_spin-1,src_spin-1)[1]+1;
	 fact = inds(snk_spin-1,src_spin-1)[2];
 }

 vector<dcmplx> buff;
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    UIntKey rkey((src_spin_tmp-1)*nEigs+srcev_ind);
    buff=DHgetPtr->getData(fkey,rkey);
    for (int n=0;n<nEigsUse;n++){
			if (g5hermMode) {
				result(srcev_ind,n)=conj(fact*buff[(snk_spin_tmp-1)*nEigs+n]);
			} else 
				result(n,srcev_ind)=buff[(snk_spin-1)*nEigs+n];
		}}
 return result;
}


 // *****************************************************


bool PerambulatorInsertionHandlerCL::queryData(int src_time, int src_spin, int snk_time, 
                                    int snk_spin, int nEigsUse)
{
 uint nEigs=getNumberOfLaplacianEigenvectors();
 int minTime=0; 
 int maxTime=uPtr->getTimeExtent();
 if ((src_time<minTime)||(src_time>maxTime)||(snk_time<minTime)||(snk_time>maxTime)
    ||(src_spin<1)||(src_spin>Nspin)||(snk_spin<1)||(snk_spin>Nspin)){
    return false;}
 FileKey fkey(snk_time,src_time);
 for (int srcev_ind=0;srcev_ind<nEigsUse;srcev_ind++){
    UIntKey rkey((src_spin-1)*nEigs+srcev_ind);
    if (!DHgetPtr->queryData(fkey,rkey)) return false;}
 return true;
}

// ***************************************************************
}
 
