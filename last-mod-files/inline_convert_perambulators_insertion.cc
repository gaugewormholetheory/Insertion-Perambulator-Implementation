#include "inline_convert_perambulators_insertion.h"
#include "perambulator_handler.h"
#include "perambulator_handler_sparse_grid.h"
#include "perambulator_insertion_handler_chroma_laph.h"
#include<array>
#include "H5Cpp.h"

using namespace std;
using namespace H5;

namespace LaphEnv {

// *********************************************************************

/*
void doPerambulatorConvert(XMLHandler& xmlr) 
{

    // read common information

 GaugeConfigurationInfo config_info(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 QuarkActionInfo qAction(xmlr);
 FileListInfo fList(xmlr);

 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>=1){
    upper_spin_only=true;}

    
 string outfile;
 xmlread(xmlr,"OutFile",outfile,"PerambulatorConvert");

 PerambulatorHandler pHand(config_info, gSmear, qSmear, qAction, fList,
		 upper_spin_only); 
 set<array<int,2> > source_sink_times;
 int Nspin=4; 
 if (upper_spin_only) 
	 Nspin=2;
 pHand.getSourceSinkTimes(source_sink_times);
 cout << "found "<<source_sink_times.size()<< " source/sink time combinations."<<endl;
 int nEv=qSmear.getNumberOfLaplacianEigenvectors();
 cout << "using "<<nEv<<" eigenvectors."<<endl;

 H5File data_file(outfile, H5F_ACC_TRUNC);
 Group simData(data_file.createGroup("/PerambulatorData"));
 DataSpace att_space(H5S_SCALAR); 

 XMLHandler hdr;
 pHand.getHeader(hdr);
 string hstr = hdr.str();
 cout << "header string: XXX"<<hstr<<"XXX"<<endl;

 Attribute xml_attr = simData.createAttribute("Header",StrType(0, H5T_VARIABLE), att_space);
 xml_attr.write(StrType(0, H5T_VARIABLE),&hstr);

 for (const auto& sst : source_sink_times) {
	 int t0 = sst[0];
	 int t = sst[1];
	 cout<<" doing (t0,t) = ("<<t0<<", "<<t<<")"<<endl;
	 stringstream sstr_t0_t;
	 sstr_t0_t<<"/PerambulatorData/srcTime"<<t0<<"_snkTime"<<t;
	 data_file.createGroup(sstr_t0_t.str());
	 for (int alpha=1; alpha<=Nspin ; alpha++) {
		 stringstream sstr_a;
		 sstr_a<<"/PerambulatorData/srcTime"<<t0<<"_snkTime"<<t
			 <<"/srcSpin"<<alpha;
		 data_file.createGroup(sstr_a.str());
		 for (int beta=1; beta<=Nspin ; beta++) {
			 //cout<<"(alpha,beta) = ("<<alpha<<", "<<beta<<"): "<<endl;
			 stringstream sstr_b;
			 sstr_b<<"/PerambulatorData/srcTime"<<t0<<"_snkTime"<<t<<
				 "/srcSpin"<<alpha<<"/snkSpin"<<beta;
			 data_file.createGroup(sstr_b.str());

			 //hsize_t dim[2];
			 //dim[0]=nEv;
			 //dim[1]=nEv;
			 hsize_t dim[1];
			 dim[0]=nEv*nEv;
			 DataSpace pspace(1, dim);

			 stringstream dname_re;
			 dname_re<<sstr_b.str()<<"/re"; 
			 DataSet pdat_re(data_file.createDataSet(dname_re.str(),PredType::NATIVE_DOUBLE, pspace));

			 stringstream dname_im;
			 dname_im<<sstr_b.str()<<"/im"; 
			 DataSet pdat_im(data_file.createDataSet(dname_im.str(),PredType::NATIVE_DOUBLE, pspace));

			 vector<double> p_re(nEv*nEv); 
			 vector<double> p_im(nEv*nEv);
			 Array<cmplx> p_c=pHand.getData(t,beta,t0,alpha,nEv);
			 for (int n=0;n<nEv;n++)
				 for (int m=0;m<nEv;m++) {
					 p_re[n*nEv+m]=real(p_c(n,m));
					 p_im[n*nEv+m]=imag(p_c(n,m));
					 //cout <<"ind = "<<n*nEv+m<<",    ("<<n<<", "<<m<<") = "<<"("<<p_re[n*nEv+m]<<", "<<
						 //p_im[n*nEv+m]<<")"<<endl;
				 }
			 pdat_re.write(&p_re[0],PredType::NATIVE_DOUBLE,pspace);
			 pdat_im.write(&p_im[0],PredType::NATIVE_DOUBLE,pspace);
		 }
	 }
	}
 }
*/
//-----------------------------------------
/*
void doPerambulatorConvertSG(XMLHandler& xmlr) 
{

    // read common information

 GaugeConfigurationInfo config_info(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 QuarkActionInfo qAction(xmlr);
 FileListInfo fList(xmlr);
 RandomSparseGrid sGrid(xmlr);

 cout << sGrid.output() << endl;

 const int nColor = 3; 
 const int nGrid = sGrid.getNGridPoints(); 
 int nTot = nColor*nGrid;
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>=1){
    upper_spin_only=true;}

    
 string outfile;
 xmlread(xmlr,"OutFile",outfile,"PerambulatorConvertSparseGrid");

 PerambulatorHandlerSparseGrid pHand(config_info, gSmear, qSmear, qAction, 
		 fList, sGrid, upper_spin_only); 
 set<array<int,2> > source_sink_times;
 int Nspin=4; 
 if (upper_spin_only) 
	 Nspin=2;
 pHand.getSourceSinkTimes(source_sink_times);
 cout << "found "<<source_sink_times.size()<< " source/sink time combinations."<<endl;
 int nEv=qSmear.getNumberOfLaplacianEigenvectors();
 cout << "using "<<nEv<<" eigenvectors."<<endl;

 H5File data_file(outfile, H5F_ACC_TRUNC);
 Group simData(data_file.createGroup("/PerambulatorSparseGridData"));
 DataSpace att_space(H5S_SCALAR); 

 XMLHandler hdr;
 pHand.getHeader(hdr);
 string hstr = hdr.str();
 cout << "header string: XXX"<<hstr<<"XXX"<<endl;

 Attribute xml_attr = simData.createAttribute("Header",StrType(0, H5T_VARIABLE), att_space);
 xml_attr.write(StrType(0, H5T_VARIABLE),&hstr);

 for (const auto& sst : source_sink_times) {
	 int t0 = sst[0];
	 int t = sst[1];
	 cout<<" doing (t0,t) = ("<<t0<<", "<<t<<")"<<endl;
	 stringstream sstr_t0_t;
	 sstr_t0_t<<"/PerambulatorSparseGridData/srcTime"<<t0<<"_snkTime"<<t;
	 data_file.createGroup(sstr_t0_t.str());
	 for (int alpha=1; alpha<=Nspin ; alpha++) {
		 stringstream sstr_a;
		 sstr_a<<"/PerambulatorSparseGridData/srcTime"<<t0<<"_snkTime"<<t
			 <<"/srcSpin"<<alpha;
		 data_file.createGroup(sstr_a.str());
		 for (int beta=1; beta<=Nspin ; beta++) {
			 //cout<<"(alpha,beta) = ("<<alpha<<", "<<beta<<"): "<<endl;
			 stringstream sstr_b;
			 sstr_b<<"/PerambulatorSparseGridData/srcTime"<<t0<<"_snkTime"<<t<<
				 "/srcSpin"<<alpha<<"/snkSpin"<<beta;
			 data_file.createGroup(sstr_b.str());

			 //hsize_t dim[2];
			 //dim[0]=nEv;
			 //dim[1]=nEv;
			 hsize_t dim[1];
			 dim[0]=nTot*nEv;
			 DataSpace pspace(1, dim);

			 stringstream dname_re;
			 dname_re<<sstr_b.str()<<"/re"; 
			 DataSet pdat_re(data_file.createDataSet(dname_re.str(),PredType::NATIVE_DOUBLE, pspace));

			 stringstream dname_im;
			 dname_im<<sstr_b.str()<<"/im"; 
			 DataSet pdat_im(data_file.createDataSet(dname_im.str(),PredType::NATIVE_DOUBLE, pspace));

			 vector<double> p_re(nTot*nEv); 
			 vector<double> p_im(nTot*nEv);
			 Array<cmplx> p_c=pHand.getData(t,beta,t0,alpha,nEv);
			 for (int n=0;n<nTot;n++)
				 for (int m=0;m<nEv;m++) {
					 p_re[n*nEv+m]=real(p_c(n,m));
					 p_im[n*nEv+m]=imag(p_c(n,m));
					 //cout <<"ind = "<<n*nEv+m<<",    ("<<n<<", "<<m<<") = "<<"("<<p_re[n*nEv+m]<<", "<<
						 //p_im[n*nEv+m]<<")"<<endl;
				 }
			 pdat_re.write(&p_re[0],PredType::NATIVE_DOUBLE,pspace);
			 pdat_im.write(&p_im[0],PredType::NATIVE_DOUBLE,pspace);
		 }
	 }
	}
 }
*/
//-----------------------------------------


void doPerambulatorInsertionConvertCL(XMLHandler& xmlr) 
{

    // read common information

 GaugeConfigurationInfo config_info(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 QuarkActionInfo qAction(xmlr);
 FileListInfo fList(xmlr);

 PerambulatorInsertionInfo insertinfo(xmlr);

 cout << "read config info"<<endl;
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>=1){
    upper_spin_only=true;}

    
 string outfile;
 xmlread(xmlr,"OutFile",outfile,"PerambulatorConvertCL");

 cout << "creating pHand"<<endl;
 PerambulatorInsertionHandlerCL pHand(config_info, gSmear, qSmear, qAction, fList, insertinfo,
		 upper_spin_only); 
 cout << "pHand created"<<endl;
 set<array<int,2> > source_sink_times;
 int Nspin=4; 
 if (upper_spin_only) 
	 Nspin=2;
 pHand.getSourceSinkTimes(source_sink_times);
 cout << "found "<<source_sink_times.size()<< " source/sink time combinations."<<endl;
 int nEv=qSmear.getNumberOfLaplacianEigenvectors();
 cout << "using "<<nEv<<" eigenvectors."<<endl;

 H5File data_file(outfile, H5F_ACC_TRUNC);
 Group simData(data_file.createGroup("/PerambulatorData"));
 DataSpace att_space(H5S_SCALAR); 

 XMLHandler hdr;
 pHand.getHeader(hdr);
 string hstr = hdr.str();
 cout << "header string: XXX"<<hstr<<"XXX"<<endl;

 Attribute xml_attr = simData.createAttribute("Header",StrType(0, H5T_VARIABLE), att_space);
 xml_attr.write(StrType(0, H5T_VARIABLE),&hstr);

 for (const auto& sst : source_sink_times) {
	 int t0 = sst[0];
	 int t = sst[1];
	 cout<<" doing (t0,t) = ("<<t0<<", "<<t<<")"<<endl;
	 stringstream sstr_t0_t;
	 sstr_t0_t<<"/PerambulatorData/srcTime"<<t0<<"_snkTime"<<t;
	 data_file.createGroup(sstr_t0_t.str());
	 for (int alpha=1; alpha<=Nspin ; alpha++) {
		 stringstream sstr_a;
		 sstr_a<<"/PerambulatorData/srcTime"<<t0<<"_snkTime"<<t
			 <<"/srcSpin"<<alpha;
		 data_file.createGroup(sstr_a.str());
		 for (int beta=1; beta<=Nspin ; beta++) {
			 //cout<<"(alpha,beta) = ("<<alpha<<", "<<beta<<"): "<<endl;
			 stringstream sstr_b;
			 sstr_b<<"/PerambulatorData/srcTime"<<t0<<"_snkTime"<<t<<
				 "/srcSpin"<<alpha<<"/snkSpin"<<beta;
			 data_file.createGroup(sstr_b.str());

			 //hsize_t dim[2];
			 //dim[0]=nEv;
			 //dim[1]=nEv;
			 hsize_t dim[1];
			 dim[0]=nEv*nEv;
			 DataSpace pspace(1, dim);

			 stringstream dname_re;
			 dname_re<<sstr_b.str()<<"/re"; 
			 DataSet pdat_re(data_file.createDataSet(dname_re.str(),PredType::NATIVE_DOUBLE, pspace));

			 stringstream dname_im;
			 dname_im<<sstr_b.str()<<"/im"; 
			 DataSet pdat_im(data_file.createDataSet(dname_im.str(),PredType::NATIVE_DOUBLE, pspace));

			 vector<double> p_re(nEv*nEv); 
			 vector<double> p_im(nEv*nEv);
			 Array<cmplx> p_c=pHand.getData(t,beta,t0,alpha,nEv);
			 for (int n=0;n<nEv;n++)
				 for (int m=0;m<nEv;m++) {
					 p_re[n*nEv+m]=real(p_c(n,m));
					 p_im[n*nEv+m]=imag(p_c(n,m));
					 //cout <<"ind = "<<n*nEv+m<<",    ("<<n<<", "<<m<<") = "<<"("<<p_re[n*nEv+m]<<", "<<
						 //p_im[n*nEv+m]<<")"<<endl;
				 }
			 pdat_re.write(&p_re[0],PredType::NATIVE_DOUBLE,pspace);
			 pdat_im.write(&p_im[0],PredType::NATIVE_DOUBLE,pspace);
		 }
	 }
	}
 }
}

