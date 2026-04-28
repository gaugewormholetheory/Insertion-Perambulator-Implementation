#ifndef __INLINE_CONVERT_PERAMBULATORS_INSERTION_H__
#define __INLINE_CONVERT_PERAMBULATORS_INSERTION_H__

#include "xml_handler.h"
#include "gauge_configuration_info.h"
#include "filelist_info.h"
#include "field_smearing_info.h"
#include "quark_action_info.h"
#include "perambulator_insertion.h"

// ************************************************************************
// *   convert perambulators to HDF5 format                               *
// *                                                                      *
// ************************************************************************


namespace LaphEnv {

// ************************************************

// void doPerambulatorConvert(XMLHandler& xmlh);
// void doPerambulatorConvertSG(XMLHandler& xmlh);
void doPerambulatorInsertionConvertCL(XMLHandler& xmlh);

// ***********************************************************
}
#endif
