#ifndef INCLUDED_basic_macros_hh
#define INCLUDED_basic_macros_hh

namespace basic {

#define idl_R_N 1.32869
#define idl_R_CA  1.458
#define idl_R_C  1.52326

#define idl_sinW_N  0.8508111094240511
#define idl_sinW_CA  0.9323238012155122
#define idl_sinW_C  0.8972583696743285

#define idl_cosW_N  -0.5254716510722679
#define idl_cosW_CA  -0.3616245700820923
#define idl_cosW_C  -0.44150585279174515

#define M_PI 3.14159265358979323846

// numeric stability

#define PSEUDO_CB_X -0.524098
#define PSEUDO_CB_Y -0.772137
#define PSEUDO_CB_Z -1.206615

#define VALINE_CB_X  -0.510472
#define VALINE_CB_Y  -0.771189
#define VALINE_CB_Z  -1.231671

#define IDL_H_X -0.491967
#define IDL_H_Y -0.882065
#define IDL_H_Z 0.0

#define IDL_O_X -0.670434
#define IDL_O_Y -1.032435
#define IDL_O_Z 0.0


#define CALC_CENTER_POSITION_PREV(WIDTH, STR) (((WIDTH + ((int)strlen(STR))) % 2) \
       ? ((WIDTH + ((int)strlen(STR)) - 1)/2) : ((WIDTH + ((int)strlen(STR)))/2))
#define CALC_CENTER_POSITION_POST(WIDTH, STR) (((WIDTH - ((int)strlen(STR))) % 2) \
       ? ((WIDTH - ((int)strlen(STR)) + 1)/2) : ((WIDTH - ((int)strlen(STR)))/2))

}

#endif
