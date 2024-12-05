#include "utils/ligand_util.hh"
#include "utils/random_util.hh"

#include <cmath>
#include <iostream>

#include "utils/utils.hh"

namespace utils {

using namespace basic;


EigenXform random_rotate_pNPA(std::string & xform_str) {
    // perturb ligand
    // shift down
    EigenXform Z_shift_down(EigenXform::Identity());
    Z_shift_down.translation() = Vec(0.0, 0.0, -1.80);

    // shift up
    EigenXform Z_shift_up(EigenXform::Identity());
    Z_shift_up.translation() = Vec(0.0, 0.0, 1.80);
    
    // rotate xoy
    Real x = utils::random_real(-1,1);
    Real y = std::sqrt(1 - x*x);
    y *= (utils::random_real(-1,1) > 0?1:-1);
    Real ang1 = utils::random_real(-25,25);
    EigenXform xform_rotate_xoy(EigenXform::Identity());
    xform_rotate_xoy.rotate( AngleAxis( ang1 /180.0*3.14159265359, Vec(x,y,0.0) ) );

    // rotate Z
    Real ang2 = utils::random_real(0,360);
    EigenXform xform_rotate_Z(EigenXform::Identity());
    xform_rotate_Z.rotate( AngleAxis( ang2 /180.0*3.14159265359, Vec(0.0,0.0,1.0) ) );

    EigenXform xform = Z_shift_up * xform_rotate_Z * xform_rotate_xoy * Z_shift_down;

    std::stringstream xform_stream;
    xform_stream << "REMARK PDBinfo-LABEL:    1 LIGAND_XFORM:";
    xform_stream << xform(0,0) << "_" << xform(0,1) << "_" << xform(0,2) << "_" << xform(0,3) << "_"
              << xform(1,0) << "_" << xform(1,1) << "_" << xform(1,2) << "_" << xform(1,3) << "_"
              << xform(2,0) << "_" << xform(2,1) << "_" << xform(2,2) << "_" << xform(2,3);

    xform_str = xform_stream.str();


    return xform;
}

Real pNPA_neighbors(scene::Pose const & pose, EigenXform const & ligand_xform) {

    // prepare the pNPA coords
    // 'O2','O0','O4','C2','C3','C4','C6','C7'

    /*
    HETATM   12  ZN  ZNX A  66       0.000   0.000   0.000  1.00  0.00           X  
    HETATM    1  C0  NPA X   1      -2.536   3.647   4.652  1.00 20.00           C
    HETATM    2  C1  NPA X   1      -1.141   3.638   4.697  1.00 20.00           C
    HETATM    3  C2  NPA X   1      -0.469   2.437   4.523  1.00 20.00           C
    HETATM    4  C3  NPA X   1      -1.184   1.249   4.293  1.00 20.00           C
    HETATM    5  C4  NPA X   1      -2.588   1.279   4.231  1.00 20.00           C
    HETATM    6  C5  NPA X   1      -3.264   2.481   4.419  1.00 20.00           C
    HETATM    7  O0  NPA X   1      -0.423   0.131   4.097  1.00 20.00           O
    HETATM    8  C6  NPA X   1      -0.573  -0.674   2.896  1.00 20.00           C
    HETATM    9  C7  NPA X   1       0.121  -1.987   3.222  1.00 20.00           C
    HETATM   10  N0  NPA X   1      -3.252   4.910   4.842  1.00 20.00           N
    HETATM   11  O1  NPA X   1      -2.580   5.927   5.030  1.00 20.00           O
    HETATM   12  O2  NPA X   1       0.000   0.000   1.800  1.00 20.00           O
    HETATM   13  O3  NPA X   1      -4.482   4.885   4.804  1.00 20.00           O
    HETATM   14  O4  NPA X   1      -1.891  -0.923   2.588  1.00 20.00           O
    HETATM   15  H0  NPA X   1      -0.617   4.512   4.857  1.00 20.00           H
    HETATM   16  H1  NPA X   1       0.562   2.416   4.563  1.00 20.00           H
    HETATM   17  H2  NPA X   1      -3.116   0.412   4.046  1.00 20.00           H
    HETATM   18  H3  NPA X   1      -4.295   2.509   4.386  1.00 20.00           H
    HETATM   19  H4  NPA X   1       0.055  -2.644   2.380  1.00 20.00           H
    HETATM   20  H5  NPA X   1      -0.354  -2.441   4.067  1.00 20.00           H
    HETATM   21  H6  NPA X   1       1.150  -1.800   3.449  1.00 20.00           H
    HETATM   22  H7  NPA X   1       0.901   0.255   1.826  1.00 20.00           H
    HETATM   23  H8  NPA X   1      -2.440  -1.369   3.203  1.00 20.00           H
    */

    Eigen::Matrix<Real, 8, 3> pNPA_coords;
    pNPA_coords <<  0.000,   0.000,   1.800,
                   -0.423,   0.131,   4.097,
                   -1.891,  -0.923,   2.588,
                   -0.469,   2.437,   4.523,
                   -1.184,   1.249,   4.293,
                   -2.588,   1.279,   4.231,
                   -0.573,  -0.674,   2.896,
                    0.121,  -1.987,   3.222;

    pNPA_coords = ( ( ligand_xform.rotation() * ( pNPA_coords.transpose() ) ).colwise() + ligand_xform.translation() ).transpose();


    return ligand_neighbors(pose, pNPA_coords);
}

EigenXform random_place_AMA(std::string & xform_str)
{

    EigenXform xform_flip(EigenXform::Identity());   
    xform_flip.rotate( AngleAxis( 3.14159265359, Vec(1.0,0.0,0.0) ) );
    Real shift = utils::random_real(-13,-7);
    EigenXform xform_translate(EigenXform::Identity());
    xform_translate.translation() = Vec(0.0, 0.0, shift);
    Real rotate = utils::random_real(0,360);
    EigenXform xform_rotate(EigenXform::Identity());
    xform_rotate.rotate( AngleAxis( rotate /180.0*3.14159265359, Vec(0.0,0.0,1.0) ) );
    
    EigenXform xform = xform_rotate * xform_translate * xform_flip;


    std::stringstream xform_stream;
    xform_stream << "REMARK PDBinfo-LABEL:    1 LIGAND_XFORM:";
    xform_stream << xform(0,0) << "_" << xform(0,1) << "_" << xform(0,2) << "_" << xform(0,3) << "_"
              << xform(1,0) << "_" << xform(1,1) << "_" << xform(1,2) << "_" << xform(1,3) << "_"
              << xform(2,0) << "_" << xform(2,1) << "_" << xform(2,2) << "_" << xform(2,3);

    xform_str = xform_stream.str();

    return xform;
}

Real AMA_neighbors(scene::Pose const & pose, EigenXform const & ligand_xform)
{
    /*
    HETATM    1  C1  AMA X   1       0.033   1.484  -0.788  1.00 20.00           C
    HETATM    2  N1  AMA X   1       0.000   0.000   2.759  1.00 20.00           N
    HETATM    3  C2  AMA X   1      -1.256   0.775  -1.310  1.00 20.00           C
    HETATM    4  C3  AMA X   1      -1.296  -0.712  -0.810  1.00 20.00           C
    HETATM    5  C4  AMA X   1      -0.028  -1.475  -1.311  1.00 20.00           C
    HETATM    6  C5  AMA X   1       1.271  -0.766  -0.788  1.00 20.00           C
    HETATM    7  C6  AMA X   1       1.303   0.712  -1.294  1.00 20.00           C
    HETATM    8  C7  AMA X   1       0.000   1.462   0.760  1.00 20.00           C
    HETATM    9  C8  AMA X   1       1.263  -0.739   0.760  1.00 20.00           C
    HETATM   10  C9  AMA X   1      -1.279  -0.735   0.739  1.00 20.00           C
    HETATM   11  C10 AMA X   1      -0.010  -0.007   1.282  1.00 20.00           C
    HETATM   12  H1  AMA X   1       0.738   0.405   3.049  1.00 20.00           H
    HETATM   13  H2  AMA X   1      -0.016  -0.837   3.059  1.00 20.00           H
    HETATM   14  H3  AMA X   1       0.056   2.401  -1.103  1.00 20.00           H
    HETATM   15  H4  AMA X   1      -1.261   0.788  -2.280  1.00 20.00           H
    HETATM   16  H5  AMA X   1      -2.097  -1.152  -1.135  1.00 20.00           H
    HETATM   17  H6  AMA X   1      -0.051  -2.388  -0.981  1.00 20.00           H
    HETATM   18  H7  AMA X   1       2.059  -1.238  -1.103  1.00 20.00           H
    HETATM   19  H8  AMA X   1       2.101   1.150  -0.957  1.00 20.00           H
    HETATM   20  H9  AMA X   1       0.784   1.919   1.104  1.00 20.00           H
    HETATM   21  H10 AMA X   1       1.268  -1.650   1.095  1.00 20.00           H
    HETATM   22  H11 AMA X   1      -2.071  -0.288   1.073  1.00 20.00           H
    HETATM   23  H12 AMA X   1      -0.716   0.437   3.059  1.00 20.00           H
    HETATM   24  H13 AMA X   1      -2.037   1.246  -0.979  1.00 20.00           H
    HETATM   25  H14 AMA X   1      -0.021  -1.484  -2.280  1.00 20.00           H
    HETATM   26  H15 AMA X   1       1.318   0.719  -2.264  1.00 20.00           H
    HETATM   27  H16 AMA X   1      -0.804   1.911   1.066  1.00 20.00           H
    HETATM   28  H17 AMA X   1       2.050  -0.268   1.074  1.00 20.00           H
    HETATM   29  H18 AMA X   1      -1.269  -1.657   1.042  1.00 20.00           H
    */


    Eigen::Matrix<Real, 11, 3> AMA_coords;
    AMA_coords <<   0.033,   1.484,  -0.788,
                    0.000,   0.000,   2.759,
                   -1.256,   0.775,  -1.310,
                   -1.296,  -0.712,  -0.810,
                   -0.028,  -1.475,  -1.311,
                    1.271,  -0.766,  -0.788,
                    1.303,   0.712,  -1.294,
                    0.000,   1.462,   0.760,
                    1.263,  -0.739,   0.760,
                   -1.279,  -0.735,   0.739,
                   -0.010,  -0.007,   1.282;

    AMA_coords = ( ( ligand_xform.rotation() * ( AMA_coords.transpose() ) ).colwise() + ligand_xform.translation() ).transpose();


    return ligand_neighbors(pose, AMA_coords);
}

Real HQA_neighbors(scene::Pose const & pose, EigenXform const & ligand_xform)
{
    /*
    HETATM    1  N   UNL 1   1      -0.973  -1.754   0.003  1.00  0.00           N  
    HETATM    2  C   UNL 1   1       1.834  -4.781   0.004  1.00  0.00           C  
    HETATM    3  C   UNL 1   1       0.498  -5.156   0.007  1.00  0.00           C  
    HETATM    4  C   UNL 1   1      -0.498  -4.146   0.007  1.00  0.00           C  
    HETATM    5  C   UNL 1   1      -0.080  -2.782   0.004  1.00  0.00           C  
    HETATM    6  C   UNL 1   1       1.310  -2.391   0.001  1.00  0.00           C  
    HETATM    7  C   UNL 1   1       2.250  -3.431   0.001  1.00  0.00           C  
    HETATM    8  C   UNL 1   1      -2.276  -1.982   0.006  1.00  0.00           C  
    HETATM    9  C   UNL 1   1      -2.783  -3.302   0.009  1.00  0.00           C  
    HETATM   10  C   UNL 1   1      -1.902  -4.369   0.009  1.00  0.00           C  
    HETATM   11  O   UNL 1   1       1.596  -1.117  -0.002  1.00  0.00           O  
    HETATM   18 CU   UNL 1   1       0.000   0.000   0.000  1.00  0.00          Cu  
    HETATM   19  N   UNL 1   1       0.973   1.754   0.003  1.00  0.00           N  
    HETATM   20  C   UNL 1   1      -1.834   4.781   0.004  1.00  0.00           C  
    HETATM   21  C   UNL 1   1      -0.498   5.156   0.007  1.00  0.00           C  
    HETATM   22  C   UNL 1   1       0.498   4.146   0.007  1.00  0.00           C  
    HETATM   23  C   UNL 1   1       0.080   2.782   0.004  1.00  0.00           C  
    HETATM   24  C   UNL 1   1      -1.310   2.391   0.001  1.00  0.00           C  
    HETATM   25  C   UNL 1   1      -2.250   3.431   0.001  1.00  0.00           C  
    HETATM   26  C   UNL 1   1       2.276   1.982   0.006  1.00  0.00           C  
    HETATM   27  C   UNL 1   1       2.783   3.302   0.009  1.00  0.00           C  
    HETATM   28  C   UNL 1   1       1.902   4.369   0.009  1.00  0.00           C  
    HETATM   29  O   UNL 1   1      -1.596   1.117  -0.002  1.00  0.00           O  
    */


    Eigen::Matrix<Real, 23, 3> HQA_coords;
    HQA_coords <<  -0.973,  -1.754,   0.003,
                    1.834,  -4.781,   0.004,
                    0.498,  -5.156,   0.007,
                   -0.498,  -4.146,   0.007,
                   -0.080,  -2.782,   0.004,
                    1.310,  -2.391,   0.001,
                    2.250,  -3.431,   0.001,
                   -2.276,  -1.982,   0.006,
                   -2.783,  -3.302,   0.009,
                   -1.902,  -4.369,   0.009,
                    1.596,  -1.117,  -0.002,
                    0.000,   0.000,   0.000,
                    0.973,   1.754,   0.003,
                   -1.834,   4.781,   0.004,
                   -0.498,   5.156,   0.007,
                    0.498,   4.146,   0.007,
                    0.080,   2.782,   0.004,
                   -1.310,   2.391,   0.001,
                   -2.250,   3.431,   0.001,
                    2.276,   1.982,   0.006,
                    2.783,   3.302,   0.009,
                    1.902,   4.369,   0.009,
                   -1.596,   1.117,  -0.002;

    HQA_coords = ( ( ligand_xform.rotation() * ( HQA_coords.transpose() ) ).colwise() + ligand_xform.translation() ).transpose();


    return ligand_neighbors(pose, HQA_coords);
}


EigenXform random_place_DTG(std::string & xform_str)
{
    // ?????????????????????????????????????????????????????
    // ?????????????????????????????????????????????????????
    // ?????????????????????????????????????????????????????
    // ??????????????    Not Done Yet    ???????????????????
    // ?????????????????????????????????????????????????????
    // ?????????????????????????????????????????????????????
    // ?????????????????????????????????????????????????????
    EigenXform xform_flip(EigenXform::Identity());   
    xform_flip.rotate( AngleAxis( 3.14159265359, Vec(1.0,0.0,0.0) ) );
    Real shift = utils::random_real(-13,-7);
    EigenXform xform_translate(EigenXform::Identity());
    xform_translate.translation() = Vec(0.0, 0.0, shift);
    Real rotate = utils::random_real(0,360);
    EigenXform xform_rotate(EigenXform::Identity());
    xform_rotate.rotate( AngleAxis( rotate /180.0*3.14159265359, Vec(0.0,0.0,1.0) ) );
    
    EigenXform xform = xform_rotate * xform_translate * xform_flip;


    std::stringstream xform_stream;
    xform_stream << "REMARK PDBinfo-LABEL:    1 LIGAND_XFORM:";
    xform_stream << xform(0,0) << "_" << xform(0,1) << "_" << xform(0,2) << "_" << xform(0,3) << "_"
              << xform(1,0) << "_" << xform(1,1) << "_" << xform(1,2) << "_" << xform(1,3) << "_"
              << xform(2,0) << "_" << xform(2,1) << "_" << xform(2,2) << "_" << xform(2,3);

    xform_str = xform_stream.str();

    return xform;

}


}
