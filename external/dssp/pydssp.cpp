#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mas.h"

#include "dssp.h"
#include "structure.h"

#include <boost/bind.hpp>
#include <boost/date_time/date_clock_device.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#if defined(_MSC_VER)
#include <conio.h>
#include <ctype.h>
#endif
#include <iostream>
#define foreach BOOST_FOREACH
namespace py=pybind11;
namespace dssp{

std::string get_dssp_from_coords(pybind11::array_t<double> coords_array,pybind11::bool_ reduece_ss=true){

    py::dict res;
    std::string dssp;
    // std::vector

    py::buffer_info buf= coords_array.request();
    assert(buf.shape[1]==4 && buf.ndim==3);
    size_t X=buf.shape[0],Y=buf.shape[1],Z=buf.shape[2];
    std::vector<std::vector<std::vector<double>>> coords;
    double *ptr = (double *) buf.ptr;
    for(size_t i=0;i<buf.shape[0];i++){
        coords.push_back(std::vector<std::vector<double>>());
        for(size_t j=0;j<buf.shape[1];j++){
            coords[i].push_back(std::vector<double>());
            for(size_t k=0;k<buf.shape[2];k++)
            coords[i][j].push_back(ptr[i*Y*Z+j*Z+k]);
        }
    }
    MProtein a;
    a.ReadCoords(coords);
    a.CalculateSecondaryStructure();
    std::vector<const MResidue*> residues;

    foreach (const MChain* chain, a.GetChains())
    {
        foreach (const MResidue* residue, chain->GetResidues())
        residues.push_back(residue);
    }

    // keep residues sorted by residue number as assigned during reading the PDB file
    sort(residues.begin(), residues.end(), boost::bind(&MResidue::GetNumber, _1) < boost::bind(&MResidue::GetNumber, _2));

    foreach (const MResidue* res, residues){
        char ss;
        const MResidue residue = *res;
        if(reduece_ss){
            switch (residue.GetSecondaryStructure())
            {
                case alphahelix:  ss = 'H'; break;
                case betabridge:  ss = 'E'; break;
                case strand:    ss = 'E'; break;
                case helix_3:    ss = 'H'; break;
                case helix_5:    ss = 'H'; break;
                case turn:      ss = 'L'; break;
                case bend:      ss = 'L'; break;
                case loop:      ss = 'L'; break;
            }
        }else{
            switch (residue.GetSecondaryStructure())
            {
                case alphahelix:  ss = 'H'; break;
                case betabridge:  ss = 'B'; break;
                case strand:    ss = 'E'; break;
                case helix_3:    ss = 'G'; break;
                case helix_5:    ss = 'I'; break;
                case turn:      ss = 'T'; break;
                case bend:      ss = 'S'; break;
                case loop:      ss = ' '; break;
            }
        }
        std::string NHO[2], ONH[2];
        int64 nNHO[2], nONH[2];
        const HBond* acceptors = residue.Acceptor();
        const HBond* donors = residue.Donor();
        for (uint32 i = 0; i < 2; ++i)
        {
            NHO[i] = ONH[i] = "0, 0.0";
            nNHO[i] = nONH[i] = 0;

            if (acceptors[i].residue != nullptr)
            {
            nNHO[i] = acceptors[i].residue->GetNumber() - residue.GetNumber();
            NHO[i] =acceptors[i].energy;
            }

            if (donors[i].residue != nullptr)
            {
            nONH[i] = donors[i].residue->GetNumber() - residue.GetNumber();
            ONH[i] = donors[i].energy;
            }
        }
        dssp.push_back(ss);
    }

    return dssp;
}



    PYBIND11_MODULE(dssp,m){
        m.def("get_dssp",&get_dssp_from_coords,"input coords of N CA C O");
    }

}


