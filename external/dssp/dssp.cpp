// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

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

#if defined(_MSC_VER)
#include <conio.h>
#include <ctype.h>
#endif
#include <iostream>


#define foreach BOOST_FOREACH


std::string FixedLengthString(const int64 number, const uint64 length)
{
    std::string s = std::to_string(number);

    if (s.length() < length)
        s = s.insert(0, length - s.length(), ' ');
    else if (s.length() > length)
    {
        // make and arrow: "--->"
        s = std::string(length - 1, '-') + ">";
    }
    return s;
}


std::string ResidueToDSSPLine(const MResidue& residue)
{
/*
  This is the header line for the residue lines in a DSSP file:

  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA           CHAIN AUTHCHAIN
 */
  boost::format kDSSPResidueLine(
  "%5.5s%5.5s%1.1s%1.1s %c  %c %c%c%c%c%c%c%c%4.4s%4.4s%c%4.4s %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f             %4.4s      %4.4s %10s %10s %10s %10s %10s %10s %10s %10s");

  const MAtom& ca = residue.GetCAlpha();

  char code = kResidueInfo[residue.GetType()].code;
  if (residue.GetType() == kCysteine and residue.GetSSBridgeNr() != 0)
    code = 'a' + ((residue.GetSSBridgeNr() - 1) % 26);

  char ss;
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

  char helix[3];
  for (uint32 stride = 3; stride <= 5; ++stride)
  {
    switch (residue.GetHelixFlag(stride))
    {
      case helixNone:      helix[stride - 3] = ' '; break;
      case helixStart:    helix[stride - 3] = '>'; break;
      case helixEnd:      helix[stride - 3] = '<'; break;
      case helixStartAndEnd:  helix[stride - 3] = 'X'; break;
      case helixMiddle:    helix[stride - 3] = '0' + stride; break;
    }
  }

  char bend = ' ';
  if (residue.IsBend())
    bend = 'S';

  double alpha;
  char chirality;
  std::tie(alpha,chirality) = residue.Alpha();

  int64 bp[2] = {};
  char bridgelabel[2] = { ' ', ' ' };
  for (uint32 i = 0; i < 2; ++i)
  {
    MBridgeParner p = residue.GetBetaPartner(i);
    if (p.residue != nullptr)
    {
      bp[i] = p.residue->GetNumber();
      bridgelabel[i] = 'A' + p.ladder % 26;
      if (p.parallel)
        bridgelabel[i] = tolower(bridgelabel[i]);
    }
  }

  char sheet = ' ';
  if (residue.GetSheet() != 0)
    sheet = 'A' + (residue.GetSheet() - 1) % 26;

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
      NHO[i] = (boost::format("%s,%3.1f") % FixedLengthString(nNHO[i], 5) % acceptors[i].energy).str();
    }

    if (donors[i].residue != nullptr)
    {
      nONH[i] = donors[i].residue->GetNumber() - residue.GetNumber();
      ONH[i] = (boost::format("%s,%3.1f") % FixedLengthString(nONH[i], 5) % donors[i].energy).str();
    }
  }

  std::string chainChar = ca.mChainID,
                          long_ChainID1 = ca.mChainID,
                          long_ChainID2 = ca.mAuthChainID;
  if (ca.mChainID.length () > 1)
  {
    // For mmCIF compatibility

    chainChar = ">";
  }

  return (kDSSPResidueLine % FixedLengthString(residue.GetNumber(), 5) % FixedLengthString(ca.mResSeq, 5) %
    ca.mICode % chainChar % code %
    ss % helix[0] % helix[1] % helix[2] % bend % chirality % bridgelabel[0] % bridgelabel[1] %
    FixedLengthString(bp[0], 4) % FixedLengthString(bp[1], 4) % sheet % FixedLengthString(floor(residue.Accessibility() + 0.5), 4) %
    NHO[0] % ONH[0] % NHO[1] % ONH[1] %
    residue.TCO() % residue.Kappa() % alpha % residue.Phi() % residue.Psi() %
    ca.mLoc.mX % ca.mLoc.mY % ca.mLoc.mZ % long_ChainID1 % long_ChainID2 %
    FixedLengthString(residue.GetNumber(), 10) % FixedLengthString(ca.mResSeq, 10) %
    FixedLengthString(bp[0], 10) % FixedLengthString(bp[1], 10) %
    FixedLengthString(nNHO[0], 10) % FixedLengthString(nONH[0], 10) % FixedLengthString(nNHO[1], 10) % FixedLengthString(nONH[1], 10)
  ).str();
}

void WriteDSSP(MProtein& protein, std::ostream& os)
{
  const std::string kFirstLine("==== Secondary Structure Definition by the program DSSP, CMBI version "  "                          ==== ");
  boost::format kHeaderLine("%1% %|127t|%2%");

  using namespace boost::gregorian;

  uint32 nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, nrOfHBonds;
  uint32 nrOfHBondsPerDistance[11] = {};

  protein.GetStatistics(nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, nrOfHBonds, nrOfHBondsPerDistance);

  date today = day_clock::local_day();

  os << kHeaderLine % (kFirstLine + "DATE=" + to_iso_extended_string(today)) % '.' << std::endl;
  os << kHeaderLine % "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637" % '.' << std::endl;
  os << kHeaderLine % protein.GetHeader() % '.' << std::endl;
  if (not protein.GetCompound().empty())
    os << kHeaderLine % protein.GetCompound() % '.' << std::endl;
  if (not protein.GetSource().empty())
    os << kHeaderLine % protein.GetSource() % '.' << std::endl;
  if (not protein.GetAuthor().empty())
    os << kHeaderLine % protein.GetAuthor() % '.' << std::endl;

  double accessibleSurface = 0;  // calculate accessibility as
  foreach (const MChain* chain, protein.GetChains())
  {
    foreach (const MResidue* residue, chain->GetResidues())
      accessibleSurface += residue->Accessibility();
  }

  os << boost::format("%5.5d%3.3d%3.3d%3.3d%3.3d TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN) %|127t|%c") %
       nrOfResidues % nrOfChains % nrOfSSBridges % nrOfIntraChainSSBridges % (nrOfSSBridges - nrOfIntraChainSSBridges) % '.' << std::endl;
     os << kHeaderLine % (boost::format("%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)") % accessibleSurface) % '.' << std::endl;

  // hydrogenbond summary

  os << kHeaderLine % (
    boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES")
      % nrOfHBonds % (nrOfHBonds * 100.0 / nrOfResidues)) % '.' << std::endl;

  uint32 nrOfHBondsInParallelBridges = protein.GetNrOfHBondsInParallelBridges();
  os << kHeaderLine % (
    boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES")
      % nrOfHBondsInParallelBridges % (nrOfHBondsInParallelBridges * 100.0 / nrOfResidues)) % '.' << std::endl;

  uint32 nrOfHBondsInAntiparallelBridges = protein.GetNrOfHBondsInAntiparallelBridges();
  os << kHeaderLine % (
    boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES")
      % nrOfHBondsInAntiparallelBridges % (nrOfHBondsInAntiparallelBridges * 100.0 / nrOfResidues)) % '.' << std::endl;

  boost::format kHBondsLine("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%1.1d), SAME NUMBER PER 100 RESIDUES");
     for (int32 k = 0; k < 11; ++k)
  {
    os << kHeaderLine % (kHBondsLine % nrOfHBondsPerDistance[k] % (nrOfHBondsPerDistance[k] * 100.0 / nrOfResidues) % (k - 5 < 0 ? '-' : '+') % abs(k - 5)) % '.' << std::endl;
  }

  // histograms...

  uint32 histogram[kHistogramSize];
  os << "  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           ." << std::endl;

  protein.GetResiduesPerAlphaHelixHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << boost::format("%3.3d") % histogram[i];
  os << "    RESIDUES PER ALPHA HELIX         ." << std::endl;

  protein.GetParallelBridgesPerLadderHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << boost::format("%3.3d") % histogram[i];
  os << "    PARALLEL BRIDGES PER LADDER      ." << std::endl;

  protein.GetAntiparallelBridgesPerLadderHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << boost::format("%3.3d") % histogram[i];
  os << "    ANTIPARALLEL BRIDGES PER LADDER  ." << std::endl;

  protein.GetLaddersPerSheetHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << boost::format("%3.3d") % histogram[i];
  os << "    LADDERS PER SHEET                ." << std::endl;

  // per residue information

  os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN AUTHCHAIN     NUMBER     RESNUM        BP1        BP2    N-H-->O    O-->H-N    N-H-->O    O-->H-N" << std::endl;
  boost::format kDSSPResidueLine(
    "%5.5s        !%c             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");

  std::vector<const MResidue*> residues;

  foreach (const MChain* chain, protein.GetChains())
  {
    foreach (const MResidue* residue, chain->GetResidues())
      residues.push_back(residue);
  }

  // keep residues sorted by residue number as assigned during reading the PDB file
  sort(residues.begin(), residues.end(), boost::bind(&MResidue::GetNumber, _1) < boost::bind(&MResidue::GetNumber, _2));

  const MResidue* last = nullptr;
  foreach (const MResidue* residue, residues)
  {
    // insert a break line whenever we detect missing residues
    // can be the transition to a different chain, or missing residues in the current chain
    if (last != nullptr and last->GetNumber() + 1 != residue->GetNumber())
    {
      char breaktype = ' ';
      if (last->GetChainID() != residue->GetChainID())
        breaktype = '*';
      os << (kDSSPResidueLine % FixedLengthString(last->GetNumber() + 1, 5) % breaktype) << std::endl;
    }
    os << ResidueToDSSPLine(*residue) << std::endl;
    last = residue;
  }
}
