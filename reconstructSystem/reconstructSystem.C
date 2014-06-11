#include "timeSelector.H"
#include "argList.H"

#include "processorMeshes.H"

#include "fvCFD.H"
#include "IOobjectList.H"

int main(int argc, char* argv[])
{
  // enable -constant ... if someone really wants it
  // enable -zeroTime to prevent accidentally trashing the initial fields
  Foam::timeSelector::addOptions(false, false);
  Foam::argList::noParallel();
  
#   include "setRootCase.H"
#   include "createTime.H"
  
  label nProcs = 0;
  // Determine the processor count directly
  while (isDir(args.path()/(word("processor") + name(nProcs))))
    {
      ++nProcs;
    }
  
  if (!nProcs)
    {
      FatalErrorIn(args.executable())
	<< "No processor* directories found"
	<< exit(FatalError);
    }
  
  // Create the processor databases
  PtrList<Time> databases(nProcs);
  
  forAll (databases, procI)
    {
      databases.set
	(
	 procI,
	 new Time
	 (
	  Time::controlDictName,
	  args.rootPath(),
	  args.caseName()/fileName(word("processor") + name(procI))
	  )
	 );
    }
  
  // use the times list from the master processor
  // and select a subset based on the command-line options
  instantList timeDirs = timeSelector::select
    (
     databases[0].times(),
     args
     );
  
  if (timeDirs.empty())
    {
      FatalErrorIn(args.executable())
	<< "No times selected"
	<< exit(FatalError);
    }
  // These 2 lines may not be needed?
#   include "createNamedMesh.H"
  fileName regionPrefix = "";
  // Set all times on processor meshes equal to reconstructed mesh
  forAll (databases, procI)
    {
      databases[procI].setTime(runTime.timeName(), runTime.timeIndex());
    }
  
  // Read all meshes and addressing to reconstructed mesh
  processorMeshes procMeshes(databases, fvMesh::defaultRegion);
  
  // Loop over all times
  forAll (timeDirs, timeI)
    {
      // Set time for global database
      runTime.setTime(timeDirs[timeI], timeI);
      
      Info << "Time = " << runTime.timeName() << endl << endl;

      // Set time for all databases
      forAll (databases, procI)
        {
	  databases[procI].setTime(timeDirs[timeI], timeI);
        }
      
      // Check if any new meshes need to be read.
      fvMesh::readUpdateState meshStat = mesh.readUpdate();
      
      fvMesh::readUpdateState procStat = procMeshes.readUpdate();
      
      if (procStat == fvMesh::POINTS_MOVED)
        {
	  // Reconstruct the points for moving mesh cases and write them out
	  procMeshes.reconstructPoints(mesh);
        }
      else if (meshStat != procStat)
        {
	  WarningIn(args.executable())
	    << "readUpdate for the reconstructed mesh:" << meshStat << nl
	    << "readUpdate for the processor meshes  :" << procStat << nl
	    << "These should be equal or your addressing"
	    << " might be incorrect."
	    << " Please check your time directories for any "
	    << "mesh directories." << endl;
        }
      
      // Get list of objects from processor0 database
      fileNameList objNames = readDir(databases[0].timePath(), fileName::FILE);
      fileNameList systemNames;
      
      for (label objI = 0, sysI = 0; objI < objNames.size(); ++objI)
	{
	  if (objNames[objI].ext() == "system")
	    {
	      systemNames.setSize(systemNames.size() + 1);
	      systemNames[sysI] = objNames[objI];
	      sysI++;
	    }
	}

      forAll (systemNames, sysI)
	{
	  Info << "Reconstruct object " << systemNames[sysI] << endl;
	  
	  PtrList<scalarField> sources(nProcs);
	  
	  PtrList<scalarField> diags(nProcs);
	  PtrList<scalarField> uppers(nProcs);
	  PtrList<scalarField> lowers(nProcs);
	  
	  PtrList<labelList> upperAddrs(nProcs);
	  PtrList<labelList> lowerAddrs(nProcs);
	  
	  List<PtrList<labelList> > localCellIds(nProcs);
	  List<PtrList<scalarField> > coeffs(nProcs);
	  
	  label totalPoints = 0,
	    totalOffDiag = 0;
	  bool anyUpper = false,
	    anyLower = false;
	  
	  for (label procI = 0; procI < nProcs; ++procI)
	    {
	      IOdictionary sysDict
		(
		 IOobject
		 (
		  systemNames[sysI],
		  databases[procI].timeName(),
		  databases[procI],
		  IOobject::MUST_READ,
		  IOobject::NO_WRITE
		  )
		 );
	      
	      // Source always present
	      sources.set(procI, new scalarField());
	      if (!sysDict.readIfPresent("source", sources[procI]))
		Info << "Source not present!" << endl;
	      
	      dictionary& matrixDict = sysDict.subDict("matrix");
	      
	      // Diag always present
	      diags.set(procI, new scalarField());
	      if (!matrixDict.readIfPresent("diag", diags[procI]))
		Info << "Diag not present!" << endl;
	      
	      // Both upper and lower are optional
	      // Also grab the number of off-diag elems
	      label nOffDiag = 0;
	      bool hasUpper = matrixDict.found("upper");
	      if (hasUpper)
		{
		  uppers.set(procI, new scalarField());
		  matrixDict.readIfPresent("upper", uppers[procI]);
		  nOffDiag = uppers[procI].size();
		  anyUpper = true;
		}
	      bool hasLower = matrixDict.found("lower");
	      if (hasLower)
		{
		  lowers.set(procI, new scalarField());
		  matrixDict.readIfPresent("lower", lowers[procI]);
		  nOffDiag = lowers[procI].size();
		  anyLower = true;
		}
	      
	      // Non-diag => must have addressing arrays
	      if (hasUpper || hasLower)
		{
		  upperAddrs.set(procI, new labelList());
		  lowerAddrs.set(procI, new labelList());
		  if (!matrixDict.readIfPresent("upperAddr", upperAddrs[procI]))
		    FatalErrorIn(args.executable())
		      << "Missing upperAddr in non-diagonal matrix"
		      << exit(FatalError);
		  if (!matrixDict.readIfPresent("lowerAddr", lowerAddrs[procI]))
		    FatalErrorIn(args.executable())
		      << "Missing lowerAddr in non-diagonal matrix"
		      << exit(FatalError);
		}
	      
	      // Get the interprocess boundary data
	      dictionary& interfaceDict = sysDict.subDict("interfaces");
	      // Resize the ptrLists
	      localCellIds[procI].resize(nProcs);
	      coeffs[procI].resize(nProcs);
	      
	      // Loop over the possible processors (as I can't figure out how to iterate over a dictionary)
	      for (label procJ = 0; procJ < nProcs; ++procJ)
		{
		  if (procJ == procI)
		    continue;
		  
		  const word key = word("processor") + name(procJ);
		  
		  if (interfaceDict.found(key))
		    {
		      dictionary& procDict = interfaceDict.subDict(key);
		      localCellIds[procI].set(procJ, new labelList());
		      if (!procDict.readIfPresent("localCellIds", localCellIds[procI][procJ]))
			FatalErrorIn(args.executable())
			  << "Missing localCellIds in system for processor " << procI
			  << " with its boundary with processor " << procJ
			  << exit(FatalError);

		      coeffs[procI].set(procJ, new scalarField());
		      if (!procDict.readIfPresent("coeffs", coeffs[procI][procJ]))
			FatalErrorIn(args.executable())
			  << "Missing coeffs in system for processor " << procI
			  << " with its boundary with processor " << procJ
			  << exit(FatalError);

		      if (coeffs[procI][procJ].size() != localCellIds[procI][procJ].size())
			FatalErrorIn(args.executable())
			  << "Size mismatch in inteface data in system for processor " << procI
			  << " with its boundary with processor " << procJ
			  << exit(FatalError);
		      
		      // Protect against double-counting of off diagonal elements.
		      if (procI > procJ)
			nOffDiag += coeffs[procI][procJ].size();
		    }
		}
	      
	      totalPoints += sources[procI].size();
	      totalOffDiag += nOffDiag;
	    }
	  
	  scalarField totalSource(totalPoints);
	  scalarField totalDiag(totalPoints);
	  
	  label totalUpperEntries = anyUpper ? totalOffDiag : 0;
	  label totalLowerEntries = anyLower ? totalOffDiag : 0;
	  scalarField totalUpper(totalUpperEntries);
	  scalarField totalLower(totalLowerEntries);
	  labelList totalUpperAddr(totalOffDiag);
	  labelList totalLowerAddr(totalOffDiag);
	  label offDiagI = 0;

	  for (label procI = 0; procI < nProcs; ++procI)
	    {
	      const labelIOList& addressesProcI = procMeshes.cellProcAddressing()[procI];
	      
	      // To deal with deleting the loaded data aggressively, we're going to
	      // NULL out the entry grab an autoPtr to the item. The item will then
	      // be freed once the autoPtr goes out of scope at the end of the loop
	      // iteration.
	      autoPtr<scalarField> src = sources.set(procI, NULL);
	      autoPtr<scalarField> dia = diags.set(procI, NULL);
	      
	      if (addressesProcI.size() != src->size())
		FatalErrorIn(args.executable())
		  << "Size mismatch"
		  << exit(FatalError);
	      
	      // Copy over the source and diagonal, translating the
	      // addressesProcI from proc-local to global.
	      for (label i = 0; i < src->size(); ++i)
		{
		  totalSource[addressesProcI[i]] = src()[i];
		  totalDiag[addressesProcI[i]] = dia()[i];
		}
	      
	      // Now deal with the processor internal off-diagonal elements.
	      // Upper, if it's there
	      autoPtr<scalarField> upp = uppers.set(procI, NULL);
	      bool hasUpper = upp.valid();
	      if (hasUpper)
		{
		  memcpy(totalUpper.data() + offDiagI,
			 upp->cdata(),
			 upp->byteSize());
		}
	      
	      // Lower, if it's there
	      autoPtr<scalarField> low = lowers.set(procI, NULL);
	      bool hasLower = low.valid();
	      if (hasLower)
		{
		  memcpy(totalLower.data() + offDiagI,
			 low->data(),
			 low->byteSize());
		}
	      
	      // Addressing if either upper or lower is present
	      autoPtr<labelList> ua = upperAddrs.set(procI, NULL);
	      autoPtr<labelList> la = lowerAddrs.set(procI, NULL);
	      if (hasUpper || hasLower)
		{
		  // Grab references to give the compiler a chance.
		  const labelList& uar = ua();
		  const labelList& lar = la();
		  // Copy the addressing over, translating from local
		  // to global. Also note that we increment both i and
		  // offDiagI.
		  for (label i = 0; i < ua->size(); ++i, ++offDiagI)
		    {
		      totalUpperAddr[offDiagI] = addressesProcI[uar[i]];
		      totalLowerAddr[offDiagI] = addressesProcI[lar[i]];
		    }
		}
	      
	      // Now onto the inter-process off-diagonal elements. We
	      // loop over all procJ > procI to ensure we don't double
	      // count. (We have to deal with each pairing exactly
	      // once.)
	      for (label procJ = procI + 1; procJ < nProcs; ++ procJ)
		{
		  autoPtr<labelList> cellIdsIJ = localCellIds[procI].set(procJ, NULL);
		  autoPtr<labelList> cellIdsJI = localCellIds[procJ].set(procI, NULL);
		  
		  autoPtr<scalarField> coeffsIJ = coeffs[procI].set(procJ, NULL);
		  autoPtr<scalarField> coeffsJI = coeffs[procJ].set(procI, NULL);
		  
		  if (cellIdsIJ.valid())
		    {
		      // Transpose must exist
		      if (!cellIdsJI.valid())
			FatalErrorIn(args.executable())
			  << "Missing expected boundary data (due to symmetry) from " 
			  << procJ << " to " << procI
			  << exit(FatalError);
		      
		      // Sizes of corresponding patches must match
		      if (cellIdsIJ->size() != cellIdsJI->size())
			FatalErrorIn(args.executable())
			  << "Patch size mismatch between processors " 
			  << procI << " and " << procJ
			  << exit(FatalError);
		      
		      const labelIOList& addressesProcJ = procMeshes.cellProcAddressing()[procJ];
		      // Get refs to help the compiler...
		      const scalarField& cIJ = coeffsIJ();
		      const scalarField& cJI = coeffsJI();
		      const labelList& idIJ = cellIdsIJ();
		      const labelList& idJI = cellIdsJI();
		      
		      
		      for (label i = 0; i < cellIdsIJ->size(); ++i, ++offDiagI)
			{
			  // Recall that we have to negate the coefficient as we copy them in.
			  // The code in processorFvPatchField.C makes clear which way round these should go.
			  if (anyUpper)
			    totalUpper[offDiagI] = -cIJ[i];
			  if (anyLower)
			    totalLower[offDiagI] = -cJI[i];
			  
			  totalUpperAddr[offDiagI] = addressesProcI[idIJ[i]];
			  totalLowerAddr[offDiagI] = addressesProcJ[idJI[i]];
			}
		    }
		}
	      
	    } // for (procI... )

	  // Prepare a dict for writing the reconstructed one back.
	  IOdictionary systemDict(
				  IOobject
				  (
				   systemNames[sysI],
				   runTime.timeName(),
				   runTime,
				   IOobject::NO_READ,
				   IOobject::NO_WRITE
				   )
				  );
	  systemDict.add("source", totalSource);
	  
	  dictionary totalMatrix;
	  totalMatrix.add("diag", totalDiag);
	  if (anyUpper)
	    totalMatrix.add("upper", totalUpper);
	  if (anyLower)
	    totalMatrix.add("lower", totalLower);
	  if (anyUpper || anyLower)
	    {
	      totalMatrix.add("upperAddr", totalUpperAddr);
	      totalMatrix.add("lowerAddr", totalLowerAddr);
	    }
	  systemDict.add("matrix", totalMatrix);

	  // Actually write!
	  systemDict.regIOobject::write();
	}
    }
}

  

