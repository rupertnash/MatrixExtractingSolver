#include "timeSelector.H"
#include "argList.H"

#include "processorMeshes.H"

#include "fvCFD.H"
#include "IOobjectList.H"

void KeyError(const Foam::fileName& file, const Foam::word& key)
{
  Foam::FatalError("main(int, char* [])", __FILE__, __LINE__)
	<< "System file '" << file << "' has no keyword '" << key << "'"
	<< Foam::exit(Foam::FatalError);
}

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
	      const List<label>* addresses = procMeshes.cellProcAddressing()(procI);
	      if (addresses->size() != sources(procI)->size())
		FatalErrorIn(args.executable())
		  << "Size mismatch"
		  << exit(FatalError);
	      for (label i = 0; i < sources(procI)->size(); ++i)
		{
		  totalSource[(*addresses)[i]] = sources[procI][i];
		  totalDiag[(*addresses)[i]] = diags[procI][i];
		}
	      
	      label nOff = 0;
	      bool hasUpper = (uppers(procI) != NULL);
	      if (hasUpper)
		nOff = uppers[procI].size();
	      bool hasLower = (lowers(procI) != NULL);
	      if (hasLower)
		nOff = lowers[procI].size();
	      
	      bool hasOffDiag = hasUpper || hasLower;
	      for (label i = 0; i < nOff; ++i)
		{
		  if (hasUpper)
		    totalUpper[offDiagI] = uppers[procI][i];
		  if (hasLower)
		    totalLower[offDiagI] = lowers[procI][i];
		  
		  totalLowerAddr[offDiagI] = (*addresses)[
							  lowerAddrs[procI][i]
							  ];
		  totalUpperAddr[offDiagI] = (*addresses)[
							  upperAddrs[procI][i]
							  ];
		  ++offDiagI;
		}
	    }

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

  /*  for (int i = 1; i < argc; ++i)
	{
	  Foam::fileName systemFile(argv[i]);
	  if (systemFile.ext() != "system")
	    {
	      Foam::Info << "Skipping file " << systemFile << Foam::endl;
	      continue;
	    }
	  
	  Foam::IFstream infile(systemFile);
	  const Foam::dictionary sysDict(infile);
	  
	  // Read the source vector
	  Foam::scalarField source;
	  if (!sysDict.readIfPresent("source", source, false, false))
	    KeyError(systemFile, "source");
	  
	  // Get the matrix sub dictionary
	  if (!sysDict.isDict("matrix"))
	    KeyError(systemFile, "matrix");
	  const Foam::dictionary& matDict = sysDict.subDict("matrix");
    
	  // Read the diagonal array - always present
	  Foam::scalarField diag;
	  if (!matDict.readIfPresent("diag", diag, false, false))
	    KeyError(systemFile, "diag");
    
	  // Upper and lower are optional
	  Foam::scalarField upper;
	  bool hasUpper = matDict.readIfPresent("upper", upper, false, false);
	  Foam::scalarField lower;
	  bool hasLower = matDict.readIfPresent("lower", lower, false, false);
    
	  // Non-diagonal matrices must have the addressing arrays.
	  Foam::labelList upperAddr;
	  Foam::labelList lowerAddr;
	  if (hasLower || hasUpper)
	    {
	      if (!matDict.readIfPresent("upperAddr", upperAddr, false, false))
		KeyError(systemFile, "upperAddr");
	      if (!matDict.readIfPresent("lowerAddr", lowerAddr, false, false))
		KeyError(systemFile, "lowerAddr");
	    }
	  
	  // Write the vector as "i val" pairs
	  {
	    Foam::fileName fileName = systemFile.lessExt() + ".vec";
	    Foam::OFstream file(fileName);
	    Foam::scalarField::const_iterator ptr = source.begin(),
	      end = source.end();
	    for (int i = 0; ptr != end; ++ptr, ++i) {
	      file << i << " " << *ptr << Foam::endl;
	    }
	  }
	  
	  // Write the matrix as "row col val" triples
	  {
	    
	    Foam::fileName fileName = systemFile.lessExt() + ".coo";
	    //Open the file.
	    Foam::OFstream file(fileName);
	    
	    file << "# row col val" << Foam::endl;
	    
	    // Diagonal first
	    Foam::scalarField::const_iterator ptr = diag.cbegin(), end = diag.cend();
	    for (int i = 0; ptr != end; ++ptr, ++i) {
	      file << i << " " << i << " " << *ptr << Foam::endl;
	    }
	    
	    // Only do more output if non-diagonal
	    if (hasUpper || hasLower)
	      {
		Foam::scalarField::const_iterator upPtr, upEnd;
		if (hasUpper)
		  {
		    upPtr = upper.cbegin();
		    upEnd = upper.cend();
		  } else {
		  upPtr = lower.cbegin();
		  upEnd = lower.cend();
		}
		Foam::scalarField::const_iterator loPtr, loEnd;
		if (hasLower)
		  {
		    loPtr = lower.cbegin();
		    loEnd = lower.cend();
		  } else {
		  loPtr = upper.cbegin();
		  loEnd = upper.cend();
		}
		
		// Get the indexing arrays
		Foam::labelList::const_iterator upAddrPtr = upperAddr.begin();
		Foam::labelList::const_iterator loAddrPtr = lowerAddr.begin();
		
		// Send to the file
		for(; loPtr != loEnd; ++loPtr, ++upPtr, ++upAddrPtr, ++loAddrPtr)
		  {
		    file << *upAddrPtr << " " << *loAddrPtr << " " << *loPtr << Foam::endl;
		    file << *loAddrPtr << " " << *upAddrPtr << " " << *upPtr << Foam::endl;
		  }
	      }
	  } // End write COO
	  
	} // End for args
  
	}*/
  
  

