#include "timeSelector.H"
#include "argList.H"
#include "fvCFD.H"
#include "IFstream.H"

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
  
  // select a subset based on the command-line options
  instantList timeDirs = timeSelector::select
    (
     runTime.times(),
     args
     );
  
  if (timeDirs.empty())
    {
      FatalErrorIn(args.executable())
	<< "No times selected"
	<< exit(FatalError);
    }
  
  // Loop over all times
  forAll (timeDirs, timeI)
    {
      // Set time for global database
      runTime.setTime(timeDirs[timeI], timeI);
      
      Info << "Time = " << runTime.timeName() << endl << endl;

      // Get list of objects from processor0 database
      fileNameList objNames = readDir(runTime.timePath(), fileName::FILE);
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
	  fileName systemFile = systemNames[sysI];
	  Info << "Reconstruct object " << systemNames[sysI] << endl;

	  // IOdictionary sysDict
	  //   (
	  //    IOobject
	  //    (
	  //     systemNames[sysI],
	  //     runTime.timeName(),
	  //     runTime,
	  //     IOobject::MUST_READ,
	  //     IOobject::NO_WRITE
	  //     )
	  //    );
	  
	  Foam::IFstream infile(runTime.timePath()/systemNames[sysI]);
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
	  
	  // Write the vector
	  {
	    Foam::fileName fileName = runTime.timePath() / systemFile.lessExt() + ".vec";
	    Foam::OFstream file(fileName);
	    for (Foam::scalarField::const_iterator ptr = source.begin();
		 ptr != source.end();
		 ++ptr) {
	      file << *ptr << Foam::endl;
	    }
	  }
	  
	  // Write the matrix as "row col val" triples
	  {
	    
	    Foam::fileName fileName = runTime.timePath() / systemFile.lessExt() + ".coo";
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
      
    }
}
