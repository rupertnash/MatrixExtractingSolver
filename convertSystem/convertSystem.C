#include "timeSelector.H"
#include "argList.H"
#include "fvCFD.H"
#include "IFstream.H"
#include <stdint.h>
#include <limits>

void KeyError(const Foam::fileName& file, const Foam::word& key)
{
  Foam::FatalError("main(int, char* [])", __FILE__, __LINE__)
	<< "System file '" << file << "' has no keyword '" << key << "'"
	<< Foam::exit(Foam::FatalError);
}

/**
 * Traits class for writing NumPy format headers
 */
template<class T>
struct numpy_traits
{
  /**
   * Return the endianness character code for this platform.
   */
  static char endianness()
  {
    // Set up the int "1"
    int32_t number = 0x1;
    char *numPtr = reinterpret_cast<char*>(&number);
    // See if the first byte is 1, if so, it's little.
    if (numPtr[0] == 1)
      return '<';
    return '>';
  }
  /**
   * Return the NumPy type character for the argument type.
   */
  static char type()
  {
    if (std::numeric_limits<T>::is_integer)
      {
	if (std::numeric_limits<T>::is_signed)
	  return 'i';
	return 'u';
      }
    return 'f';
  }
  
  /**
   * Return the size in bytes of the type.
   */
  static char size()
  {
    switch(sizeof(T))
      {
      case 1:
	return '1';
      case 2:
	return '2';
      case 4:
	return '4';
      case 8:
	return '8';
      }
  }
  
  /**
   * Return the full, quoted datatype string
   */
  static std::string dtype()
  {
    std::string ans(5, '\0');
    ans[0] = '\'';
    ans[1] = endianness();
    ans[2] = type();
    ans[3] = size();
    ans[4] = '\'';
    return ans;
  }
};

/**
 * Construct the full numpy header, including padding.
 * See:
 *   https://github.com/numpy/numpy/blob/master/numpy/lib/format.py
 */
std::string MakeNumpyHeader(const std::string& dtypeExpression, int n)
{
  std::ostringstream headerStream;
  // Magic number
  headerStream << "\x93\x4E\x55\x4D\x50\x59"
    // Version
	       << std::string("\x01\x00", 2)
    // Dummy HEADER_LEN
	       << std::string("\x00\x00", 2)
    // Dict describing the data
	       << "{'descr': " << dtypeExpression << ", "
	       << "'fortran_order': False, "
	       << "'shape': (" << n << ",), }";
  
  // +1 for the required terminating newline		
  while ((headerStream.str().size() + 1) % 16)
    {
      headerStream << " ";
    }
  headerStream << "\n";
  uint16_t HEADER_LEN = headerStream.str().size() - 10; // 6 for magic, 2 for version, 2 for HEADER_LEN
  std::string header = headerStream.str();
  *(reinterpret_cast<uint16_t*>(&header[8])) = HEADER_LEN;
  return header;
}

// Typedefs to save a lot of typing!
typedef Foam::labelList::value_type idx;
typedef Foam::scalarField::value_type flt;

/**
 * Write a vector in platform native binary format.
 */
void WriteVecBinary(const Foam::scalarField& source, std::ostream& file)
{
  const char* charPtr = reinterpret_cast<const char*>(source.begin());
  std::streamsize sizeBytes = source.size() * sizeof(flt);
  file.write(charPtr, sizeBytes);
}

// Size of a single record for the COO files
const size_t recordLen = 2 * sizeof(idx) + sizeof(flt);
/**
 * Write the COO data in platform native binary format.
 */
void WriteCooBinary(const Foam::scalarField& diag,
		    const Foam::scalarField& upper,
		    const Foam::scalarField& lower,
		    const Foam::labelList& upperAddr,
		    const Foam::labelList& lowerAddr,
		    std::ostream& file)
{
  // Diagonal first
  char* diagBuf = new char[diag.size() * recordLen];
  for (idx i = 0; i < diag.size(); ++i)
    {
      // idx = (row, col)
      idx* idxPtr = reinterpret_cast<idx*>(diagBuf + (recordLen * i));
      // Row
      idxPtr[0] = i;
      // Col
      idxPtr[1] = i;
      // Val
      flt* valPtr = reinterpret_cast<flt*>(idxPtr + 2);
      *valPtr = diag[i];
    }
  file.write(diagBuf, diag.size() * recordLen);
  delete[] diagBuf;
  
  bool hasUpper = upper.size();
  bool hasLower = lower.size();
  
  // Only do more output if non-diagonal
  if (hasUpper || hasLower)
    {
      char* offDiagBuf = new char[2 * upperAddr.size() * recordLen];
      
      Foam::scalarField::const_iterator upPtr, upEnd;
      if (hasUpper)
	{
	  upPtr = upper.cbegin();
	  upEnd = upper.cend();
	}
      else
	{
	  upPtr = lower.cbegin();
	  upEnd = lower.cend();
	}
      Foam::scalarField::const_iterator loPtr, loEnd;
      if (hasLower)
	{
	  loPtr = lower.cbegin();
	  loEnd = lower.cend();
	}
      else
	{
	  loPtr = upper.cbegin();
	  loEnd = upper.cend();
	}
      
      // Get the indexing arrays
      Foam::labelList::const_iterator upAddrPtr = upperAddr.begin();
      Foam::labelList::const_iterator loAddrPtr = lowerAddr.begin();
      
      // Send to the buffer
      for(size_t i = 0; loPtr != loEnd; ++loPtr, ++upPtr, ++upAddrPtr, ++loAddrPtr)
	{
	  idx* idxPtr = reinterpret_cast<idx*>(offDiagBuf + (recordLen * i));
	  idxPtr[0] = *upAddrPtr;
	  idxPtr[1] = *loAddrPtr;
	  // Val
	  flt* valPtr = reinterpret_cast<flt*>(idxPtr + 2);
	  *valPtr = *loPtr;
	  ++i;
	  
	  idxPtr = reinterpret_cast<idx*>(offDiagBuf + (recordLen * i));
	  idxPtr[0] = *loAddrPtr;
	  idxPtr[1] = *upAddrPtr;
	  // Val
	  valPtr = reinterpret_cast<flt*>(idxPtr + 2);
	  *valPtr = *upPtr;
	  ++i;
	}
      file.write(offDiagBuf, 2 * upperAddr.size() * recordLen);
      delete[] offDiagBuf;
    }
}

int main(int argc, char* argv[])
{
  // enable -constant ... if someone really wants it
  // enable -zeroTime to prevent accidentally trashing the initial fields
  Foam::timeSelector::addOptions(false, false);
  Foam::argList::noParallel();
  Foam::argList::validOptions.set("format", "outputFormat");

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
  
  const string& format = args.optionFound("format") ? args.option("format") : "ascii";

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
	  Info << "Convert object " << systemNames[sysI] << endl;
	  
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
	    if (format == "ascii")
	      {
		std::ofstream file(fileName.c_str(), std::ios_base::out|std::ios_base::trunc);
		for (Foam::scalarField::const_iterator ptr = source.begin();
		     ptr != source.end();
		     ++ptr) {
		  file << *ptr << std::endl;
		}
	      }
	    else if (format == "binary")
	      {
		std::ofstream file(fileName.c_str(), std::ios_base::trunc | std::ios_base::binary);
		WriteVecBinary(source, file);
	      }
	    else if (format == "numpy")
	      {
		std::ofstream file(fileName.c_str(), std::ios_base::trunc | std::ios_base::binary);
		std::string header = MakeNumpyHeader(numpy_traits<flt>::dtype(), source.size());
		file.write(header.c_str(), header.size());
		WriteVecBinary(source, file);
	      }
	    else
	      {
		// ERROR
		FatalErrorIn(args.executable())
		  << "Unknown format: " << format
		  << exit(FatalError);
	      }
	  }
	  
	  // Write the matrix as "row col val" triples
	  {
	    Foam::fileName fileName = runTime.timePath() / systemFile.lessExt() + ".coo";

	    if (format == "ascii")
	      {
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
	      }
	    else if (format == "binary")
	      {
		std::ofstream file(fileName.c_str(), std::ios_base::trunc | std::ios_base::binary);
		WriteCooBinary(diag, upper, lower, upperAddr, lowerAddr, file);
	      }
	    else if (format == "numpy")
	      {
		size_t nRecords = diag.size();
		if (hasUpper || hasLower)
		  nRecords += 2 * upperAddr.size();
		
		std::ostringstream dtype;
		
		dtype << "["
		      << "('row', " << numpy_traits<idx>::dtype() << "), "
		      << "('col', " << numpy_traits<idx>::dtype() << "), "
		      << "('val', " << numpy_traits<flt>::dtype() << ")"
		      << "]";
		
		std::string header = MakeNumpyHeader(dtype.str(), nRecords);
		std::ofstream file(fileName.c_str(), std::ios_base::trunc | std::ios_base::binary);
		file.write(header.c_str(), header.size());
		
		WriteCooBinary(diag, upper, lower, upperAddr, lowerAddr, file);
	      }
	    else
	      {
		// ERROR
		FatalErrorIn(args.executable())
		  << "Unknown format: " << format
		  << exit(FatalError);
	      }
	    
	  } // End write COO
	  
	} // End for args
      
    }
}
