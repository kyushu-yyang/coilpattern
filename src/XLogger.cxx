#include <iostream>
#include "XLogger.h"

XLogger* XLogger::instance = NULL;

XLogger :: XLogger()
  : fActive(false)
{}

XLogger* XLogger :: GetInstance()
{
  if (instance==NULL)
    instance = new XLogger();

  return instance;
}

void XLogger :: Start(Level minpriority, const string& filename)
{
  instance->fActive = true;
  instance->fPriority = minpriority;

  if ( filename != " " )
    instance->fFile.open( filename.c_str() );
}

void XLogger :: Stop()
{
  instance->fActive = false;

  if ( instance->fFile.is_open() )
    instance->fFile.close();
}

void XLogger :: Write(Level priority, const string& message)
{
  if ( instance->fActive && priority>=instance->fPriority ) {
    std::ostream& stream = instance->fFile.is_open() ? instance->fFile : std::cout;

    stream << "["
           << instance->GetPriorityName(priority)
           << "]"
           << ": "
           << message;
  }
}

std::ostream& XLogger :: GetLogStream(Level priority)
{
  if (!instance->fFile.is_open())
    return std::cout;

  instance->fFile << "["
                  << instance->GetPriorityName(priority)
                  << "]"
                  << ": ";

  return instance->fFile;
}

const string XLogger :: GetPriorityName(Level level)
{
  const std::string priority[5] = {
      "DEBUG",
      "CONFIG",
      "INFO",
      "WARNING",
      "ERROR"
  };

  switch (level) {
    case   DEBUG: return priority[0]; break;
    case  CONFIG: return priority[1]; break;
    case    INFO: return priority[2]; break;
    case WARNING: return priority[3]; break;
    case   ERROR: return priority[4]; break;
    default: return ""; break;
  }
}

