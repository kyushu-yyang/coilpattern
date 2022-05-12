#ifndef XLogger_HH
#define XLogger_HH

#include <iomanip>
#include <string>
#include <fstream>

using std::string;

/// @file   XLogger.h
/// @author Ye Yang (QST)
/// @date   05.12.2022

class XLogger
{
  public:
    /// @enum Level
    enum Level { DEBUG, CONFIG, INFO, WARNING, ERROR };

    /// @brief retutn the instance
    static XLogger* GetInstance();

    /// @brief start logging and write message to file
    static void Start(Level minpriority, const string& filename);

    /// @brief stop logging
    static void Stop();

    /// @brief write a message
    static void Write(Level priority, const string& message);

    /// @brief return the logging stream
    static std::ostream& GetLogStream(Level priority);

  private:
    /// @brief constructor
    XLogger();

    /// @brief deconstructor
    XLogger(const XLogger& log) {}

    /// @brief static instance of this class
    static XLogger* instance;

    /// @brief priority names
    static const string GetPriorityName(Level level);

  private:
    bool          fActive;
    std::ofstream fFile;
    Level         fPriority;
};

/// @def macro to handle the output of error message
#ifndef _OUTPUT_ERROR
#define _OUTPUT_ERROR(level, message)                \
  do {                                               \
    XLogger* log = XLogger::GetInstance();           \
    log->GetLogStream(level) << message << " "       \
      << std::setprecision(6) << std::setw(6)        \
      << std::setfill(' ') << std::endl;             \
  } while(0)
#endif

#ifndef ERROR_OUTPUT
  #define ERROR_OUTPUT true
#endif

#ifndef Error
#define Error(level, message)                        \
  do {                                               \
    if (ERROR_OUTPUT){                               \
      _OUTPUT_ERROR(level, message);                 \
    }} while(0)
#endif

#ifndef Debug
#define Debug(message)                               \
  do {                                               \
    Error(XLogger::DEBUG, message);                  \
  } while(0)
#endif

#ifndef Fatal
#define Fatal(message)                               \
  do {                                               \
    Error(XLogger::ERROR, message);                  \
  } while(0)
#endif

#ifndef Warning
#define Warning(message)                             \
  do {                                               \
    Error(XLogger::WARNING, message);                \
  } while(0)
#endif

#ifndef Info
#define Info(message)                                \
  do {                                               \
    Error(XLogger::INFO, message);                   \
  } while(0)
#endif
  
#endif
