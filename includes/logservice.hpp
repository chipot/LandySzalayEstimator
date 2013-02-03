#pragma once

#include <fstream>

// USE THIS COMMAND TO ADD A NEW MESSAGE
// X = log type (see static const (NOTICE|WARNING|FATAL) below)
// Y = module's name (string)
// Z = log message (string)
#define LS_ADDMSG(x,y,z) LogService::GetInstance()->AddMessage(x, y, z)
#define MSG_SQUELLETON(x) x << date << " : " << module << " => " << message;

// LogService manages all log system for the BLINK project.
// For now it is able to wirte logs on a file or on the term.
// Todo : also manages log using type and date (priority queue)
// Todo : Storing messages on dB ?
class LogService
{
 public:
  // Add a message to the log service
  void	AddMessage(short int type, std::string module, std::string message);

  /// Set a new configuration
  void	SetConfiguration(short int config);

 private:
  // Check if file is correctly setted.
  bool CheckFile(void);

 public:
  // Statics that defines the log type.
  // Log type determines the level of the message
  // In configuration allow to determine the log's types to print or write down.
  static const short int NOTICE = 1 << 0;
  static const short int WARNING = 1 << 1;
  static const short int FATAL = 1 << 2;

  // Statics defines that defines configuration options.
  static const short int LS_WRITE_NOTICE_TO_FILE = 1 << 0; // write notices to log file.
  static const short int LS_PRINT_NOTICE_ON_COUT = 1 << 1; // print notices to log on the display
  static const short int LS_WRITE_WARNING_TO_FILE = 1 << 2; // Guess what ?
  static const short int LS_PRINT_WARNING_ON_COUT = 1 << 3; // Obvious names isn't it !
  static const short int LS_WRITE_FATAL_TO_FILE = 1 << 4; // Let's keep it simple
  static const short int LS_PRINT_FATAL_ON_COUT = 1 << 5; // thanks for using it...
  static const short int LS_WRITE_TO_FILE = LS_WRITE_NOTICE_TO_FILE | LS_WRITE_WARNING_TO_FILE | LS_WRITE_FATAL_TO_FILE; // ...

  static const short int LS_PRINT_ON_COUT = LS_PRINT_NOTICE_ON_COUT | LS_PRINT_WARNING_ON_COUT | LS_PRINT_FATAL_ON_COUT; //

 private:
  short int	            _config; //< Configuration table

  std::fstream          _logFile; //< file stream

  std::string           _fileName; //< log file's name

 private:
  static LogService    *_singleton; //< singleton ptr for the service.

 public:
  /// Return a pointer the the log service instance.
  static LogService* GetInstance(void);

  /// Delete the singleton instance
  static void Delete(void);

 private:
  // DEFAULT CTOR
  LogService(void);

  // DEFAULT DTOR
  ~LogService(void);
};
