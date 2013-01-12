#include "../includes/logservice.hpp"

using namespace ICoDF;

LogService* ICoDF::LogService::_singleton = NULL; // initialize the singleton pointer to NULL

/// ADDMESSAGE
/// Send a new message to the log service.
/// Message will be display, write or saved depending on configuration 
void ICoDF::LogService::AddMessage(short int msgType, std::string module, std::string message)
{
    std::stringstream msg;
    time_t now = time(0);
    struct tm *ts = localtime(&now);
    char date[80];
    strftime(date, sizeof(date), "%c", ts);

    // Create the message
    if (msgType & LogService::NOTICE)
    {
        msg << MSG_SQUELLETON("(o)");
    }
    else if (msgType & LogService::WARNING)
    {
        msg << MSG_SQUELLETON("/!\\");
    }
    else if (msgType & LogService::FATAL)
    {
        msg << MSG_SQUELLETON("[X]");
    }

    // Handle it according to the configuration
    if (this->_config & LS_WRITE_TO_FILE)
        // Create and check the file
        if (this->CheckFile())
            this->_logFile << msg.str() << std::endl;
    if (this->_config & LS_PRINT_ON_COUT)
        // Display message on cout
        std::cout << msg.str() << std::endl;
}

// SETCONFIGURATION
// Set a new configuration for log management
void ICoDF::LogService::SetConfiguration(short int config)
{
    this->_config = config;
}

// TODO : Secure file loading (try,catch)
bool ICoDF::LogService::CheckFile()
{
    if (!this->_logFile.good() | this->_fileName.empty())
    {
        if (this->_config & LS_WRITE_TO_FILE)
        {
            if (this->_fileName.empty())
            {
                // Create the file name
                std::ostringstream fileName;
                time_t now = time(0);
                struct tm *ts = localtime(&now);
                char buf[80];
                strftime(buf, sizeof(buf), "%c", ts);
                fileName << "BLINK_LOG_[" << buf << "].txt";
                this->_fileName = fileName.str();
            }
            // Open the file
            this->_logFile.open(this->_fileName.c_str(), std::fstream::app | std::fstream::out | std::fstream::ate);
            // Check if the stream has been correctly open
            if (!this->_logFile.good())
            {
                std::cerr << "(WARNING) BLINK : Unable to access logFile [" << this->_fileName << "]" << std::endl;
                return false;
            }
        }
    }
    return true;
}

/// GETINSTANCE
/// Create a singleton instance if applicable and/or return the pointer
LogService* ICoDF::LogService::GetInstance()
{
    if (LogService::_singleton == NULL)
    {
        LogService::_singleton = new LogService;
    }
    return LogService::_singleton;
}

// DELETE
// Delete the singleton instance.
void ICoDF::LogService::Delete()
{
    delete LogService::_singleton;
}

// DEFAULT CTOR
ICoDF::LogService::LogService()
{
}

// DEFAULT DTOR
ICoDF::LogService::~LogService()
{
    this->AddMessage(LogService::NOTICE, "LogService", "shutting down log service...");
    if (this->_logFile.good())
    {
        this->_logFile.flush();
        this->_logFile.close();
    }
}
