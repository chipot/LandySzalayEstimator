#pragma once
#ifndef LOG_VQWPMFGE
#define LOG_VQWPMFGE

#include <ctime>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>

namespace llog {

enum class level : int
{
    DEBUG,
    NOTICE,
    WARNING,
    FATAL,
};

struct level_base
{
    static std::string now()
    {
        std::time_t now = std::time(nullptr);
        std::tm *ts = std::localtime(&now);
        char date[80];

        strftime(date, sizeof(date), "%c", ts);
        return std::string(date);
    }
};

struct level_debug : private level_base
{
    static level const L = level::NOTICE;
    static std::string const name() {
        return "debug";
    }
    static std::string const prefix() {
        return "(&)" + now() + " : ";
    }
};

struct level_notice : private level_base
{
    static level const L = level::NOTICE;
    static std::string const name() {
        return "notice";
    }
    static std::string const prefix() {
        return "(o)" + now() + " : ";
    }
};

struct level_warning : private level_base
{
    static level const L = level::WARNING;
    static std::string const name() {
        return "warning";
    }
    static std::string const prefix() {
        return "/!\\" + now() + " : ";
    }
};

struct level_fatal : private level_base
{
    static level const L = level::FATAL;
    static std::string const name() {
        return "fatal";
    }
    static std::string const prefix() {
        return "[X]" + now() + " : ";
    }
};

template <class Level>
class notifier
{
 public:
    class notifier_forward
    {
     public:
        notifier_forward(notifier &n)
            : _notifier(n)
        {}
        notifier_forward(notifier_forward &&f)
            : _notifier(f._notifier)
        {}

        notifier_forward(notifier_forward const &f) = delete;

        ~notifier_forward()
        {
            this->_notifier._out << this->_cache.str();
        }

        template <typename T>
        notifier_forward const & operator << (T &&type) const
        {
            this->_cache << std::forward<T>(type);
            return *this;
        }
        notifier_forward const & operator << (std::ostream & (*func)(std::ostream &)) const
        {
            func(this->_cache);
            return *this;
        }
     private:
        notifier &_notifier;
        std::stringstream mutable _cache;
    };

    notifier()
        : notifier(std::clog)
    {}

    notifier(std::ostream &out)
        : _out(out)
    {}

    template <typename T>
    notifier_forward operator << (T &&type)
    {
        auto p = notifier_forward(*this);
        p << Level::prefix() << std::forward<T>(type);
        return p;
    }
    notifier_forward operator << (std::ostream & (*func)(std::ostream &))
    {
        auto p = notifier_forward(*this);
        p << Level::prefix() << func;
        return p;
    }

    notifier_forward operator [] (std::string const &str)
    {
        auto p = notifier_forward(*this);
        p << Level::prefix() << str << " => ";
        return p;
    }
 private:
    std::ostream &_out;
    friend class notifier_forward;
};

extern notifier<level_notice> notice;
extern notifier<level_warning> warn;
extern notifier<level_fatal> fatal;
extern notifier<level_debug> debug;

} /* log */


#endif /* end of include guard: LOG_VQWPMFGE */
