export module utils.logger;

import std;

export enum class LogLevel {
    Debug = 0,
    Info  = 1,
    Warn  = 2,
    Error = 3
};

export class Logger {
public:
    static auto instance() -> Logger& {
        static Logger logger;
        return logger;
    }

    auto setLevel(LogLevel level) -> void { level_ = level; }
    auto level() const -> LogLevel { return level_; }
    auto isEnabled(LogLevel lv) const -> bool { return lv >= level_; }

    template<typename... Args>
    auto log(LogLevel lv, std::format_string<Args...> fmt, Args&&... args) const -> void {
        if (lv < level_) return;
        std::println(fmt, std::forward<Args>(args)...);
    }

private:
    LogLevel level_ = LogLevel::Info;
    Logger() = default;
};
