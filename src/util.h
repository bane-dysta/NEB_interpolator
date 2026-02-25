#ifndef UTIL_H
#define UTIL_H

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <string>
#include <sstream>
#include <cstdlib>
#include <sys/wait.h>

namespace util {

// Quote a string for safe use as a single shell argument (POSIX sh).
// Uses single-quote wrapping and escapes embedded single quotes.
static inline std::string shellQuote(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('\'');
    for (char c : s) {
        if (c == '\'') {
            // close, escape, reopen: ' -> '\''
            out.append("'\\''");
        } else {
            out.push_back(c);
        }
    }
    out.push_back('\'');
    return out;
}

struct CommandResult {
    int raw_status = -1;      // raw return from system(); -1 means system() failed
    bool exited = false;      // true if process exited normally
    int exit_code = -1;       // valid when exited==true
    bool signaled = false;    // true if terminated by signal
    int term_signal = -1;     // valid when signaled==true
};

static inline CommandResult runCommand(const std::string& cmd) {
    CommandResult r;
    r.raw_status = std::system(cmd.c_str());

    if (r.raw_status == -1) {
        return r;
    }

    if (WIFEXITED(r.raw_status)) {
        r.exited = true;
        r.exit_code = WEXITSTATUS(r.raw_status);
    }

    if (WIFSIGNALED(r.raw_status)) {
        r.signaled = true;
        r.term_signal = WTERMSIG(r.raw_status);
    }

    return r;
}

static inline std::string formatCommandFailure(const CommandResult& r) {
    std::ostringstream oss;
    if (r.raw_status == -1) {
        oss << "failed to invoke system()";
    } else if (r.exited) {
        oss << "exit code " << r.exit_code;
    } else if (r.signaled) {
        oss << "terminated by signal " << r.term_signal;
    } else {
        oss << "non-zero status " << r.raw_status;
    }
    return oss.str();
}

// Create a unique temporary file path in the system temp directory.
//
// Notes:
// - This does not create the file; it only returns a path that is very unlikely
//   to collide.
// - If the system temp directory cannot be determined, it falls back to the
//   current working directory.
//
// Example:
//   auto p = util::makeTempFilePath("neb_final_", ".xyz");
//
static inline std::filesystem::path makeTempFilePath(const std::string& prefix,
                                                     const std::string& suffix) {
    std::error_code ec;
    std::filesystem::path dir = std::filesystem::temp_directory_path(ec);
    if (ec) {
        dir = std::filesystem::path(".");
    }

    // Monotonic counter to reduce collision probability even if called multiple
    // times in the same clock tick.
    static std::atomic<std::uint64_t> counter{0};

    for (int attempt = 0; attempt < 128; ++attempt) {
        const std::uint64_t c = counter.fetch_add(1, std::memory_order_relaxed);
        const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();

        std::filesystem::path p = dir / (prefix + std::to_string(static_cast<long long>(now)) + "_" +
                                         std::to_string(static_cast<unsigned long long>(c)) + suffix);

        ec.clear();
        if (!std::filesystem::exists(p, ec)) {
            return p;
        }
    }

    // Fallback: still return something sensible.
    return dir / (prefix + "fallback" + suffix);
}

} // namespace util

#endif // UTIL_H
