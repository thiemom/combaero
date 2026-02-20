#include "correlation_status.h"
#include <iostream>

namespace combaero {

namespace {

void default_warning_handler(const std::string& msg) {
    std::cerr << "[combaero] " << msg << "\n";
}

WarningHandler& global_handler() {
    static WarningHandler handler = default_warning_handler;
    return handler;
}

}  // namespace

void set_warning_handler(WarningHandler fn) {
    global_handler() = fn ? std::move(fn) : default_warning_handler;
}

WarningHandler get_warning_handler() {
    return global_handler();
}

void warn(const std::string& msg) {
    global_handler()(msg);
}

}  // namespace combaero
