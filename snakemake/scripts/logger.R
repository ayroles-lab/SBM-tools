library(log4r)

my_console_appender = console_appender(layout = default_log_layout())
my_file_appender = file_appender(my_logfile, append = TRUE, 
                            layout = default_log_layout())

my_logger <- log4r::logger(threshold = "INFO", 
                appenders= list(my_console_appender,my_file_appender))

log4r_info <- function(msg) {
  log4r::info(my_logger, msg)
}

log4r_error <- function(error) {
  log4r::error(my_logger, error)
}

log4r_debug <- function(debug) {
  log4r::debug(my_logger, debug)
}